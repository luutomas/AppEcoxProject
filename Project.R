#install.packages('devtools')
#install.packages("quantmod")
#install.packages("fGarch")
#install.packages("repr")
#install.packages("forecast")
#install.packages("aTSA")
#install.packages("rugarch") 
#install.packages("tseries")

library(quantmod)
library(fGarch)
library(repr)
library(forecast)
library(aTSA)
library(rugarch)
library(tseries)
library(MLmetrics)
library(stargazer)

bestARIMA <- function(dat,max_p,max_q, d, const){
  model_config <-c()
  model_AIC <-c()
  model_BIC <-c()
  Ljung_Box_pval <-c()
  RMSE <- c()
  for (i in 1:(max_p+1)){
    for (j in 1:(max_q+1)){
      arima<- tryCatch(Arima(dat,  # variable
                             order = c(i-1,d,j-1),  # (p,d,q) parameters
                             include.constant = const),
                       warning = function(w) {print(paste("non-finite finite-difference value", ''));
                         NaN},
                       error = function(e) {print(paste("non-finite finite-difference value", ''));
                         NaN})
      
      if (!is.numeric(arima)){
        model_config <-c(model_config, 
                         paste(i-1,d,j-1,sep=',')
        )
        model_AIC <- c(model_AIC, (AIC(arima)))
        model_BIC <- c(model_BIC, (BIC(arima)))
        Ljung_Box_pval <-c(Ljung_Box_pval,(as.numeric(sub(".*p-value = ","", 
                                                          capture.output(checkresiduals(arima, plot=FALSE)[5]))[5])))
        RMSE_tmp <-round(sqrt(mean((as.vector(arima$fitted) - as.vector(dat))^2)),1)
        RMSE <- c(RMSE, RMSE_tmp)
      }
    }
  }
  df <- as.data.frame(cbind(model_config, 
                            model_AIC,
                            as.numeric(as.character(model_BIC)),
                            as.numeric(as.character(Ljung_Box_pval)),
                            RMSE))
  names(df) <- c('model_config','model_AIC','model_BIC','Ljung_Box_pval','RMSE')
  df$model_config <- as.character(df$model_config)
  df$model_AIC <- as.numeric(as.character(df$model_AIC))
  df$model_BIC <- as.numeric(as.character(df$model_BIC))
  df$Ljung_Box_pval <- as.numeric(as.character(df$Ljung_Box_pval))
  return(df)
}


##### Getting the stock data for the chosen stock ####
getSymbols('BABA', src = 'yahoo', return.class = 'xts',from = "2014-09-20",to="2019-12-31")
head(BABA)
BABA <- BABA[,"BABA.Close"]
plot.xts(BABA, ylab = NA) #clearly non-stationary

BABA_logret <- na.omit(diff(log(BABA)))

plot.xts(BABA_logret, ylab = NA)
hist(BABA_logret,freq=FALSE,breaks=100)
curve(dnorm(x, mean=mean(BABA_logret), sd=sd(BABA_logret)), add=TRUE, col="red") #seems stationary

adf.test(BABA_logret)
Box.test(BABA_logret, type = "Ljung-Box")
kpss.test(BABA_logret)
adf <- data.frame("lags"=1:10,"p-val"=NA)
for (i in 1:10){
  adf[i,"p.val"] =as.numeric(adf.test(BABA_logret, k = i)$p.val)  
}
adf #all lags below confidence level, seems stationary

#examining ACF and PACF
par(mfrow = c(1, 2))
acf(BABA_logret)
pacf(BABA_logret)
par(mfrow = c(1,1))

Box.test(randomNormTS, type = "Ljung")
# testing joint significance of lags
LB <- data.frame("lags"=1:10,"p-val"=NA)
for (i in 1:10){
  LB[i,"p.val"] = Box.test(BABA_logret, type = "Ljung-Box", lag = i)$p.val
}
LB
# null hypothesis of autocorrelations up to lag k equal zero is very likely to be rjected as shown by LB text
# From ACF and PACF, LB test, it appears that there is a dependence on lags in both subsamples => proceed to GARCH
# Unfortunately from ACF and PACF, it is unclear which order of AR or MA should we choose
# To address this issue, we will use find the minimzed AIC and BIC

bestARIMA_results <- bestARIMA(BABA_logret, 10, 10, 0, FALSE)

head(bestARIMA_results[order(bestARIMA_results$model_AIC),],5)
head(bestARIMA_results[order(bestARIMA_results$model_BIC),],5)
# AIC (2,0,2), BIC (0,0,0)

auto.arima(BABA_logret,ic ="aic", stepwise = FALSE)
auto.arima(BABA_logret,ic ="bic", stepwise = FALSE)
# We also used out-of-the-box function, which finds the best order by minimizing information criterion
# AIC (2,0,3), BIC (0,0,0)

arima202 <- arima(BABA_logret, order = c(2,0,2))
checkresiduals(arima202)
arch.test(arima202)

arima000 <- arima(BABA_logret, order = c(0,0,0))
checkresiduals(arima000)
arch.test(arima000)
#Ljung Box test not sufficient therefore model rejected

# Homoskedastic residuals rejected => focus on conditional volatility => GARCH family models

#### Conditional Volatility ####
# Save last 183 observations (6 months) for out of sample forecasting
arma202_garch11_spec = ugarchspec(mean.model = list(armaOrder = c(2, 2)),
                                  variance.model = list(garchOrder = c(1, 1)))
arma202_garch11 = ugarchfit(spec = arma202_garch11_spec,data = residuals(arima202),
                            out.sample = 183,solver = "hybrid"
                          )

arma202_garch11_fc <- ugarchforecast(arma202_garch11, n.roll = 183, n.ahead = 1)
str(arma202_garch11_fc)

# From the results, we can see that all coefficients except mean are significant 
# and null hypothesis of no autocorrelation can not be rejected (LB test)

bestGARCH_AIC <- function(arima_model) {
  p.max <- 5
  q.max <- 5
  aic.min <- Inf
  best.p <- 0
  best.q <- 0
  inf.crit <- 0
  for (i1 in 1:p.max) {
    for (i2 in 1:q.max) {
      ourSpec <-
        ugarchspec(
          mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
          variance.model = list(garchOrder = c(i1, i2))
        )
      fit <- ugarchfit(spec = ourSpec, data = resid(arima_model), solver = "hybrid")
      inf.crit <- infocriteria(fit)[1]
      aic.min <- ifelse(inf.crit < aic.min, inf.crit, aic.min)
      
      best.p <- ifelse(inf.crit == aic.min, i1, best.p)
      best.q <- ifelse(inf.crit == aic.min, i2, best.q)
    }
  }
  return(c(best.p, best.q))
}
#print(bestGARCH_AIC(arima000))
#print(bestGARCH_AIC(arima202))

bestGARCH_BIC <- function(arima_model) {
  p.max <- 5
  q.max <- 5
  bic.min <- Inf
  best.p <- 0
  best.q <- 0
  inf.crit <- 0
  for (i1 in 1:p.max) {
    for (i2 in 1:q.max) {
      ourSpec <-
        ugarchspec(
          mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
          variance.model = list(garchOrder = c(i1, i2))
        )
      fit <- ugarchfit(spec = ourSpec, data = resid(arima_model), solver = "hybrid")
      inf.crit <- infocriteria(fit)[2]
      bic.min <- ifelse(inf.crit < bic.min, inf.crit, bic.min)
      
      best.p <- ifelse(inf.crit == bic.min, i1, best.p)
      best.q <- ifelse(inf.crit == bic.min, i2, best.q)
    }
  }
  return(c(best.p, best.q))
}
#print(bestGARCH_BIC(arima000))
#print(bestGARCH_BIC(arima202))

bestGJRGARCH_AIC <- function(arima_model) {
  p.max <- 5
  q.max <- 5
  aic.min <- Inf
  best.p <- 0
  best.q <- 0
  inf.crit <- 0
  for (i1 in 1:p.max) {
    for (i2 in 1:q.max) {
      ourSpec <-
        ugarchspec(
          mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
          variance.model = list(garchOrder = c(i1, i2), model = "gjrGARCH")
        )
      fit <- ugarchfit(spec = ourSpec, data = resid(arima_model), solver = "hybrid")
      inf.crit <- infocriteria(fit)[1]
      aic.min <- ifelse(inf.crit < aic.min, inf.crit, aic.min)
      
      best.p <- ifelse(inf.crit == aic.min, i1, best.p)
      best.q <- ifelse(inf.crit == aic.min, i2, best.q)
    }
  }
  return(c(best.p, best.q))
}
bestGJRGARCH_BIC <- function(arima_model) {
  p.max <- 5
  q.max <- 5
  bic.min <- Inf
  best.p <- 0
  best.q <- 0
  inf.crit <- 0
  for (i1 in 1:p.max) {
    for (i2 in 1:q.max) {
      ourSpec <-
        ugarchspec(
          mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
          variance.model = list(garchOrder = c(i1, i2), model = "gjrGARCH")
        )
      fit <- ugarchfit(spec = ourSpec, data = resid(arima_model), solver = "hybrid")
      inf.crit <- infocriteria(fit)[2]
      bic.min <- ifelse(inf.crit < bic.min, inf.crit, bic.min)
      
      best.p <- ifelse(inf.crit == bic.min, i1, best.p)
      best.q <- ifelse(inf.crit == bic.min, i2, best.q)
    }
  }
  return(c(best.p, best.q))
}

#print(bestGJRGARCH_AIC(arima202)) #(1,1)
#print(bestGJRGARCH_AIC(arima000)) #(1,1)
#print(bestGJRGARCH_BIC(arima202)) #(1,1)
#print(bestGJRGARCH_BIC(arima000)) #(1,1)

bestEGARCH_AIC <- function(arima_model) {
  p.max <- 5
  q.max <- 5
  aic.min <- Inf
  best.p <- 0
  best.q <- 0
  inf.crit <- 0
  for (i1 in 1:p.max) {
    for (i2 in 1:q.max) {
      ourSpec <-
        ugarchspec(
          mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
          variance.model = list(garchOrder = c(i1, i2), model = "eGARCH")
        )
      fit <- ugarchfit(spec = ourSpec, data = resid(arima_model), solver = "hybrid")
      inf.crit <- infocriteria(fit)[1]
      aic.min <- ifelse(inf.crit < aic.min, inf.crit, aic.min)
      
      best.p <- ifelse(inf.crit == aic.min, i1, best.p)
      best.q <- ifelse(inf.crit == aic.min, i2, best.q)
    }
  }
  return(c(best.p, best.q))
}
bestEGARCH_BIC <- function(arima_model) {
  p.max <- 5
  q.max <- 5
  bic.min <- Inf
  best.p <- 0
  best.q <- 0
  inf.crit <- 0
  for (i1 in 1:p.max) {
    for (i2 in 1:q.max) {
      ourSpec <-
        ugarchspec(
          mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
          variance.model = list(garchOrder = c(i1, i2), model = "eGARCH")
        )
      fit <- ugarchfit(spec = ourSpec, data = resid(arima_model), solver = "hybrid")
      inf.crit <- infocriteria(fit)[2]
      bic.min <- ifelse(inf.crit < bic.min, inf.crit, bic.min)
      
      best.p <- ifelse(inf.crit == bic.min, i1, best.p)
      best.q <- ifelse(inf.crit == bic.min, i2, best.q)
    }
  }
  return(c(best.p, best.q))
}

#print(bestEGARCH_AIC(arima202)) #(5,3)
#print(bestEGARCH_AIC(arima000)) #(4,4)
#print(bestEGARCH_BIC(arima202)) #(5,3)
#print(bestEGARCH_BIC(arima000)) #(4,4)



bestARMA_GARCH_AIC <- function(ts_data) {
  arma_p.max <- 5
  arma_q.max <- 5
  garch_p.max <- 5
  garch_q.max <- 5
  aic.min <- Inf
  arma_best.p <- 0
  arma_best.q <- 0
  garch_best.p <- 0
  garch_best.q <- 0
  inf.crit <- 0
  for (i1 in 1:arma_p.max) {
    for (i2 in 1:arma_q.max) {
      for (i3 in 1:garch_p.max) {
        for (i4 in 1:garch_q.max) {
          ourSpec <-ugarchspec(mean.model = list(armaOrder = c(i1, i2),include.mean = FALSE),
                              variance.model = list(garchOrder = c(i3, i4))
            )
          fit <-ugarchfit(spec = ourSpec, data = ts_data, solver = "hybrid")
          inf.crit <- infocriteria(fit)[1]
          aic.min <- ifelse(inf.crit < aic.min, inf.crit, aic.min)
          arma_best.p <- ifelse(inf.crit == aic.min, i3, arma_best.p)
          arma_best.q <- ifelse(inf.crit == aic.min, i4, arma_best.q)
          garch_best.p <- ifelse(inf.crit == aic.min, i3, garch_best.p)
          garch_best.q <- ifelse(inf.crit == aic.min, i4, garch_best.q)
        }
      }
    }
  }
  return(c(arma_best.p, arma_best.q, garch_best.p, garch_best.q))
}
bestARMA_GARCH_BIC <- function(ts_data) {
  arma_p.max <- 5
  arma_q.max <- 5
  garch_p.max <- 5
  garch_q.max <- 5
  bic.min <- Inf
  arma_best.p <- 0
  arma_best.q <- 0
  garch_best.p <- 0
  garch_best.q <- 0
  inf.crit <- 0
  for (i1 in 1:arma_p.max) {
    for (i2 in 1:arma_q.max) {
      for (i3 in 1:garch_p.max) {
        for (i4 in 1:garch_q.max) {
          ourSpec <-ugarchspec(mean.model = list(armaOrder = c(i1, i2),include.mean = FALSE),
                               variance.model = list(garchOrder = c(i3, i4))
          )
          fit <-ugarchfit(spec = ourSpec, data = ts_data, solver = "hybrid")
          inf.crit <- infocriteria(fit)[2]
          bic.min <- ifelse(inf.crit < bic.min, inf.crit, bic.min)
          arma_best.p <- ifelse(inf.crit == bic.min, i3, arma_best.p)
          arma_best.q <- ifelse(inf.crit == bic.min, i4, arma_best.q)
          garch_best.p <- ifelse(inf.crit == bic.min, i3, garch_best.p)
          garch_best.q <- ifelse(inf.crit == bic.min, i4, garch_best.q)
        }
      }
    }
  }
  return(c(arma_best.p, arma_best.q, garch_best.p, garch_best.q))
}

#bestARMA_GARCH_AIC(BABA_logret) #1 5 1 5
#bestARMA_GARCH_BIC(BABA_logret) #1 1 1 1

bestARMA_GJRGARCH_AIC <- function(ts_data) {
  arma_p.max <- 5
  arma_q.max <- 5
  garch_p.max <- 5
  garch_q.max <- 5
  aic.min <- Inf
  arma_best.p <- 0
  arma_best.q <- 0
  garch_best.p <- 0
  garch_best.q <- 0
  inf.crit <- 0
  for (i1 in 1:arma_p.max) {
    for (i2 in 1:arma_q.max) {
      for (i3 in 1:garch_p.max) {
        for (i4 in 1:garch_q.max) {
          ourSpec <-ugarchspec(mean.model = list(armaOrder = c(i1, i2),include.mean = FALSE),
                               variance.model = list(garchOrder = c(i3, i4), model = "gjrGARCH")
          )
          fit <-ugarchfit(spec = ourSpec, data = ts_data, solver = "hybrid")
          inf.crit <- infocriteria(fit)[1]
          aic.min <- ifelse(inf.crit < aic.min, inf.crit, aic.min)
          arma_best.p <- ifelse(inf.crit == aic.min, i3, arma_best.p)
          arma_best.q <- ifelse(inf.crit == aic.min, i4, arma_best.q)
          garch_best.p <- ifelse(inf.crit == aic.min, i3, garch_best.p)
          garch_best.q <- ifelse(inf.crit == aic.min, i4, garch_best.q)
        }
      }
    }
  }
  return(c(arma_best.p, arma_best.q, garch_best.p, garch_best.q))
}
bestARMA_GJRGARCH_BIC <- function(ts_data) {
  arma_p.max <- 5
  arma_q.max <- 5
  garch_p.max <- 5
  garch_q.max <- 5
  bic.min <- Inf
  arma_best.p <- 0
  arma_best.q <- 0
  garch_best.p <- 0
  garch_best.q <- 0
  inf.crit <- 0
  for (i1 in 1:arma_p.max) {
    for (i2 in 1:arma_q.max) {
      for (i3 in 1:garch_p.max) {
        for (i4 in 1:garch_q.max) {
          ourSpec <-ugarchspec(mean.model = list(armaOrder = c(i1, i2),include.mean = FALSE),
                               variance.model = list(garchOrder = c(i3, i4), model = "gjrGARCH")
          )
          fit <-ugarchfit(spec = ourSpec, data = ts_data, solver = "hybrid")
          inf.crit <- infocriteria(fit)[2]
          bic.min <- ifelse(inf.crit < bic.min, inf.crit, bic.min)
          arma_best.p <- ifelse(inf.crit == bic.min, i3, arma_best.p)
          arma_best.q <- ifelse(inf.crit == bic.min, i4, arma_best.q)
          garch_best.p <- ifelse(inf.crit == bic.min, i3, garch_best.p)
          garch_best.q <- ifelse(inf.crit == bic.min, i4, garch_best.q)
        }
      }
    }
  }
  return(c(arma_best.p, arma_best.q, garch_best.p, garch_best.q))
}

#bestARMA_GJRGARCH_AIC(BABA_logret) # 1 1 1 1
#bestARMA_GJRGARCH_BIC(BABA_logret) # 1 1 1 1

bestARMA_eGARCH_AIC <- function(ts_data) {
  arma_p.max <- 5
  arma_q.max <- 5
  garch_p.max <- 5
  garch_q.max <- 5
  aic.min <- Inf
  arma_best.p <- 0
  arma_best.q <- 0
  garch_best.p <- 0
  garch_best.q <- 0
  inf.crit <- 0
  for (i1 in 1:arma_p.max) {
    for (i2 in 1:arma_q.max) {
      for (i3 in 1:garch_p.max) {
        for (i4 in 1:garch_q.max) {
          ourSpec <-ugarchspec(mean.model = list(armaOrder = c(i1, i2),include.mean = FALSE),
                               variance.model = list(garchOrder = c(i3, i4), model = "eGARCH")
          )
          fit <-ugarchfit(spec = ourSpec, data = ts_data, solver = "hybrid")
          inf.crit <- infocriteria(fit)[1]
          aic.min <- ifelse(inf.crit < aic.min, inf.crit, aic.min)
          arma_best.p <- ifelse(inf.crit == aic.min, i3, arma_best.p)
          arma_best.q <- ifelse(inf.crit == aic.min, i4, arma_best.q)
          garch_best.p <- ifelse(inf.crit == aic.min, i3, garch_best.p)
          garch_best.q <- ifelse(inf.crit == aic.min, i4, garch_best.q)
        }
      }
    }
  }
  return(c(arma_best.p, arma_best.q, garch_best.p, garch_best.q))
}
bestARMA_eGARCH_BIC <- function(ts_data) {
  arma_p.max <- 5
  arma_q.max <- 5
  garch_p.max <- 5
  garch_q.max <- 5
  bic.min <- Inf
  arma_best.p <- 0
  arma_best.q <- 0
  garch_best.p <- 0
  garch_best.q <- 0
  inf.crit <- 0
  for (i1 in 1:arma_p.max) {
    for (i2 in 1:arma_q.max) {
      for (i3 in 1:garch_p.max) {
        for (i4 in 1:garch_q.max) {
          ourSpec <-ugarchspec(mean.model = list(armaOrder = c(i1, i2),include.mean = FALSE),
                               variance.model = list(garchOrder = c(i3, i4), model = "eGARCH")
          )
          fit <-ugarchfit(spec = ourSpec, data = ts_data, solver = "hybrid")
          inf.crit <- infocriteria(fit)[2]
          bic.min <- ifelse(inf.crit < bic.min, inf.crit, bic.min)
          arma_best.p <- ifelse(inf.crit == bic.min, i3, arma_best.p)
          arma_best.q <- ifelse(inf.crit == bic.min, i4, arma_best.q)
          garch_best.p <- ifelse(inf.crit == bic.min, i3, garch_best.p)
          garch_best.q <- ifelse(inf.crit == bic.min, i4, garch_best.q)
        }
      }
    }
  }
  return(c(arma_best.p, arma_best.q, garch_best.p, garch_best.q))
}

bestARMA_eGARCH_AIC(BABA_logret) # 5 3 5 3
bestARMA_eGARCH_BIC(BABA_logret) # 4 3 4 3


#
arma202_garch11_spec <- ugarchspec(mean.model = list(armaOrder = c(2, 2)),
                                  variance.model = list(garchOrder = c(1, 1)))
arma202_garch11 <- ugarchfit(spec = arma202_garch11_spec, data = BABA_logret,
                             out.sample = 183,solver = "hybrid")
Box.test(residuals(arma202_garch11))
jarque.bera.test(residuals(arma202_garch11))
par(mfrow = c(1, 2))
hist(residuals(arma202_garch11), breaks = 30, main ='Histogram', cex.main = 0.8, cex.lab = 0.8, xlab = NA,
     cex.axis = 0.8)
box()
qqnorm(residuals(arma202_garch11), cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8) 
qqline(residuals(arma202_garch11), lwd = 2)
par(mfrow = c(1,1))

arma202_garch14_spec <- ugarchspec(mean.model = list(armaOrder = c(2, 2)),
                                   variance.model = list(garchOrder = c(1, 4)))
arma202_garch14 <- ugarchfit(spec = arma202_garch14_spec, data = BABA_logret,
                             out.sample = 183,solver = "hybrid")
Box.test(residuals(arma202_garch14))
jarque.bera.test(residuals(arma202_garch14))
par(mfrow = c(1, 2))
hist(residuals(arma202_garch14), breaks = 30, main ='Histogram', cex.main = 0.8, cex.lab = 0.8, xlab = NA,
     cex.axis = 0.8)
box()
qqnorm(residuals(arma202_garch14), cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8) 
qqline(residuals(arma202_garch14), lwd = 2)
par(mfrow = c(1,1))

arma000_garch11_spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                                   variance.model = list(garchOrder = c(1, 1)))
arma000_garch11 <- ugarchfit(spec = arma000_garch11_spec, data = BABA_logret,
                             out.sample = 183,solver = "hybrid")
Box.test(residuals(arma000_garch11))
jarque.bera.test(residuals(arma000_garch11))
par(mfrow = c(1, 2))
hist(residuals(arma000_garch11), breaks = 30, main ='Histogram', cex.main = 0.8, cex.lab = 0.8, xlab = NA,
     cex.axis = 0.8)
box()
qqnorm(residuals(arma000_garch11), cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8) 
qqline(residuals(arma000_garch11), lwd = 2)
par(mfrow = c(1,1))

arma000_garch14_spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                                   variance.model = list(garchOrder = c(1, 4)))
arma000_garch14 <- ugarchfit(spec = arma000_garch14_spec, data = BABA_logret,
                             out.sample = 183,solver = "hybrid")
Box.test(residuals(arma000_garch14))
jarque.bera.test(residuals(arma000_garch14))
par(mfrow = c(1, 2))
hist(residuals(arma000_garch14), breaks = 30, main ='Histogram', cex.main = 0.8, cex.lab = 0.8, xlab = NA,
     cex.axis = 0.8)
box()
qqnorm(residuals(arma000_garch14), cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8) 
qqline(residuals(arma000_garch14), lwd = 2)
par(mfrow = c(2,2))
Acf(residuals(arma000_garch14))
Pacf(residuals(arma000_garch14))

arma105_garch15_spec <- ugarchspec(mean.model = list(armaOrder = c(1, 5)),
                                   variance.model = list(garchOrder = c(1, 5)))
arma105_garch15 <- ugarchfit(spec = arma105_garch15_spec, data = BABA_logret,
                             out.sample = 183,solver = "hybrid")
Box.test(residuals(arma105_garch15))
jarque.bera.test(residuals(arma105_garch15))
par(mfrow = c(1, 2))
hist(residuals(arma105_garch15), breaks = 30, main ='Histogram', cex.main = 0.8, cex.lab = 0.8, xlab = NA,
     cex.axis = 0.8)
box()
curve(dnorm(x, mean=mean(residuals(arma105_garch15)), sd=sd(residuals(arma105_garch15))), add=TRUE, col="red")
qqnorm(residuals(arma105_garch15), cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8) 
qqline(residuals(arma105_garch15), lwd = 2)
par(mfrow = c(1,1))

arma101_garch11_spec <- ugarchspec(mean.model = list(armaOrder = c(1, 1)),
                                   variance.model = list(garchOrder = c(1, 1)))
arma101_garch11 <- ugarchfit(spec = arma101_garch11_spec, data = BABA_logret,
                             out.sample = 183,solver = "hybrid")
Box.test(residuals(arma101_garch11))
jarque.bera.test(residuals(arma101_garch11))

arma202_gjrgarch11_spec <- ugarchspec(mean.model = list(armaOrder = c(2, 2)),
                                      variance.model = list(garchOrder = c(1, 1), model  = "gjrGARCH"))
arma202_gjrgarch11 <- ugarchfit(spec = arma202_gjrgarch11_spec, data = BABA_logret,
                                out.sample = 183)
Box.test(residuals(arma202_gjrgarch11))
jarque.bera.test(residuals(arma202_gjrgarch11))

arma000_gjrgarch11_spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                                      variance.model = list(garchOrder = c(1, 1), model  = "gjrGARCH"))
arma000_gjrgarch11 <- ugarchfit(spec = arma000_gjrgarch11_spec, data = BABA_logret,
                                out.sample = 183)
Box.test(residuals(arma202_gjrgarch11))
jarque.bera.test(residuals(arma202_gjrgarch11))

arma101_gjrgarch11_spec <- ugarchspec(mean.model = list(armaOrder = c(1, 1)),
                                      variance.model = list(garchOrder = c(1, 1), model  = "gjrGARCH"))
arma101_gjrgarch11 <- ugarchfit(spec = arma101_gjrgarch11_spec, data = BABA_logret,
                             out.sample = 183)
Box.test(residuals(arma101_gjrgarch11))
jarque.bera.test(residuals(arma101_gjrgarch11))

arma503_egarch53_spec <- ugarchspec(mean.model = list(armaOrder = c(5, 3)),
                                      variance.model = list(garchOrder = c(5, 3), model  = "eGARCH"))
arma503_egarch53 <- ugarchfit(spec = arma503_egarch53_spec, data = BABA_logret,
                                out.sample = 183, solver = "hybrid")
Box.test(residuals(arma503_egarch53))
jarque.bera.test(residuals(arma503_egarch53))

arma403_egarch43_spec <- ugarchspec(mean.model = list(armaOrder = c(4, 3)),
                                    variance.model = list(garchOrder = c(4, 3), model  = "eGARCH"))
arma403_egarch43 <- ugarchfit(spec = arma403_egarch43_spec, data = BABA_logret,
                              out.sample = 183, solver = "hybrid")
Box.test(residuals(arma403_egarch43))
jarque.bera.test(residuals(arma403_egarch43))

arma202_egarch44_spec <- ugarchspec(mean.model = list(armaOrder = c(2, 2)),
                                    variance.model = list(garchOrder = c(4, 4), model  = "eGARCH"))
arma202_egarch44 <- ugarchfit(spec = arma202_egarch44_spec, data = BABA_logret,
                              out.sample = 183, solver = "hybrid")
Box.test(residuals(arma202_egarch44))
jarque.bera.test(residuals(arma202_egarch44))
par(mfrow = c(1, 2))
hist(residuals(arma105_garch15), breaks = 30, main ='Histogram', cex.main = 0.8, cex.lab = 0.8, xlab = NA,
     cex.axis = 0.8)
box()
qqnorm(residuals(arma105_garch15), cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8) 
qqline(residuals(arma105_garch15), lwd = 2)
par(mfrow = c(1,1))

arma202_egarch53_spec <- ugarchspec(mean.model = list(armaOrder = c(2, 2)),
                                    variance.model = list(garchOrder = c(5, 3), model  = "eGARCH"))
arma202_egarch53 <- ugarchfit(spec = arma202_egarch53_spec, data = BABA_logret,
                              out.sample = 183, solver = "hybrid")
Box.test(residuals(arma202_egarch53))
jarque.bera.test(residuals(arma202_egarch53))

arma000_egarch44_spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                                    variance.model = list(garchOrder = c(4, 4), model  = "eGARCH"))
arma000_egarch44 <- ugarchfit(spec = arma000_egarch44_spec, data = BABA_logret,
                              out.sample = 183, solver = "hybrid")
Box.test(residuals(arma000_egarch44))
jarque.bera.test(residuals(arma000_egarch44))
par(mfrow = c(2, 1))
hist(residuals(arma000_egarch44), breaks = 30, main ='Histogram', cex.main = 0.8, cex.lab = 0.8, xlab = NA,
     cex.axis = 0.8)
box()
qqnorm(residuals(arma000_egarch44), cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8) 
qqline(residuals(arma000_egarch44), lwd = 2)
par(mfrow = c(1,1))

arma000_egarch53_spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                                    variance.model = list(garchOrder = c(5, 3), model  = "eGARCH"))
arma000_egarch53 <- ugarchfit(spec = arma000_egarch53_spec, data = BABA_logret,
                              out.sample = 183, solver = "hybrid")
Box.test(residuals(arma000_egarch44))
jarque.bera.test(residuals(arma000_egarch44))

arma000_garch00_spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                                    variance.model = list(garchOrder = c(0, 0), model  = "sGARCH"))
arma000_garch00 <- ugarchfit(spec = arma000_egarch00_spec, data = BABA_logret,
                              out.sample = 183, solver = "hybrid")
Box.test(residuals(arma000_egarch00))
jarque.bera.test(residuals(arma000_egarch00))

#

ic_comp_table <- data.frame(cbind(infocriteria(arma000_garch11), infocriteria(arma000_garch14), 
                                  infocriteria(arma202_garch11), infocriteria(arma202_garch14), 
                                  infocriteria(arma105_garch15), infocriteria(arma101_garch11),
                                  infocriteria(arma101_gjrgarch11),
                                  infocriteria(arma000_gjrgarch11),infocriteria(arma202_gjrgarch11),
                                  infocriteria(arma503_egarch53), infocriteria(arma403_egarch43),
                                  infocriteria(arma202_egarch44), infocriteria(arma202_egarch53),
                                  infocriteria(arma000_egarch44), infocriteria(arma000_egarch53)))
names(ic_comp_table) <- c('ARMA(0,0,0)-GARCH(1,1)','ARMA(0,0,0)-GARCH(1,4)',
                          'ARMA(2,0,2)-GARCH(1,1)','ARMA(2,0,2)-GARCH(1,4)',
                          'ARMA(1,0,5)-GARCH(1,5)','ARMA(1,0,1)-GARCH(1,1)',
                          'ARMA(1,0,1)-gjrGARCH(1,1)',
                          'ARMA(0,0,0)-gjrGARCH(1,1)','ARMA(2,0,2)-gjrGARCH(1,1)',
                          'ARMA(5,0,3)-eGARCH(5,3)', 'ARMA(4,0,3)-eGARCH(4,3)',
                          'ARMA(2,0,2)-eGARCH(4,4)', 'ARMA(2,0,2)-eGARCH(5,3)',
                          'ARMA(0,0,0)-eGARCH(4,4)', 'ARMA(0,0,0)-eGARCH(5,3)')
ic_comp_table[1:2,]
min(ic_comp_table['Akaike',])
min(ic_comp_table['Bayes',])
stargazer(t(ic_comp_table[1:2,]), summary = FALSE)

#### Forecast ####
forecast_arma000_garch11 = ugarchforecast(arma000_garch11, n.ahead = 1, n.roll = 183)
forecast_arma000_garch14 = ugarchforecast(arma000_garch14, n.ahead = 1, n.roll = 183)
forecast_arma202_garch11 = ugarchforecast(arma202_garch11, n.ahead = 1, n.roll = 183)
forecast_arma202_garch14 = ugarchforecast(arma202_garch14, n.ahead = 1, n.roll = 183)
forecast_arma105_garch15 = ugarchforecast(arma105_garch15, n.ahead = 1, n.roll = 183)
forecast_arma101_garch11 = ugarchforecast(arma101_garch11, n.ahead = 1, n.roll = 183)

forecast_arma101_gjrgarch11  = ugarchforecast(arma101_gjrgarch11 , n.ahead = 1, n.roll = 183)
forecast_arma000_gjrgarch11  = ugarchforecast(arma000_gjrgarch11, n.ahead = 1, n.roll = 183)
forecast_arma202_gjrgarch11  = ugarchforecast(arma202_gjrgarch11 , n.ahead = 1, n.roll = 183)

forecast_arma503_egarch53  = ugarchforecast(arma503_egarch53, n.ahead = 1, n.roll = 183)
forecast_arma403_egarch43  = ugarchforecast(arma403_egarch43, n.ahead = 1, n.roll = 183)
forecast_arma202_egarch44  = ugarchforecast(arma202_egarch44, n.ahead = 1, n.roll = 183)
forecast_arma202_egarch53  = ugarchforecast(arma202_egarch53, n.ahead = 1, n.roll = 183)
forecast_arma000_egarch44  = ugarchforecast(arma000_egarch44, n.ahead = 1, n.roll = 183)
forecast_arma000_egarch53  = ugarchforecast(arma000_egarch53, n.ahead = 1, n.roll = 183)

forecast_arma000_garch00 = ugarchforecast(arma000_garch00, n.ahead = 1, n.roll = 183)
# Bind the prediction vectors
predictions  <- cbind(forecast_arma000_garch11@forecast$sigmaFor[1:183],
                      forecast_arma000_garch14@forecast$sigmaFor[1:183],
                      forecast_arma202_garch11@forecast$sigmaFor[1:183],
                      forecast_arma202_garch14@forecast$sigmaFor[1:183],
                      forecast_arma105_garch15@forecast$sigmaFor[1:183],
                      forecast_arma101_garch11@forecast$sigmaFor[1:183],
                      forecast_arma101_gjrgarch11@forecast$sigmaFor[1:183],
                      forecast_arma000_gjrgarch11@forecast$sigmaFor[1:183],
                      forecast_arma202_gjrgarch11@forecast$sigmaFor[1:183],
                      forecast_arma503_egarch53@forecast$sigmaFor[1:183],
                      forecast_arma403_egarch43@forecast$sigmaFor[1:183],
                      
                      forecast_arma202_egarch44@forecast$sigmaFor[1:183],
                      forecast_arma202_egarch53@forecast$sigmaFor[1:183],
                      forecast_arma000_egarch44@forecast$sigmaFor[1:183],
                      forecast_arma000_egarch53@forecast$sigmaFor[1:183],
                      forecast_arma000_garch00@forecast$sigmaFor[1:183])



# Calculate volatility proxy
vol <- (tail(BABA_logret,183))^2
# Proceeds to MSE
MSE_results <- data.frame(cbind(MSE(y_pred = predictions[1], y_true = vol),
                                MSE(y_pred = predictions[2], y_true = vol),
                                MSE(y_pred = predictions[3], y_true = vol),
                                MSE(y_pred = predictions[4], y_true = vol),
                                MSE(y_pred = predictions[5], y_true = vol),
                                MSE(y_pred = predictions[6], y_true = vol),
                                MSE(y_pred = predictions[7], y_true = vol),
                                MSE(y_pred = predictions[8], y_true = vol),
                                MSE(y_pred = predictions[9], y_true = vol),
                                MSE(y_pred = predictions[10], y_true = vol),
                                MSE(y_pred = predictions[11], y_true = vol),
                                MSE(y_pred = predictions[12], y_true = vol),
                                MSE(y_pred = predictions[13], y_true = vol),
                                MSE(y_pred = predictions[14], y_true = vol),
                                MSE(y_pred = predictions[15], y_true = vol),
                                MSE(y_pred = predictions[16], y_true = vol)))


colnames(MSE_results) <- c('ARMA(0,0,0)-GARCH(1,1)',
                           'ARMA(0,0,0)-GARCH(1,4)',
                           'ARMA(2,0,2)-GARCH(1,1)',
                           'ARMA(2,0,2)-GARCH(1,4)',
                           'ARMA(1,0,5)-GARCH(1,5)',
                           'ARMA(1,0,1)-GARCH(1,1)',
                           'ARMA(1,0,1)-gjrGARCH(1,1)',
                           'ARMA(0,0,0)-gjrGARCH(1,1)',
                           'ARMA(2,0,2)-gjrGARCH(1,1)',
                           'ARMA(5,0,3)-eGARCH(5,3)',
                           'ARMA(4,0,3)-eGARCH(4,3)',
                           'ARMA(2,0,2)-eGARCH(4,4)',
                           'ARMA(2,0,2)-eGARCH(5,3)',
                           'ARMA(0,0,0)-eGARCH(4,4)',
                           'ARMA(0,0,0)-eGARCH(5,3)',
                           'ARMA(0,0,0ú-GARCH(0,0)')

MAPE_results <- data.frame(cbind(MAPE(y_pred = predictions[1], y_true = vol),
                                MAPE(y_pred = predictions[2], y_true = vol),
                                MAPE(y_pred = predictions[3], y_true = vol),
                                MAPE(y_pred = predictions[4], y_true = vol),
                                MAPE(y_pred = predictions[5], y_true = vol),
                                MAPE(y_pred = predictions[6], y_true = vol),
                                MAPE(y_pred = predictions[7], y_true = vol),
                                MAPE(y_pred = predictions[8], y_true = vol),
                                MAPE(y_pred = predictions[9], y_true = vol),
                                MAPE(y_pred = predictions[10], y_true = vol),
                                MAPE(y_pred = predictions[11], y_true = vol),
                                MAPE(y_pred = predictions[12], y_true = vol),
                                MAPE(y_pred = predictions[13], y_true = vol),
                                MAPE(y_pred = predictions[14], y_true = vol),
                                MAPE(y_pred = predictions[15], y_true = vol),
                                MAPE(y_pred = predictions[16], y_true = vol)))

colnames(MAPE_results) <- c('ARMA(0,0,0)-GARCH(1,1)',
                           'ARMA(0,0,0)-GARCH(1,4)',
                           'ARMA(2,0,2)-GARCH(1,1)',
                           'ARMA(2,0,2)-GARCH(1,4)',
                           'ARMA(1,0,5)-GARCH(1,5)',
                           'ARMA(1,0,1)-GARCH(1,1)',
                           'ARMA(1,0,1)-gjrGARCH(1,1)',
                           'ARMA(0,0,0)-gjrGARCH(1,1)',
                           'ARMA(2,0,2)-gjrGARCH(1,1)',
                           'ARMA(5,0,3)-eGARCH(5,3)',
                           'ARMA(4,0,3)-eGARCH(4,3)',
                           'ARMA(2,0,2)-eGARCH(4,4)',
                           'ARMA(2,0,2)-eGARCH(5,3)',
                           'ARMA(0,0,0)-eGARCH(4,4)',
                           'ARMA(0,0,0)-eGARCH(5,3)',
                           'ARMA(0,0,0)-GARCH(0,0)')
MAPE_results <- MAPE_results[,order(MAPE_results)]
t(MAPE_results)
#reordering (lowest RMSE on the left)
MSE_results <- MSE_results[,order(MSE_results)]
t(MSE_)
stargazer(t(MSE_results), summary  = FALSE, digits = 10)
t(MSE_results)[,1] * 1000

stargazer::stargazer(arma000_egarch44@fit$matcoef, 
                     title = "Parameter Estimates of the GARCH(1, 1)") %>% 
  gsub("Std. Error", "Rob. Std. Error", .) %>%  
  gsub("t value", "Rob. t value", .) %>%  
  gsub("mu", "$\\\\mu$", .) %>%
  gsub("alpha1", "$\\\\alpha$", .) %>%
  gsub("omega", "$\\\\omega$", .) %>%  
  gsub("beta1", "$\\\\beta$", .) %>%
  gsub("shape", "$\\\\nu$", .)  %>%
  writeLines("arch_output.tex")
