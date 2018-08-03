#####     Agustin Palao
#####     Policy Intervention - Interrupted Time Series
#####     Credit Guarantees


require(foreign)
library(data.table)
require(tseries)
require(forecast)
require(TSA)
library(urca)

####(1)Get the data ####
count<-read.csv("Count.csv")
count<-as.data.frame(count)

####(2) Interrupted time series analysis ####
#Dataset for Number of Guarantees and time series plot
Gts <- ts(count$Freq, start=c(2004, 1), end=c(2013, 12), frequency=12)
plot(Gts)
plot(decompose(Gts))
# #Transform the data if necessary
LGts<-log(Gts)
plot(LGts)
plot(decompose(LGts))
# #Lets go with the transformation and figure out how to interpret later
# #We can see a non-stationary behavior so lets difference the ts
DLGts<-diff(LGts)
plot(DLGts)
# #running a Dickey-Fuller test for stationarity
adf.test(DLGts)
# #We can see that now the series is stationary
# #Next, lets do the model specification. Run ACF and PACF
acf(DLGts)
pacf(DLGts)
# #There is seasonal component in lag 12 and we can include MA(1)
fit1<-arima(Gts, order = c(0,1,1), seasonal = list(order=c(0,0,2),period=12),include.drift=TRUE)
fit1
acf(fit1$residuals)
pacf(fit1$residuals)
tsdiag(fit1,gof.lag = 36)
# 
fit2<-Arima(LGts, order = c(0,1,2), seasonal = list(order=c(0,0,2),period=12))
fit2
acf(fit2$residuals)
pacf(fit2$residuals)
tsdiag(fit2,gof.lag = 42)
# 
fit3<-arima(DGts, order = c(0,0,1), seasonal = list(order=c(0,0,2),period=12))
fit3
acf(fit3$residuals)
pacf(fit3$residuals)
tsdiag(fit3,gof.lag = 36)

# I will keep fit1 as the model specification


#Model de preintervention series
preGts<-ts(Gts[1:51], start=c(2004, 1), end=c(2008, 03), frequency=12)
plot(preGts)
acf(preGts, main="ACF Number of Guaranteed Credits")
pacf(preGts, main="PACF Number of Guaranteed Credits")
summary(ur.df(preGts, type="none", selectlags = "AIC"))

#difference seasonality
plot(diff(preGts, 12))
#difference 1st period
plot(diff(diff(preGts,12)))
acf(diff(diff(preGts,12)),main="ACF Number of Guaranteed Credits \n(differenced and seasonal (12) differenced)")
pacf(diff(diff(preGts,12)),main="PACF Number of Guaranteed Credits \n(differenced and seasonal (12) differenced)")
DpreGts<-(diff(diff(preGts,12)))
Box.test(preGts, lag=20, type="Ljung-Box")
summary(ur.df(DpreGts, type = "none", selectlags = "AIC"))
adf.test(DpreGts)
#looks like it needs a ma term
fitpreGts<-Arima(preGts, order=c(0,1,1), seasonal=list(order=c(0,1,0), period=12))
fitpreGts
acf(fitpreGts$residuals)
pacf(fitpreGts$residuals)
par(oma=c(0,0,2,0))
tsdiag(fitpreGts)
title("Model Residual Diagnostics", outer = TRUE)

tsdiag.Arima(fitpreGts)

fitpreGts2<-auto.arima(preGts)
tsdiag(fitpreGts2)

#add an intervention variable
#for a step funtion
myvec<-c(rep(0,51),rep(1,69))
I <- ts(myvec, start=c(2004, 1), end=c(2013, 12), frequency=12)
#for a pulse function
myvec2<-c(rep(0,51),1,rep(0,68))
I2 <- ts(myvec2, start=c(2004, 1), end=c(2013, 12), frequency=12)
#or
I2<-diff(I)

#Estimate intervention analysis
# mod.1 <- arimax(Gts, order=c(0,1,1), seasonal=list(order=c(0,0,2),period=12), xtransf=I, transfer=list(c(1,0)))
# mod.1
# DGts<-diff(Gts)
# mod.1a <- arimax(DGts, order=c(0,0,1), seasonal=list(order=c(0,0,2),period=12), xtransf=I2, transfer=list(c(0,0)))
# mod.1a
# DLGts<-log(DGts)
# mod.1b <- arimax(DLGts, order=c(0,0,2), seasonal=list(order=c(0,0,2),period=12), xtransf=I2, transfer=list(c(0,0)))
# mod.1b
# mod.2 <- arimax(LGts, order=c(0,1,2), seasonal=list(order=c(0,0,2),period=12), xtransf=I, transfer=list(c(0,0)))
# mod.2
mod.fit <- arimax(Gts, order=c(0,1,1), seasonal=list(order=c(0,0,2),period=12), xtransf=I, transfer=list(c(1,0)))
mod.fit
tsdiag(mod.fit)
plot(Gts)
lines(fitted(mod.fit), col="blue")

y.pred <- 505.41*I
  -462.746*(0.0101^(seq(Gts)-52)*as.numeric(seq(Gts)>52))
plot(y=pre,x=seq(pre), type="l", xlab="Time")
plot(y=Gts,x=seq(Gts),type="l",xlab="Time")
plot(y=diff(Gts),x=c(1:119),type="l",xlab="Time")
lines(y=y.pred,x=seq(Gts),col="blue", lwd=2)


# TIME SERIES FOR THE MONEY AMOUNT
#Dataset for Balances
Balance<-read.csv("Balance.csv")
Balance<-as.data.frame(Balance)
Balts <- ts(Balance$x, start=c(2004, 1), end=c(2013, 12), frequency=12)
plot(Balts)

#Deflate Balances
Inflation<-read.csv("Inflation.csv")
GIF<- Inflation$GIF_12.2011[13:132]
Gts<-(Balts*GIF)/1000000
plot(Gts)

#Create the pre-intervention series
preBal<-ts(Gts[1:51], start=c(2004, 1), end=c(2008, 03), frequency=12)
plot(preBal)
acf(preBal, main="ACF Guaranteed Credit Lines Balances")
pacf(preBal, main="PACF Guaranteed Credit Lines Balances")
summary(ur.df(preBal, type = "none", selectlags = "AIC"))

acf(diff(diff(preBal,12)),main="ACF Guaranteed Credit Lines Balances \n(differenced and seasonal (12) differenced)")
pacf(diff(diff(preBal,12)),main="PACF Guaranteed Credit Lines Balances \n(differenced and seasonal (12) differenced)")
DpreBal<-(diff(diff(preBal,12)))
summary(ur.df(DpreBal, type = "none", selectlags = "AIC"))

fitpreBal<-Arima(preBal, order=c(0,1,1), seasonal=list(order=c(0,1,1), period=12))
fitpreBal
tsdiag(fitpreBal, gof.lag = 36)
title("Model Residual Diagnostics", outer = TRUE)

fitpreBal2<-auto.arima(preBal, allowdrift = FALSE)
fitpreBal2
tsdiag(fitpreBal, gof.lag = 36)


#using forecast
summary(forecast(fitpreBal, h=69,level=c(50,80,99)))
ffitpreBal<-forecast(fitpreBal, h=69,level=c(50,80,99))
plot.forecast(ffitpreBal)
lines(Gts)

# TIME SERIES FOR MONEY PER CREDIT LINE
Gramt<-Gts/Gts
preGramt<-ts(Gramt[1:51], start=c(2004, 1), end=c(2008, 03), frequency=12)
plot(preGramt)
acf(preGramt, main="ACF Amount Guaranteed per Credit Line")
pacf(preGramt, main="PACF Amount Guaranteed per Credit Line")
summary(ur.df(preGramt, type = "none", selectlags = "AIC"))

acf(diff(diff(preGramt,12)),main="ACF Amount Guaranteed per Credit Line \n(differenced and seasonal (12) differenced)")
pacf(diff(diff(preGramt,12)),main="PACF Amount Guaranteed per Credit Line \n(differenced and seasonal (12) differenced)")
DpreGramt<-(diff(diff(preGramt,12)))
summary(ur.df(DpreGramt, type = "none", selectlags = "AIC"))

fitpreGramt<-Arima(preGramt, order=c(1,1,1), seasonal=list(order=c(0,1,1), period=12))
fitpreGramt
tsdiag(fitpreGramt, gof.lag = 36)
title("Model Residual Diagnostics", outer = TRUE)


fitpreGramt<-auto.arima(preGramt)
fitpreGramt

#Using forecast for the three models
summary(forecast(fitpreGts, h=69,level=c(50, 95,99)))
summary(forecast(fitpreBal, h=69,level=c(50, 95,99)))

ffitpreGts<-forecast(fitpreGts, h=69, level=c(50, 95,99))
ffitpreBal<-forecast(fitpreBal, h=69, level=c(50, 95,99))
ffitpreGramt<-forecast(fitpreGramt, h=69, level=c(50, 95,99))
par(mfrow=c(3,1))
plot(ffitpreGts, shaded=TRUE, pi.lty=c(3,4,5), fcol='white', shadecols='oldstyle', pi.col=c(2,3,4), main="Forecast for the Number of Guaranteed Credit Lines")
plot(ffitpreBal, shaded=TRUE, pi.lty=c(3,4,5), fcol='white', shadecols='oldstyle', pi.col=c(2,3,4), main="Forecast for the Amount Guaranteed")
plot(ffitpreGramt, shaded=TRUE, pi.lty=c(3,4,5), fcol='white', shadecols='oldstyle', pi.col=c(2,3,4), main="Forecast for the Amount Guaranteed per Guaranteed Credit Line")


#INTERVENTION TIME SERIES ANALYSIS
dev.off()
plot(ffitpreGts, shaded=TRUE, pi.lty=c(3,4,5), fcol='white', shadecols='oldstyle', pi.col=c(2,3,4), main="Intervention Analysis \nNumber of Guaranteed Credit Lines")
postGts<-ts(Gts[52:120], start=c(2008, 4), end=c(2013, 12), frequency=12)
lines(postGts, col="blue", lwd=2)

plot(ffitpreBal, shaded=TRUE, pi.lty=c(3,4,5), fcol='white', shadecols='oldstyle', pi.col=c(2,3,4), main="Intervention Analysis \nBalances of Guaranteed Credit Lines",ylab="MXN millions",xlab="year")
postBal<-ts(Gts[52:120], start=c(2008, 4), end=c(2013, 12), frequency=12)
lines(postBal, col="blue", lwd=2)

plot(ffitpreGramt, shaded=TRUE, pi.lty=c(3,4,5), fcol='white', shadecols='oldstyle', pi.col=c(2,3,4), main="Intervention Analysis \nAverage Amount per Guaranteed Credit Line",ylab="MXN millions",xlab="year")
postGramt<-ts(Gramt[52:120], start=c(2008, 4), end=c(2013, 12), frequency=12)
lines(postGramt, col="blue", lwd=2)



year<-c("2008-2009","2009-2010","2010-2011","2011-2012","2012-2013")
GtsGwthf<-c((ffitpreGts$mean[13]-ffitpreGts$mean[1])/ffitpreGts$mean[1],
           (ffitpreGts$mean[25]-ffitpreGts$mean[13])/ffitpreGts$mean[13],
           (ffitpreGts$mean[37]-ffitpreGts$mean[25])/ffitpreGts$mean[25],
           (ffitpreGts$mean[49]-ffitpreGts$mean[37])/ffitpreGts$mean[37],
           (ffitpreGts$mean[61]-ffitpreGts$mean[49])/ffitpreGts$mean[49]
           #,(ffitpreGts$mean[69]-ffitpreGts$mean[61])/ffitpreGts$mean[61]
           )
GtsGwthr<-c((Gts[64]-Gts[52])/Gts[52],
            (Gts[76]-Gts[64])/Gts[64],
            (Gts[88]-Gts[76])/Gts[76],
            (Gts[100]-Gts[88])/Gts[88],
            (Gts[112]-Gts[100])/Gts[100]
            # ,(Gts[120]-Gts[112])/Gts[112]
            )

BalGwthf<-c((ffitpreBal$mean[13]-ffitpreBal$mean[1])/ffitpreBal$mean[1],
            (ffitpreBal$mean[25]-ffitpreBal$mean[13])/ffitpreBal$mean[13],
            (ffitpreBal$mean[37]-ffitpreBal$mean[25])/ffitpreBal$mean[25],
            (ffitpreBal$mean[49]-ffitpreBal$mean[37])/ffitpreBal$mean[37],
            (ffitpreBal$mean[61]-ffitpreBal$mean[49])/ffitpreBal$mean[49]
            # ,(ffitpreBal$mean[69]-ffitpreBal$mean[61])/ffitpreBal$mean[61]
            )
BalGwthr<-c((GIFBalances[64]-GIFBalances[52])/GIFBalances[52],
            (GIFBalances[76]-GIFBalances[64])/GIFBalances[64],
            (GIFBalances[88]-GIFBalances[76])/GIFBalances[76],
            (GIFBalances[100]-GIFBalances[88])/GIFBalances[88],
            (GIFBalances[112]-GIFBalances[100])/GIFBalances[100]
            # ,(GIFBalances[120]-GIFBalances[112])/GIFBalances[112]
            )

GraGwthf<-c((ffitpreGramt$mean[13]-ffitpreGramt$mean[1])/ffitpreGramt$mean[1],
            (ffitpreGramt$mean[25]-ffitpreGramt$mean[13])/ffitpreGramt$mean[13],
            (ffitpreGramt$mean[37]-ffitpreGramt$mean[25])/ffitpreGramt$mean[25],
            (ffitpreGramt$mean[49]-ffitpreGramt$mean[37])/ffitpreGramt$mean[37],
            (ffitpreGramt$mean[61]-ffitpreGramt$mean[49])/ffitpreGramt$mean[49]
            # ,(ffitpreGramt$mean[69]-ffitpreGramt$mean[61])/ffitpreGramt$mean[61]
)
GraGwthr<-c((Gramt[64]-Gramt[52])/Gramt[52],
            (Gramt[76]-Gramt[64])/Gramt[64],
            (Gramt[88]-Gramt[76])/Gramt[76],
            (Gramt[100]-Gramt[88])/Gramt[88],
            (Gramt[112]-Gramt[100])/Gramt[100]
            # ,(Gramt[120]-Gramt[112])/Gramt[112]
)


GrowthR<-data.frame(year,GtsGwthf,GtsGwthr, BalGwthf, BalGwthr, GraGwthf, GraGwthr)

library(xtable)
xtable(GrowthR)


steps<-forecast(fittest, h=69)
real<-ts(Gts[52:120], start=c(2008, 4), end=c(2013, 12), frequency=12)
Dsr<-real-steps$mean
Dsr
plot(Dsr)
mean(Dsr)
StepsUp95<-steps$upper[,2]
Dsrup<-real-StepsUp95
plot(Dsrup)
mean(Dsrup)
StepsUp80<-steps$upper[,1]
Dsrup80<-real-StepsUp80
plot(Dsrup80)
mean(Dsrup80)
