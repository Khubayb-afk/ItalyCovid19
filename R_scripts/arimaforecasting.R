library(forecast)
library(MLmetrics)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)


#Read in Data Set#
italycovid <- read.csv("data/timeseriesdata.csv")

#Create 2 Subset Data sets#
names(italycovid)
#All records from 2020
italycovid20 <- subset(italycovid, y<=311)
#All records from 2021
italycovid21 <- subset(italycovid, y>311)


#Plot time series for all records, cases and deaths#
tsitalycovid19d <-ts (italycovid20$Deaths, frequency = 1,start = c(25/02/2020,1))
tsitalycovid19c <-ts (italycovid20$Cases, frequency = 1,start = c(25/02/2020,1))
#plot(tsitalycovid19c)
#plot(tsitalycovid19d)

#Registered cases forcasting
autoarima1 <- auto.arima(tsitalycovid19c)
forecast1 <- forecast(autoarima1, h=60)
#calculate MAPE
cm <- MAPE(forecast1$mean, italycovid21$Cases)
plot(forecast1,ylab='Cumalative number of cases',xlab='Days from 25/02/2020')
plot(forecast1,ylab='Cumalative number of cases',xlab='Days from 25/02/2020')
     lines(italycovid$Cases, col='purple')
     text(100,y=3000000, paste0('MAPE=',round(cm*100,3),'%'),adj=1)
#plot(forecast1$residuals)
#qqnorm(forecast1$residuals)

#Registered deaths forcasting
autoarima2 <- auto.arima(tsitalycovid19d)
forecast2 <- forecast(autoarima2, h=60)
#calculate MAPE
dm <- MAPE(forecast2$mean, italycovid21$Death)
plot(forecast2,ylab='Cumalative number of deaths',xlab='Days from 25/02/2020')
plot(forecast2,ylab='Cumalative number of deaths',xlab='Days from 25/02/2020')
     lines(italycovid$Deaths, col='purple')
     text(100,y=100000, paste0('MAPE=',round(dm*100,3),'%'),adj=1)
#plot(forecast1$residuals)
#qqnorm(forecast1$residuals)

#Registered cases forcasting for March and April
tsitalycovid19dd <-ts (italycovid$Deaths, frequency = 1,start = c(25/02/2020,1))
tsitalycovid19cc <-ts (italycovid$Cases, frequency = 1,start = c(25/02/2020,1))

#Registered cases forcasting
autoarima3 <- auto.arima(tsitalycovid19cc)
forecast3 <- forecast(autoarima3, h=60)
#calculate MAPE
plot(forecast3,ylab='Cumalative number of cases',xlab='Days from 25/02/2020')
#plot(forecast1$residuals)
#qqnorm(forecast1$residuals)

#Registered deaths forcasting
autoarima4 <- auto.arima(tsitalycovid19dd)
forecast4 <- forecast(autoarima4, h=60)
#calculate MAPE
plot(forecast4,ylab='Cumalative number of deaths',xlab='Days from 25/02/2020')
#plot(forecast1$residuals)
#qqnorm(forecast1$residuals)

write.csv(forecast3, file = "output/arimacases.csv")
write.csv(forecast4, file = "output/arimadeaths.csv")