##Aditya Prabaswara Mardjikoen (S2264710)

##Overview:
##------------------------------------------------------------------------------------------------------------------------------------
##The aims of this project is:
##1. Observe and interpret the trend of the daily new infections and expected deaths of COVID-19 in England hospital at 2020.
##2. Perform diagnostic analysis on the new infections and expected deaths
##------------------------------------------------------------------------------------------------------------------------------------
library(rjags)
library(coda)

##Data on daily hospital deaths with COVID-19 in England for the first 100 days from March 2nd 2020,from NHS England.

y <-c(1,2,0,2,2,0,4,4,1,9,14,20,22,27,40,46,66,63,105,103,149,159,204,263,326,353,360,437,498,576,645,647,
      700,778,743,727,813,900,792,740,779,718,699,646,686,639,610,571,524,566,486,502,451,438,386,380,345,341,
      325,313,306,268,251,259,252,266,256,215,203,196,166,183,162,180,172,167,138,160,144,153,150,122,128,116,
      133,139,122,124,117,92,83,94,109,111,83,86,83,80,73,69)

##Recall that there were no COVID-19 deaths prior to March 2nd 2020. So for modelling purposes it makes sense to extend the data
## by adding 20 days of zero deaths to the start of the data in order to observe the data from February 11, 2020.

y_new <- c(rep(0,20),y) ##Daily deaths data from February 11, 2020

B <- matrix(0,nrow=length(y_new),ncol=length(y_new))

for (i in 1:nrow(B)){
  for (j in 1:ncol(B)){
    if (j > i){
      B[i,j] <- 0
    }
    else{
      B[i,j] <- dlnorm(i-j, meanlog = 3.235, sdlog = 0.4147)
    }
  }
}

##Question: Do we include m and n in MCMC monitoring diagnostic like tau?

##input to jags
setwd('C:/Users/Aditya Prabaswara/Bayesian-inference-on-Covid-19-incidence-in-England')
mod <- jags.model("model.jags",data=list(y=y_new,N=length(y_new),B=B))
sam.coda <- coda.samples(mod,c("m","n"),n.iter=10000)

##diagnostic plots
par(mfrow=c(2,3))
traceplot(sam.coda[[1]][,c('n[2]','n[50]','n[100]','m[2]','m[50]','m[120]')])
acfplot(sam.coda[[1]][,c('n[2]','n[50]','n[100]','m[2]','m[50]','m[120]')],aspect=1,type='l')

##effective sample size
recommended_sample <- effectiveSize(sam.coda[[1]])
max(round(recommended_sample,0))

##After we run the code, we find a burn-in period from the first iterations until around the
##4000 iterations. So, we recommend to do a burn-in for 4000 sample in the first iteration
##for the simulations. Because the sampling is quite random, we would suggest that the effective
##sample size should be around 400 until 800 sample.

##Find sample in MCMC sample list for expected death and new infections
death_index <- grep("m", colnames(sam.coda[[1]]))
infection_index <- grep("n", colnames(sam.coda[[1]]))[1:length(y)]

##Generate credible interval for new infections
credible_interval <- apply(sam.coda[[1]][,infection_index],2,quantile,prob=(c(0.025,0.975)))
upper_bound <- credible_interval[2,]
lower_bound <- credible_interval[1,]

##Posterior mean for new infections and expected death
expect_death <- colMeans(sam.coda[[1]][,death_index])
new_infect <- colMeans(sam.coda[[1]][,infection_index])
max_vertical <- max(new_infect,expect_death,upper_bound,lower_bound)
lockdown_day <- julian(as.Date("2020-3-24"),origin=as.Date("2019-12-31"))
date_vector <- seq(as.Date("2020-2-11"), by = "days", length.out = 120)
julian_day <- julian(date_vector,origin=as.Date("2019-12-31"))

##plot the data
#Reference: 
#https://stackoverflow.com/questions/14069629/how-can-i-plot-data-with-confidence-intervals
#https://statisticsglobe.com/r-polygon-function-plot/
#https://statisticsglobe.com/rev-r-function
par(mfrow=c(1,1))
plot(x=julian_day,y=y_new,xlab='Day of the year',ylab='Number of Individuals',ylim=c(0,max_vertical),col='grey',pch=20)
polygon(x=c(rev(julian_day[1:length(y)]), y=julian_day[1:length(y)]), c(rev(lower_bound),upper_bound), col = 'gray91', border = NA)
lines(x=julian_day[1:length(y)],y=lower_bound,col='red',lty=2)
lines(x=julian_day[1:length(y)],y=upper_bound,col='red',lty=2)
lines(x=julian_day[1:length(y)],y=new_infect,col='green')
lines(x=julian_day,y=expect_death,col='blue')

abline(v=lockdown_day,lty=2)
legend(x='topright',legend=c("Mean for New Infections", "Mean for Expected Death", "95% Credible Interval (New Infections)","Actual Daily Death"),
       col=c("green", "blue",NA,'grey'), lty=c(1,1,NA,NA),fill = c(NA,NA, 'gray91',NA),border = c(NA,NA,'gray91',NA), pch=c(NA,NA,NA,20),
       cex=0.55,bty='n',inset=0.05,x.intersp=c(2,2,1.5,2.2))
text(x=lockdown_day+4,y=1100,label='lockdown I',pos=3,cex=0.7,srt=90)
title('Daily Death and New Infections from COVID-19 at England (2020)',cex.main=0.8)
