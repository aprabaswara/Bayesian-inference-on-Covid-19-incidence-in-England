##Aditya Prabaswara Mardjikoen (S2264710)

##Overview:
##------------------------------------------------------------------------------------------------------------------------------------
##The aims of this project is:
##1. Observe and interpret the trend of the daily new infections and expected deaths of COVID-19 in England hospital at 2020.
##2. Perform diagnostic analysis on the new infections and expected deaths
##------------------------------------------------------------------------------------------------------------------------------------
##
##Model assumption:
##------------------------------------------------------------------------------------------------------------------------------------
##1. The observed number of deaths in day i has Poisson distribution.
##2. The 
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

##Create matrix B
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


##Generate posterior random sample from 10000 iterations using JAGS
mod <- jags.model("model.jags",data=list(y=y_new,N=length(y_new),B=B))
sam.coda <- coda.samples(mod,c("m","n"),n.iter=10000)

##Diagnostic plots
par(mfrow=c(4,2),mar=c(2,2,2,2))
traceplot(sam.coda[[1]][,c('n[2]','n[50]','n[60]','n[100]','m[2]','m[50]','m[60]','m[120]')])
acfplot(sam.coda[[1]][,c('n[2]','n[50]','n[60]','n[100]','m[2]','m[50]','m[60]','m[120]')],aspect=1,type='l')

##From the trace plot we can see that the chain for m and n have a good mixing near its peaks compare near the
##beginning and the end. In addition, we also found out that the chain autocorrelation reduces faster near
##its peaks compare to near the beginning and the end. 

##Recommended sample for iterations
recommended_sample <- max(effectiveSize(sam.coda[[1]]))
recommended_sample

##After we run the code with different choice of m and n, we find a burn-in period from the first iterations until around the 4000 iterations. So, we recommend 
##to do a burn-in for 4000 sample in the first iteration for the simulations. Because the sampling is quite random, we would suggest that the effective sample size 
##should be around 400 until 800 sample. Furthermore, regarding the recommended number of sampling iteration to run,we found out that when we run for 30000 iterations 
##we get the effective sample size is not too small and the autocorrelation reduces faster compare to using 10000 iterations.

##Find sample for expected death and new infections from the sample list

death_index <- grep("m", colnames(sam.coda[[1]]))

infection_index <- grep("n", colnames(sam.coda[[1]]))[1:length(y)]

##Get credible interval for new infections

credible_interval <- apply(sam.coda[[1]][,infection_index],2,quantile,prob=(c(0.025,0.975)))

upper_bound <- credible_interval[2,]

lower_bound <- credible_interval[1,]

##Posterior mean for new infections and expected death
expect_death <- colMeans(sam.coda[[1]][,death_index])

new_infect <- colMeans(sam.coda[[1]][,infection_index])

max_vertical <- max(max(new_infect),max(expect_death),max(upper_bound),max(lower_bound),max(y_new))
lockdown_day <- julian(as.Date("2020-3-24"),origin=as.Date("2019-12-31"))
date_vector <- seq(as.Date("2020-2-11"), by = "days", length.out = 120)
julian_day <- julian(date_vector,origin=as.Date("2019-12-31"))
infection_day <- julian_day[1:(length(y_new)-20)]

##In the single summary plot we will plot the 95% credible interval for new infections, the
##Reference: 
#https://stackoverflow.com/questions/14069629/how-can-i-plot-data-with-confidence-intervals
#https://statisticsglobe.com/r-polygon-function-plot/
#https://statisticsglobe.com/rev-r-function

##Plot single summary plot
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))

plot(x=julian_day,y=y_new,xlab='Day of the year',ylab='Number of Individuals',ylim=c(0,max_vertical),col='grey',pch=20)

polygon(x=c(rev(infection_day), infection_day), y=c(rev(lower_bound),upper_bound), col = 'gray91', border = NA) ##Draw credible interval region

lines(x=infection_day,y=lower_bound,col='grey',lty=2)

lines(x=infection_day,y=upper_bound,col='grey',lty=2)

lines(x=infection_day,y=new_infect,col='green') ##Plot posterior mean for new infections

lines(x=julian_day,y=expect_death,col='blue') ##Plot posterior mean for expected deaths


##Highlight the first day of UK lockdown (24th March 2020)

abline(v=lockdown_day,lty=2) 

text(x=lockdown_day+4,y=1100,label='lockdown I',pos=3,cex=0.6,srt=90) 

##Add legend and title to the single summary plot

legend(x='topright',legend=c("Mean for New Infections", "Mean for Expected Death", "95% Credible Interval (New Infections)","Actual Daily Death"),
       col=c("green", "blue",NA,'grey'), lty=c(1,1,NA,NA),fill = c(NA,NA, 'gray91',NA),border = c(NA,NA,'gray91',NA), pch=c(NA,NA,NA,20),
       cex=0.55,bty='n',x.intersp=c(2,2,1.5,2.2))

title('Daily Death and New Infections from COVID-19 at England (2020)',cex.main=0.8) ##Add plot title

##Conclusion:
##From the single summary plot we can see that the posterior mean for the expected death fits well with the daily
##death data from NHS. 