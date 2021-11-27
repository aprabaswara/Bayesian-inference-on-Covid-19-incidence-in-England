##Aditya Prabaswara Mardjikoen (S2264710)

library(rjags)
library(coda)

##data on daily hospital deaths with COVID-19 in England for the first 100 days from March 2nd 2020,from NHS England.

y <-c(1,2,0,2,2,0,4,4,1,9,14,20,22,27,40,46,66,63,105,103,149,159,204,263,326,353,360,437,498,576,645,647,
      700,778,743,727,813,900,792,740,779,718,699,646,686,639,610,571,524,566,486,502,451,438,386,380,345,341,
      325,313,306,268,251,259,252,266,256,215,203,196,166,183,162,180,172,167,138,160,144,153,150,122,128,116,
      133,139,122,124,117,92,83,94,109,111,83,86,83,80,73,69)

##new data starting from February 11, 2020
y_new <- c(rep(0,20),y)

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
sam <- jags.samples(mod,c("m","n"),n.iter=10000)
str(sam)
sam.coda <- coda.samples(mod,c("tau"),n.iter=10000)
str(sam.coda)
plot(sam.coda)
acfplot(sam.coda,aspect=1)
effectiveSize(sam.coda)
HPDinterval(sam.coda[[1]])
apply(sam.coda[[1]],2,quantile,prob=(c(0.025,0.975)))

##generate credible interval
credible.int <- apply(sam$n,1,quantile,prob=(c(0.025,0.975)))
max_vertical <- max(sam$n[1:length(y)],credible.int[1,1:length(y)],credible.int[2,1:length(y)],sam$m[1:length(y_new)])
lockdown_day <- julian(as.Date("2020-3-24"),origin=as.Date("2019-12-31"))
date_vector <- seq(as.Date("2020-2-11"), by = "days", length.out = 120)
julian_day <- julian(date_vector,origin=as.Date("2019-12-31"))

##generate start date
start_date <- c()
observed_month <- unique(months(date_vector))
for (i in observed_month){
  first_date <- min(date_vector[which(months(date_vector)==i)])
  first_date <- as.character(first_date)
  start_date <- c(start_date,first_date)
}
start_date <- as.Date(start_date)
julian_start <- julian(start_date,origin=as.Date("2019-12-31"))

##plot the data
plot(julian_day,y_new,xlab='Day of the year',ylab='Incidence',ylim=c(0,max_vertical),col='grey')
lines(julian_day[1:length(y)],credible.int[1,1:length(y)],col='red',lty=2)
lines(julian_day[1:length(y)],credible.int[2,1:length(y)],col='red',lty=2)
lines(julian_day[1:length(y)],sam$n[1:length(y)],col='green')
lines(julian_day,sam$m[1:length(y_new)],col='blue')

##Question: After running MCMC do we take 120 data in the final iteration?

##lines(julian_day[1:length(y)],sam$n[(length(sam$n)-119):length(sam$n)][1:100],col='green')
##lines(julian_day,sam$m[(length(sam$m)-119):length(sam$m)],col='blue')

abline(v=lockdown_day,lty=2)
axis(3,at=julian_start,labels=observed_month,cex.axis=0.55,lwd.ticks=0.5,padj=1)
legend(x='topright',legend=c("New Infections", "Expected Death", "95% Credible Interval"),col=c("green", "blue","red"), lty=c(1,1,2), cex=0.55,bty='o',inset=0.05)
text(x=lockdown_day+4,y=1100,label='UK Lockdown',pos=3,cex=0.7,srt=90)
title('Daily Death and New Infections from COVID-19 at England (2020)',cex.main=0.8)

##Question: Does new infection decrease sharply before lockdown?