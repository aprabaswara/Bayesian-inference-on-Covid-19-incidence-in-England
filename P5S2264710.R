##Aditya Prabaswara Mardjikoen (S2264710)

library(rjags)
library(coda)

##data on daily hospital deaths with COVID-19 in England for the first 100 days from March 2nd 2020,from NHS England.

y <-c(1,2,0,2,2,0,4,4,1,9,14,20,22,27,40,46,66,63,105,103,149,159,204,263,326,353,360,437,498,576,645,647,
      700,778,743,727,813,900,792,740,779,718,699,646,686,639,610,571,524,566,486,502,451,438,386,380,345,341,
      325,313,306,268,251,259,252,266,256,215,203,196,166,183,162,180,172,167,138,160,144,153,150,122,128,116,
      133,139,122,124,117,92,83,94,109,111,83,86,83,80,73,69)

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

##input to jags
setwd('C:/Users/Aditya Prabaswara/Bayesian-inference-on-Covid-19-incidence-in-England')
mod <- jags.model("model.jags",data=list(y=y_new,N=length(y_new),B=B))
sam <- jags.samples(mod,c("m","n"),n.iter=10000)
str(sam)
sam.coda <- coda.samples(mod,c("tau"),n.iter=10000)
str(sam.coda)
plot(sam.coda)
acfplot(sam.coda,aspect=1)
crosscorr(sam.coda)
effectiveSize(sam.coda)
HPDinterval(sam.coda[[1]])
apply(sam.coda[[1]],2,quantile,prob=(c(0.025,0.975)))

##plot for death and incidence
credible.int <- apply(sam$n,1,quantile,prob=(c(0.025,0.975)))
max_vertical <- max(sam$n[1:length(y)],credible.int[1,1:length(y)],credible.int[2,1:length(y)],sam$m[1:length(y_new)])
plot(sam$n[1:length(y)],xlab='Day of the year',ylab='Number of patient',type='l',xlim=c(0,length(y_new)),ylim=c(0,max_vertical),col='green',xaxt='n')
lines(credible.int[1,1:length(y)],col='red',lty=2)
lines(credible.int[2,1:length(y)],col='red',lty=2)
points(y_new,pch=1,col='grey')
lines(sam$m[1:length(y_new)],col='blue')
lockdown_day <- julian(as.Date("2020-3-24"),origin=as.Date("2019-12-31"))
horiz_axis <- c(1,20,40,60,lockdown_day,length(y),length(y_new))
abline(v=lockdown_day,lty=2)
axis(1,at=sort(unique(horiz_axis)))
legend(x='top',legend=c("New Infections", "Death", "Credible Interval"),col=c("green", "blue","red"), lty=c(1,1,2), cex=0.55,bty='n',lwd=1)
text(x=lockdown_day+1,y=1500,label='UK Lockdown \n (24 March 2020)',pos=1,cex=0.7,srt=90)
title('Daily Death and New Infections at England',cex.main=1)