##Aditya Prabaswara Mardjikoen (S2264710)

##Overview:
##------------------------------------------------------------------------------
##We know that COVID-19 disease is an pandemic that has threaten our life from
##around December 2019 and it start to spread widely in 2020. It is caused by a
##virus called SARS-COV. If we talk about COVID-19 we should talk a little about 
##death and infections. Infections is the process of invasion and multiplication 
##of microorganism such as bacteria, virus, and parasite that are not normally 
##present in our body, while death is  the permanent cessation of all vital 
##functions of the body including the heartbeat, brain activity (including the 
##brain stem), and breathing. 
##
##This R code is mainly about Bayesian data analysis using the Markov Chain 
##Monte Carlo sampling in R with the help of computation in JAGS by using the 
##daily death data of COVID-19 at a hospital which observed from March 2nd, 2020. 
##The data source is from NHS England. Overall, this code is mostly about 
##The aims of this project is:
##1. Observe and interpret the trend of the daily new infections and expected 
##   deaths of COVID-19 in England hospital at 2020.
##
##2. Perform diagnostic analysis on the new infections and expected deaths.
##
##Overall, the final output of this project is the single summary plot 
##containing the plot of actual daily death of the year, posterior mean of 
##expected daily death, posterior mean of daily new infections, and the 95% 
##credible interval for daily new infections which observed daily in terms 
##of day of the year.
##------------------------------------------------------------------------------
##
##Notation:
##------------------------------------------------------------------------------
##m[i] = expected deaths in day i, y[i] = the observed number of deaths in day i , 
##x[i] = the log of the number of new infections in day i, n[i] = new infections in day i
##------------------------------------------------------------------------------
##
##Model assumption:
##------------------------------------------------------------------------------
##1. y[i] has Poisson distribution with rate m[i].
##
##2. We assume that x[1] has Normal distribution with mean 0 and precision 0.01, 
##   x[2] has Normal distribution with mean x[1] and precision tau, and x[i] has
##   Normal distribution with mean 2*x[i-1]-x[i-2] and precision tau. The 
##   precision tau here has Gamma distribution with shape 4 and rate 0.04.
##
##3. Let D be the fatal disease duration. The log of D has Normal distribution 
##   with mean 3.235 and standard deviation 0.4147, where the mean and standard
##   deviation already in log scale. This distribution is based on data from
##   the ISARIC study of 24000+ hospitalized patients who died from COVID-19.
##   We will use this distribution to build matrix B.
##
##4. n changes fairly smoothly from day to day.
##------------------------------------------------------------------------------
##
##Algorithm:
##------------------------------------------------------------------------------
##1. Add 20 zeros at the start of the data so we can observe the data from 
##   February 11, 2020.
##
##2. Let i and j be the index for row and columns of square matrix B respectively.
##   Create matrix B by set the entry where j > i be 0 while the other is 
##   the lognormal density values of i-j which parameterization follows log D.
##
##3. Pass the data, number of observations, and the matrix B into our JAGS model 
##   from R so the 'model.jags' can generate posterior sample for m and n from 
##   10000 sampling iterations into R. Note that computation for m and n is done
##   in JAGS. Further detail about how the sampling and computation works can be
##   seen in the 'model.jags' file.
##
##4. Perform diagnostic analysis for m and n.
##
##5. Create a single summary plot containing the plot of actual daily death,
##   posterior mean of expected daily death, posterior mean of daily new 
##   infections, and the 95% credible interval for daily new infections which
##   observed in day of the year.
##------------------------------------------------------------------------------

##Load library to perform posterior sampling using JAGS

library(rjags)

library(coda)

##Data on daily hospital deaths with COVID-19 in England for the first 100 days 
##from March 2nd 2020,from NHS England.

y <-c(1,2,0,2,2,0,4,4,1,9,14,20,22,27,40,46,66,63,105,103,149,159,204,263,326,
      353,360,437,498,576,645,647,700,778,743,727,813,900,792,740,779,718,699,
      646,686,639,610,571,524,566,486,502,451,438,386,380,345,341,325,313,306,
      268,251,259,252,266,256,215,203,196,166,183,162,180,172,167,138,160,144,
      153,150,122,128,116,133,139,122,124,117,92,83,94,109,111,83,86,83,80,73,69)

##Recall that there were no COVID-19 deaths prior to March 2nd 2020. So for 
##modelling purposes it makes sense to extend the data by adding 20 days of zero
##deaths to the start of the data in order to observe the data from February 11, 
##2020.

y_new <- c(rep(0,20), y) ##Daily deaths data from February 11, 2020

##Create matrix B
B <- matrix(0, nrow = length(y_new), ncol = length(y_new))

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

mod <- jags.model("model.jags", data = list(y = y_new, N = length(y_new), B = B))

sam.coda <- coda.samples(mod, c("m","n"), n.iter = 10000)


##Diagnostic plots (trace plot and autocorrelation plots)

par(mfrow = c(4,2), mar = c(2,2,2,2))

traceplot(sam.coda[[1]][,c('n[2]','n[50]','n[60]','n[100]','m[2]','m[50]','m[60]','m[120]')])

acfplot(sam.coda[[1]][,c('n[2]','n[50]','n[60]','n[100]','m[2]','m[50]','m[60]','m[120]')], aspect = 1, type = 'l')


##From the trace plot we can see that the chain for m and n have a good mixing 
##near its peaks compare near the beginning and the end. In addition, we also 
##found out that the chain autocorrelation reduces faster near its peaks 
##compare to near the beginning and the end. 


##Recall that the effective sample size is compute daily because we have 120 
##m and n respectively. For the recommended sample size, we take the maximum
##value of the effective sample size of all m and n without distinct them 
##because it would be difficult to work with distinct sample size and sometime
##too small suitable for some nodes in m and n.


##Recommended sample for iterations

recommended_sample <- max(effectiveSize(sam.coda[[1]]))
recommended_sample

##After we run the code with different choice of m and n, we find a burn-in 
##period from the first iterations until around the 4000 iterations. So, we 
##recommend to do a burn-in for 4000 sample in the first iteration for the 
##simulations. Because the sampling is quite random, we would suggest that 
##the effective sample size should be around 400 until 800 sample. 


##Regarding the recommended number of sampling iteration to run,we 
##found out that when we run for 30000 iterations we get the effective sample 
##size is not too small and the autocorrelation reduces faster compare to using
##10000 iterations. Thus, we would recommend to run 30000 sampling iterations.


##Find sample for expected death and new infections from the sample list

death_index <- grep("m", colnames(sam.coda[[1]]))

infection_index <- grep("n", colnames(sam.coda[[1]]))[1:length(y)]


##Get credible interval for new infections along with its bounds

credible_interval <- apply(sam.coda[[1]][,infection_index], 2, quantile, prob = (c(0.025,0.975)))

upper_bound <- credible_interval[2,] 

lower_bound <- credible_interval[1,]


##Posterior mean for new infections and expected death

expect_death <- colMeans(sam.coda[[1]][,death_index])

new_infect <- colMeans(sam.coda[[1]][,infection_index])


##In order to avoid more than one label in x axis when plotting using a time 
##period and to give the axis a meaningful interpretation, we will use day of
##the year, which we can get by using a Julian days (days since some origin).
##Let day 1 be January 1, 2020. Because Julian function in R count its origin 
##day as day 0, then we set the origin date as December 31, 2019. 


##Find the Julian day for the first UK lockdown (March 24, 2020)

lockdown_day <- julian(as.Date("2020-3-24"), origin = as.Date("2019-12-31")) 


##Generate Julian day for observing daily death and new infections

date_vector <- seq(as.Date("2020-2-11"), by = "days", length.out = length(y_new))


julian_day <- julian(date_vector, origin = as.Date("2019-12-31")) 


infection_day <- julian_day[1:(length(y_new)-20)] ##New infections observation days


##In order to get better visualization of credible interval plot, it would be
##better if we could visualize the interval region. We could use Polygon and
##rev function in R to create a closed region of the interval with shading
##color to its region, which reference of how to use it can be seen in this 
##following links: 
#1. https://stackoverflow.com/questions/14069629/how-can-i-plot-data-with-confidence-intervals
#2. https://statisticsglobe.com/r-polygon-function-plot/
#3. https://statisticsglobe.com/rev-r-function


##Plot actual daily death

par(mfrow = c(1,1), mar = c(5.1,4.1,4.1,2.1))

max_vertical <- max(max(new_infect),max(expect_death),max(upper_bound),max(lower_bound),max(y_new))

plot(x = julian_day, y = y_new, xlab = 'Day of the year', ylab = 'Number of Individuals', ylim = c(0,max_vertical), col = 'grey', pch = 20)


##Draw credible interval region along with its boundary

polygon(x = c(rev(infection_day), infection_day), y = c(rev(lower_bound),upper_bound), col = 'gray91', border = NA) 

lines(x = infection_day, y = lower_bound, col = 'grey', lty = 2) 

lines(x = infection_day, y = upper_bound, col = 'grey', lty = 2) 


##Plot posterior mean for new infections and expected death

lines(x = infection_day, y = new_infect, col = 'green') 

lines(x = julian_day, y = expect_death, col = 'blue') 


##Highlight the first day of UK lockdown (24th March 2020)

abline(v = lockdown_day,lty = 2) 

text(x = lockdown_day+4, y = 1100, label = 'lockdown I', pos= 3, cex = 0.6, srt = 90) 


##Add legend and title to the single summary plot

text_legend <- c("Mean for New Infections", "Mean for Expected Death", "95% Credible Interval (New Infections)","Actual Daily Death")

legend(x = 'topright', legend = text_legend, col = c("green", "blue",NA,'grey'), lty = c(1,1,NA,NA), fill = c(NA,NA, 'gray91',NA),
       border = c(NA,NA,'gray91',NA), pch = c(NA,NA,NA,20),cex = 0.55,bty = 'n', x.intersp = c(2,2,1.5,2.2))

title('Daily Death and New Infections from COVID-19 at England (2020)', cex.main = 0.8) 

##Conclusion:
##From the single summary plot we can see that the posterior mean for the 
##expected death fits well with the daily death data from NHS. In addition, 
##we can see that the posterior mean for new infections lies inside the 95%
##credible interval region for new infections. According to the final plot, we
##can see that the estimated mean for daily new infections already decline 
##before lockdown started and almost remained constant after the lockdown 
##started, but the timing of lockdown appears to coincide when the expected 
##deaths rises before reaching its peaks.