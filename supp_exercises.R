#used in testing the effect of data resolution on fit for exponential and Weibull models
library(plotrix)

#functions
survivorship <- function(dur){
  out <- c()
  for(i in sort(unique(dur))){
    out <- c(out,sum(dur >= i)/length(dur))}
  return(out)}

#==empirical==
dat <- read.csv('cramptonetal2016data.csv')
dat <- dat[dat$FA_age!=dat$LA_age,] #remove 0 duration
dur_my <- dat$FA_age-dat$LA_age
unitlength = 5000
J = 18 #entire dataset
dur_lifespans <- dur_my*1E6/unitlength
dur_timesteps <- dur_lifespans*J
dur_timesteps_r <- ceiling(dur_timesteps)
resolution_y <- 50000
#resolution_timesteps <- resolution_y/unitlength*J #resolution in timesteps = resolution in years/unitlength*J
resolution_timesteps <- 1
dur_timesteps_r_coarse <- ceiling(dur_timesteps_r/resolution_timesteps)*resolution_timesteps

#==theoretical==
#params
lambdaexp <- 1/6214
lambdaweibull <- 1/4341
shapeweibull <- 0.59
#processing
delt <- sort(unique(dur_timesteps_r_coarse))
counts <- as.numeric(table(delt[match(dur_timesteps_r_coarse, delt)]))
#expprob_precumulated <- lambdaexp*exp(-lambdaexp*delt_precumulated)

#fitting exponential
# explogliks <- c();for(i in 5500:7500){
#   lambdaexp <- 1/i
#   expprobs <- (1-exp(-lambdaexp*delt)) - (1-exp(-lambdaexp*(delt-resolution_timesteps)))
#   explogliks <- c(explogliks,sum(log(expprobs)*counts))
# };plot(5500:7500,explogliks)
# (5500:7500)[which(explogliks == max(explogliks), arr.ind = TRUE)]
#expprobs <- (1-exp(-lambdaexp*delt)) - (1-exp(-lambdaexp*(delt-resolution_timesteps)))
#print(paste("exponential loglik:",sum(log(expprobs)*counts)))
#fitexp <- fitdist(dur_timesteps_r,"weibull",fix.arg = list(shape=1))
#weibprobs <- (1-exp(-(lambdaweibull*delt)^shapeweibull)) - (1-exp(-(lambdaweibull*(delt-resolution_timesteps))^shapeweibull))
#print(paste("Weibull loglik:",sum(log(weibprobs)*counts)))

#fitting weibull
# scales <- 4000:6000; shapes <- seq(0.57,0.8,0.01)
# weiblogliks <- matrix(NA,nrow=length(scales),ncol=length(shapes));for(i in seq(length(scales))){for(j in seq(length(shapes))){
#   lambdaweibull <- 1/scales[i]
#   weibprobs <- (1-exp(-(lambdaweibull*delt)^shapes[j])) - (1-exp(-(lambdaweibull*(delt-resolution_timesteps))^shapes[j]))
#   weiblogliks[i,j] <- sum(log(weibprobs)*counts)
# }};image(scales,shapes,weiblogliks,col=heat.colors(100));contour(scales,shapes,weiblogliks,add=TRUE)
# max(weiblogliks)
# indices <- which(weiblogliks == max(weiblogliks), arr.ind = TRUE); print(paste(scales[indices[1]],shapes[indices[2]]))


#==plotting==
#setup
plot(delt,survivorship(dur_timesteps_r_coarse),log='y',type='n',xlim=c(0,57600),ylim=c(1E-4,1),ylab='Survivorship',xlab='Species age (My)',xaxt='n')
axis(1,at=c(0,3600*5,3600*10,3600*15),labels=c(0,5,10,15))
#theoretical
delt_theoretical <- seq(resolution_timesteps,max(dur_timesteps_r_coarse),resolution_timesteps)
expsurv <- exp(-lambdaexp*delt_theoretical)
lines(c(0,head(rep(delt_theoretical,each=2),-1)),head(rep(c(1,expsurv),each=2),-2),col='red',lwd=3)
weibsurv <- exp(-(lambdaweibull*delt_theoretical)^shapeweibull)
lines(c(0,head(rep(delt_theoretical,each=2),-1)),head(rep(c(1,weibsurv),each=2),-2),col='orange',lwd=3)
#empirical
lines(c(0,head(rep(delt,each=2),-1)),c(1,rep(survivorship(dur_timesteps_r_coarse),each=2)[-1]),lwd=2)

#pies
pie(c(0,0,1),col=c("red","orange","blue"),init.angle=90,angle=90,labels='')
pie(c(2.48356E-22,0.858835991,0.141164009),col=c("red","orange","blue"),init.angle=90,labels='')
