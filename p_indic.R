
# This file contains analysis for lgs passage indicator

rm(list=ls())
wd<- "G:/STAFF/Bobby/adult passage"
#

# simulate expected lgs counts using observed lmn counts ----
library(arm)

indc_func<- function(da_yr, n_sim=100, ran='n'){
  if (ran=='n') {
    mdl1<- lm(log(lgs)~ log(lmn_1), data= subset(adop, lmn_1>100, obs_yr!=da_yr))
  } else {
    mdl1<- lmer(log(lgs)~ log(lmn_1)+ (log(lmn_1)|obs_yr), data= subset(adop, lmn_1>100 & obs_yr!=da_yr))
  }
  
  sim1<- sim(mdl1, n.sims=n_sim)
  # colnames(sim1@coef)
  
  with(ad_dat[ad_dat$obs_yr==da_yr,], 
    plot(obs_date, cumsum(lgs), pch=20, main=da_yr, ylim=c(0, sum(lgs)+10000),
      xlab=NA, ylab='Cumulative Counts at LGS', ty='n'))
  for(i in 1:n_sim){
    if (ran=='n') {
      with(subset(ad_dat, obs_yr==da_yr & lmn_1>0),
        lines(obs_date, cumsum(exp(sim1@coef[i, 1]+ sim1@coef[i, 2]*log(lmn_1)) ),
          col='grey50'))
    } else{
      with(subset(ad_dat, obs_yr==da_yr & lmn_1>0),
        lines(obs_date, cumsum(exp(sim1@fixef[i, 1]+ sim1@fixef[i, 2]*log(lmn_1)) ),
          col='grey50'))
    }
  }
  with(ad_dat[ad_dat$obs_yr==da_yr,],
    lines(obs_date, cumsum(lmn_1), pch=20, lwd=2, col='blue'))
  with(ad_dat[ad_dat$obs_yr==da_yr,],
    lines(obs_date, cumsum(lgs), pch=20, lwd=2, col='red'))
}

windows()
par(mfrow=c(3,3))
for(yr in 1991:1999){
  indc_func(yr, ran='y')
}


# using gls
library(MASS)

indc_ar<- function(da_yr, n_sim=100, ran='n'){
  sig<- rep(NA, n_sim)
  betty<- matrix(NA, nrow=n_sim, ncol=2)
  
  if(ran=='n'){
    mdl2<- gls(log(lgs)~ log(lmn_1), correlation= corAR1(), data= subset(adop, lmn_1>100 & obs_yr!=da_yr))
  } else {
    mdl2<- lme(log(lgs)~ log(lmn_1), random= ~1|obs_yr, correlation= corAR1(), data= subset(adop, lmn_1>100 & obs_yr!=da_yr))
  }
  
  with(ad_dat[ad_dat$obs_yr==da_yr,], 
    plot(obs_date, cumsum(lgs), pch=20, main=da_yr, ylim=c(0, sum(lgs)+10000),
      xlab=NA, ylab='Cumulative Counts at LGS', ty='n'))
  
  for(s in 1:n_sim){
    if(ran=='n'){
      sig[s]<- with(mdl2, sigma * sqrt((dims$N-dims$p)/rchisq(1, dims$N-dims$p)))
      betty[s,]<- mvrnorm(1, mdl2$coefficients, mdl2$varBeta*sig[s]^2)
    } else{
      sig[s]<- with(mdl2, sigma * sqrt((dims$N-dims$Q)/rchisq(1, dims$N-dims$Q)))
      betty[s,]<- mvrnorm(1, mdl2$coefficients$fixed, mdl2$varFix*sig[s]^2)
    }
    
    with(subset(ad_dat, obs_yr==da_yr & lmn_1>0),
      lines(obs_date, cumsum(exp(betty[s, 1]+ betty[s, 2]*log(lmn_1)) ),
        col='grey50'))
  }
  
  with(ad_dat[ad_dat$obs_yr==da_yr,],
    lines(obs_date, cumsum(lmn_1), pch=20, lwd=2, col='blue'))
  with(ad_dat[ad_dat$obs_yr==da_yr,],
    lines(obs_date, cumsum(lgs), pch=20, lwd=2, col='red'))
}

for(yr in 2009:2017){
  indc_ar(yr)
}

# ----

# ratio estimate, using LMN counts to predict LGS counts ----
ratio_est<- function(x, y, mux= NA, N= NA){
  if(length(x)!= length(y)) stop('x and y must have same length')
  n<- length(x)
  fpc<- 1
  if(!is.na(N)){
    fpc<- (N-n)/N
  }
  r<- sum(y)/sum(x)
  sr2<- (1/(n-1))* sum((y- r*x)^2)
  if(is.na(mux)){
    mx<- mean(x)
  } else {mx<- mux}
  cat('r=', r, ' SE=', sqrt((fpc*sr2)/n), '\n')
  if(!is.na(N) & !is.na(mux)){
    cat('tau-hat=', N* r* mux, 'SE=', N* sqrt((fpc*sr2)/n), '\n')
  }
}
ratio_est<- function(x, y){
  if(length(x)!= length(y)) stop('x and y must have same length')
  n<- length(x)
  fpc<- 1
  r<- sum(y)/sum(x)
  sr2<- (1/(n-1))* sum((y- r*x)^2)
  mx<- mean(x)
  cat('r=', r, ' SE=', sqrt((fpc*sr2)/n), '\n')
} # simplified version

ratio_est(adop$lmn_1, adop$lgs)
# ----

# travel time and flow and stuff
# using lmn to lgs (2014- 2017) ----
# loading and format data
load(file= paste0(wd, '/data compile/pit_flow.Rdata'))

# there seems to be more (in portion) longer ftt in 2017
# longer ftt in 2017 creates a bimodel distribution in data
pit_flow$ftt<- pit_flow$ftt_l2g
windows()
par(mfrow=c(2,2))
for(i in unique(pit_flow$yr)){
  hist(subset(pit_flow, yr==i)$ftt, breaks=500, main=i,
    xlim=c(0, max(pit_flow$ftt, na.rm=TRUE)), ylim=c(0, 120))
  print(summary(subset(pit_flow, yr==i)$ftt))
}
for(i in unique(pit_flow$yr)){
  hist(1/subset(pit_flow, yr==i)$ftt, breaks=50, main=i)
  print(summary(subset(pit_flow, yr==i)$ftt))
} # 1/ftt

# what outliers to exclude?
hist(pit_flow$ftt, breaks=100)
quantile(pit_flow$ftt, c(90, 95, 97.5)/100, na.rm=TRUE)
mean(pit_flow$ftt>600, na.rm=TRUE) # <1%

# linear model
year<- 2018
fmdl1<- lm(log(ftt)~ jday+ dis_lgs+ km+ yr, data=subset(pit_flow, ftt<600 & yr!=year))
fmdl1<- lm(I(ftt^(-0.95))~ jday+ dis_lgs+ km+ yr, data=subset(pit_flow, ftt<600 & yr!=year))
summary(fmdl1)
par(mfrow=c(2,2))
plot(fmdl1)

# fmdl2<- glm(ftd~ jday+ dis_lgs+ +km+ yr, data=subset(pit_flow, ftt<600 & yr!=year), family=quasipoisson) # nice but ftt is not count
# summary(fmdl2)

# might consider Box-Cox transformation
MASS::boxcox(fmdl1, lambda = seq(-1.1, -0.8, 0.05)) # lambda=-0.95 without yr 2017

# higher flow lower ftt? (can't really tell)
par(mfrow=c(2,2))
for(i in unique(pit_flow$yr)){
  with(subset(pit_flow, yr==i&ftt<600), 
    plot(dis_lgs, log(ftt), pch=20, main=i) )
}

for(i in unique(pit_flow$yr)){
  with(subset(pit_flow, yr==i&ftt<600), 
    plot(jday, dis_lgs, pch=20, main=i) )
} # seasonality in flow

for(i in unique(pit_flow$yr)){
  with(subset(pit_flow, yr==i&ftt<600), 
    plot(jday, log(ftt), pch=20, main=i) )
} # don't see seasonality in ftt
meds<- tapply(pit_flow$ftt, pit_flow$jday, function(x) median(x, na.rm=TRUE))
plot(as.numeric(dimnames(meds)[[1]]), as.numeric(meds), pch=20) # here you see seasonalisy, but also coincide with flow
ms<- tapply(pit_flow$ftt, pit_flow$jday, function(x) mean(x, na.rm=TRUE))
plot(as.numeric(dimnames(ms)[[1]]), as.numeric(ms), pch=20)

# release site and ftt
boxplot(log(ftt)~km, data=pit_flow)

# arrival distribution
for(i in unique(pit_flow$yr)){
  with(subset(pit_flow, yr==i), 
    hist(jday, breaks=50, xlim=c(90, 200), main=i) )
} 

# mixed models
require(arm)
year<- 2009
# pit_flow$jdays<- scale(pit_flow$jday)
# pit_flow$dislgss<- scale(pit_flow$dis_lgs)
# fmmdl1<- lmer(log(ftt)~ jday+ dis_lgs+ km+ (1|yr),
#   data=subset(pit_flow, !yr%in%c(year,2017) & ftt<=600))
fmmdl1<- lmer(I(1/ftt)~ jday+ dis_lgs+ km+ (1|yr),
  data=subset(pit_flow, !yr%in%c(year,2017) & ftt<=600))
plot(fmmdl1)
summary(fmmdl1)

nsim<- 100
mdlsim<- sim(fmmdl1, n.sims=nsim)
betties<- mdlsim@fixef
da_yr<- subset(pit_flow, yr==year)
fttsim<- matrix(NA, nrow=nsim, ncol=nrow(da_yr))
for(i in 1:nsim){
  for(j in 1:nrow(da_yr)){
  fttsim[i,j]<- 1/(betties[i,1]+ 
      betties[i,2]*da_yr$jday[j]+ 
      betties[i,3]*da_yr$dis_lgs[j]+
      betties[i,4]*da_yr$km[j])
  }
}

meds<- apply(fttsim, 1, median)
quantile(meds, c(.025, .975))
summary(meds)
summary(subset(da_yr, jday<=190)$ftt)

require(nlme)
year<- 2017
# fmmdl1<- lme(log(ftt)~ jday+ dis_lgs+ km, random= ~1|yr, method='REML', data=subset(pit_flow, ftt<600 & yr!=year))
# fmmdl2<- gls(log(ftt)~ jday+ dis_lgs, correlation= corAR1(), method='REML', data=subset(pit_flow, ftt<600 & yr!=year))
# AIC(fmmdl1, fmmdl2) # corAR1 is not better 
fmmdl1<- lme(I(1/ftt)~ jday+ dis_lgs+ km, random= ~1|yr, method='REML', data=subset(pit_flow, ftt<600 & yr!=year))

windows()
par(mfrow=c(2,1))
acf(resid(fmmdl1, type='normalized'))
pacf(resid(fmmdl1, type='normalized'))
car::vif(fmmdl1)
plot(fmmdl1, pch=20); par(mfrow=c(1,1))
qqnorm(resid(fmmdl1, type='normalized'))
qqline(resid(fmmdl1, type='normalized'))
hist(resid(fmmdl1, type='normalized'), breaks=50, freq=F); curve(dnorm(x), add=TRUE)

summary(fmmdl1)

betty<- coef(fmdl1)

# using data for ice to granite (2005-2017) ----
load(file= paste0(wd, '/data compile/pitflow2.Rdata'))

windows()
par(mfrow=c(3,5))
for(i in unique(pitflow2$yr)){
  hist(1/(subset(pitflow2, yr==i)$ftt), breaks=500, main=i)
  print(summary(subset(pit_flow, yr==i)$ftt))
}

# mixed models
require(arm)
year<- 2016
fmmdl1<- lmer(I(1/ftt)~ jday+ dis_ihr+ km+ ihr_temp+ (1|mig_his)+ (1|yr),
  data=subset(pitflow2, !yr %in% c(2011, year)) )
plot(fmmdl1)
summary(fmmdl1)

getCI <- function(ms) {
  result <- rep(NA, 6)
  result[1] <- min(ms, na.rm = TRUE)
  result[2] <- quantile(ms, 0.025)
  result[3] <- median(ms, na.rm = TRUE)
  result[4] <- mean(ms, na.rm = TRUE)
  result[5] <- quantile(ms, 0.975)
  result[6] <- max(ms, na.rm = TRUE)
  return(result)
}

tt_func<- function(dat, year, strt='04-01', cutoff='06-30', nsim=100, use_median='t', mdl_summary='f'){
  mmdl<- lmer(I(1/ftt)~ jday+ dis_ihr+ km+ ihr_temp+ (1|mig_his)+ (1|yr),
    data=subset(dat, !yr %in% c(2011, year)) )
  if(mdl_summary=='t') {print(summary(mmdl)); cat('\n')}
  betties<- sim(mmdl, n.sims=nsim)@fixef
  da_yr<- subset(dat, yr==year)
  vmtx<- cbind(1, da_yr[,c('jday','dis_ihr','km','ihr_temp')])
  fttsim<- 1/(betties %*% t(vmtx))

  if(use_median=='t'){
    ms<- apply(fttsim, 1, median)
  } else {
    ms<- apply(fttsim, 1, mean)
  }
  
  Predicted<- getCI(ms)
  st<- as.numeric(format(as.Date(paste0(year, '-', strt)), format='%j'))
  en<- as.numeric(format(as.Date(paste0(year, '-', cutoff)), format='%j'))
  Observed<- getCI(subset(da_yr, jday>=st & jday<=en)$ftt)
  sumtab<- cbind(Predicted, Observed)
  row.names(sumtab) <- c("Min.", "2.5%", "Median","Mean", "97.5%", "Max.")
  output<- list()
  output$sumtab<- sumtab
  output$ms<- ms
  return(output)
}

year<- 2016
out<- tt_func(pitflow2, year, strt='04-01', cutoff='06-30', nsim=500)
out$sumtab
hist(out$ms, breaks=100, xlim= c(out$sumtab[1,2]+1, out$sumtab[4,2]))
abline(v=out$sumtab[3,2], col='red', lwd=2)

windows()
par(mfrow=c(3,5))
for (i in 2005:2017){
  out<- tt_func(pitflow2, year=i, strt='04-01', cutoff='06-30', nsim=500)
  obs_med<- out$sumtab[3,2]
  hist(out$ms, breaks=100,
    xlim= c(min(min(out$ms),obs_med)-1, max(max(out$ms),obs_med)+1), main=i)
  abline(v=obs_med, col='red', lwd=2)
}









