
# This file contains analysis for lgs passage indicator

rm(list=ls())
wd<- "G:/STAFF/Bobby/adult passage"

#
# simulate expected lgs counts using observed lmn counts ----
load(file= paste0(wd, '/data compile/ad_dat.Rdata'))
load(file= paste0(wd, '/data compile/adop.Rdata'))

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
# plotting dist/ftt
for(i in unique(pitflow2$yr)){
  hist(157/(subset(pitflow2, yr==i)$ftt), breaks=500, xlab='Velocity (km/day)', main=i)
  print(summary(subset(pitflow2, yr==i)$ftt))
}

# mixed models
require(arm)
year<- 2018
# looking at model manually ----
# inverse gaussian is much slower in comupting, but mathematically more legit
mmdl1<- glmer(ftt~ sjday+ sdis_ihr+ skm+ sihr_temp+ mig_his+ (1|yr),
  data=subset(pitflow2, !yr %in% c(2011, year)),
  # family= Gamma(link='identity') )
  family=inverse.gaussian(link='identity') )
mmdl1<- lmer(I(157/ftt)~ sjday+ sihr_temp+ sdis_ihr+ skm
  + mig_his+ (1|yr), data=subset(pitflow2, !yr %in% c(2011, year)) )
#
vif.lme <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v }
vif.lme(mmdl1)

# m_res<- residuals(mmdl1)
y<- mmdl1@resp$y
V<- mmdl1@pp$`.->V`
H <- V %*% solve(t(V) %*% V) %*% t(V)
hat_val <- diag(H)
m_fit<- fitted(mmdl1)

m_res<- (y+ 2*(m_fit^-0.5))/ sqrt(1-hat_val) # constant info scale
s<- sd(y- m_fit)
m_res<- (y- m_fit)/ (s*sqrt(1-hat_val))
# m_dr<- mmdl1@resp$devResid()
plot(m_fit, m_res) # not sure resid plot is helpful
plot(mmdl1)
summary(mmdl1)

plot(mmdl2)
summary(mmdl2)
# ----

getCI <- function(ms) {
  result <- rep(NA, 6)
  result[1] <- min(ms, na.rm = TRUE)
  result[2] <- quantile(ms, 0.025, na.rm = TRUE)
  result[3] <- median(ms, na.rm = TRUE)
  result[4] <- mean(ms, na.rm = TRUE)
  result[5] <- quantile(ms, 0.975, na.rm = TRUE)
  result[6] <- max(ms, na.rm = TRUE)
  return(result)
}
tt_func<- function(dat, betties, year, strt, cutoff, realtime_adj='f'){
  # specify start and cutoff dates
  st<- as.numeric(format(as.Date(paste0(year, '-', strt)), format='%j'))
  en<- as.numeric(format(as.Date(paste0(year, '-', cutoff)), format='%j'))
  
  # subset data based on cutoff dates
  subdat<- subset(dat, yr==year & jday>=st & jday<=en)
  subdat[subdat$jday2>en|is.na(subdat$jday2), c('jday2','ftt')]<- NA
  
  # simulate mortality
  subdat$surv<- rbinom(nrow(subdat), 1, 0.96)
  subdat2<- subset(subdat, surv==1)
  
  # for inverse gaussian
  vmtx<- cbind(1, subset(subdat2,
    select=c('sjday','sdis_ihr','skm','sihr_temp','mig_his') ))
  vmtx$mig_his<- ifelse(vmtx$mig_his=='trans', 1, 0) 
  fttsim<- betties %*% t(vmtx)

  # conversion rate
  jday2sim<- t(t(fttsim)+ subdat2$jday)
  con_all<- apply(jday2sim, 1, function(x) sum(ceiling(x)<=en)/nrow(subdat) )
  con_pre<- cbind(getCI(con_all))
  row.names(con_pre) <- c("Min.", "2.5%", "Median", "Mean", "97.5%", "Max.")
  # travel time
  if(realtime_adj=='f'){
    Observed<- getCI(subdat$ftt)
  } else {
    subadj<- subset(dat, yr==year & jday>=st & jday<=en)
    Observed<- getCI(subadj$ftt)
  }
  ms<- apply(fttsim, 1, median)
  Predicted<- getCI(ms)
  
  
  sumtab<- cbind(Predicted, Observed)
  row.names(sumtab) <- c("Min.", "2.5%", "Median", "Mean", "97.5%", "Max.")
  # pit-tag counts
  ihr_ct<- table(subdat$jday)
  gra_ct<- with(subdat, tapply(jday2, jday, function(x) sum(!is.na(x))))
  gra_pre<- apply(jday2sim, 1,
    function(x) tapply(x, subdat2$jday, function(x) sum(ceiling(x)<=en) ))
  pit_ct<- cbind(ihr_ct, gra_ct, gra_pre)
  
  output<- list()
  output$n<- nrow(subdat)
  output$sumtab<- sumtab
  output$ms<- ms
  output$con_obs<- mean(!is.na(subdat$jday2))
  output$con_all<- con_all
  output$con_pre<- con_pre
  output$pit_ct<- pit_ct
  return(output)
}
#
# snap shot ----
year<- 2017
nsim<- 200
# specify model
mmdl<- lmer(I(1/ftt)~ jday+ dis_ihr+ km+ ihr_temp+ mig_his+ (1|yr),
  data=subset(pitflow2, !yr %in% c(2011, year)) )
mmdl<- glmer(ftt~ sjday+ sdis_ihr+ skm+ sihr_temp+ mig_his+ (1|yr),
  data=subset(pitflow2, !yr %in% c(2011, year)),
  family= Gamma('identity') )
mmdl<- glmer(ftt~ sjday+ sdis_ihr+ skm+ sihr_temp+ mig_his+ (1|yr),
  data=subset(pitflow2, !yr %in% c(2011, year)),
  family=inverse.gaussian(link='identity') )
betties<- sim(mmdl, n.sims=nsim)@fixef # simulation

out<- tt_func(pitflow2, betties, year, strt='04-01', cutoff='06-30')
out$sumtab
cat('n = ', out$n)
cat('Expected conversion rate = '); out$con_pre
cat('Observed conversion rate = ', out$con_obs)
obs_med<- out$sumtab[3,2]
# histogram for ftt median
hist(out$ms, breaks=100,
  xlim= c(min(min(out$ms),obs_med)-1, max(max(out$ms),obs_med)+1), main=year)
abline(v=obs_med, col='red', lwd=2)
# cumulative counts using pit-tags
plot(cumsum(out$pit_ct[,1]), pch=20, main=year, ylim=c(0, sum(out$pit_ct[,1])+100), xlab=NA, ylab='Cumulative PIT-tag Counts', ty='n')
for(s in 1:nsim){
    lines(cumsum(out$pit_ct[,s+2]), col='grey50')
}
lines(cumsum(out$pit_ct[,1]), pch=20, lwd=2, col='blue')
lines(cumsum(out$pit_ct[,2]), pch=20, lwd=2, col='red')
# ----

# summary look (for individual year)
prep_it<- function(dat, year, nsim, allDates, inverse_gau='t', adj){
  # specify model
  if(inverse_gau=='f'){
    mmdl<- glmer(ftt~ sjday+ sdis_ihr+ skm+ sihr_temp+ mig_his+ (1|yr),
      data=subset(dat, !yr %in% c(2011, year)),
      family=Gamma(link='identity') )
  } else{
    mmdl<- glmer(ftt~ sjday+ sdis_ihr+ skm+ sihr_temp+ mig_his+ (1|yr),
      data=subset(dat, !yr %in% c(2011, year)),
      family=inverse.gaussian(link='identity') )
  }
  # betties<- sim(mmdl, n.sims=nsim)@fixef # simulation
  coefs<- fixef(mmdl)
  vcdf<- as.data.frame(VarCorr(mmdl))
  varesp<- vcdf[vcdf$grp=='Residual', 'vcov']
  betties<- as.matrix(cbind(statmod::rinvgauss(nsim, mean=coefs[1], dispersion=varesp),
    coefs[2], coefs[3], coefs[4], coefs[5], coefs[6]))
  #
  prep_out<- list()
  prep_out$ftt_pre<- matrix(NA, nrow=3, ncol=length(allDates))
  # prep_out$ftt_pre<- matrix(NA, nrow=nsim, ncol=length(allDates))
  prep_out$ftt_obs<- rep(NA, length(allDates))
  prep_out$conv_pre<- matrix(NA, nrow=3, ncol=length(allDates))
  # prep_out$conv_pre<- matrix(NA, nrow=nsim, ncol=length(allDates))
  prep_out$conv_obs<- rep(NA, length(allDates))
  
  for(i in 1:length(allDates)){
    out<- tt_func(dat, betties, year, strt='04-01', cutoff=allDates[i],
      realtime_adj=adj)
    # travel time
    prep_out$ftt_pre[,i]<- out$sumtab[c(1,3,6),1]
    # prep_out$ftt_pre[,i]<- out$ms
    prep_out$ftt_obs[i]<- out$sumtab[3,2]
    # conversion
    prep_out$conv_pre[,i]<- out$con_pre[c(1,4,6),1]
    # prep_out$conv_pre[,i]<- out$con_all
    prep_out$conv_obs[i]<- out$con_obs
  }
  return(prep_out)
}
plot_it<- function(a, b, ftt_obs, ftt_pre, conv_obs, conv_pre){
  # ftt
  plot(a,0, xlim=c(a,b),
    ylim=c(min(min(ftt_obs, na.rm=TRUE), min(ftt_pre, na.rm=TRUE))-0.5,
      max(max(ftt_obs, na.rm=TRUE), max(ftt_pre, na.rm=TRUE))+0.5),
    main=year, xlab=NA, ylab='Fish Travel Time (day)', ty='n')
  # apply(ftt_pre, 1, function(x) lines(seq(a, b, 'day'), x, col='grey30'))
  # lines(seq(a, b, 'day'), apply(ftt_pre, 2, median), lty=2, col='grey80')
  lines(seq(a, b, 'day'), ftt_pre[1,], lty=2, lwd=2)
  lines(seq(a, b, 'day'), ftt_pre[2,], lty=1, lwd=2)
  lines(seq(a, b, 'day'), ftt_pre[3,], lty=2, lwd=2)
  lines(seq(a, b, 'day'), ftt_obs, lwd=3, col='red')
  # conv
  plot(a,0, xlim=c(a,b),
    ylim=c(min(min(conv_obs, na.rm=TRUE), min(conv_pre, na.rm=TRUE))-0.05, 1),
    main=year, xlab=NA, ylab='Conversion', ty='n')
  # apply(conv_pre, 1, function(x) lines(seq(a, b, 'day'), x, col='grey30'))
  # lines(seq(a, b, 'day'), apply(conv_pre, 2, mean), lty=2, col='grey80')
  lines(seq(a, b, 'day'), conv_pre[1,], lty=2, lwd=2)
  lines(seq(a, b, 'day'), conv_pre[2,], lty=1, lwd=2)
  lines(seq(a, b, 'day'), conv_pre[3,], lty=2, lwd=2)
  lines(seq(a, b, 'day'), conv_obs, lwd=3, col='coral')
}
#
# summary look (for individual year) ----
year<-2017; start<-'05-10'; end<-'06-30'
a<- as.Date(paste0(year, '-', start))
b<-as.Date(paste0(year, '-', end))
allDates <- format(seq(a, b, 'day'), format='%m-%d')
nsim<- 300

windows()
par(mfrow=c(2,1))
prep_out<- prep_it(pitflow2, year, nsim, allDates, inverse_gau='f', adj='f')
plot_it(a, b, prep_out$ftt_obs, prep_out$ftt_pre, prep_out$conv_obs, prep_out$conv_pre)
#
# summary look (for multiple years) ----
par(mfrow=c(3,4))
years2c<- 2005:2010
years2c<- 2012:2017
start<-'05-05'; end<-'06-30'
nsim<- 100
for(year in years2c){
  a<- as.Date(paste0(year, '-', start))
  b<-as.Date(paste0(year, '-', end))
  allDates <- format(seq(a, b, 'day'), format='%m-%d')
  
  prep_out<- prep_it(pitflow2, year, nsim, allDates, inverse_gau='t', adj='t')
  plot_it(a, b, prep_out$ftt_obs, prep_out$ftt_pre,
    prep_out$conv_obs, prep_out$conv_pre)
}

# quick look on conversion ----
pitflow2$surv<- ifelse(is.na(pitflow2$ftt), 0, 1)
smd<- glm(surv~ 1, data=pitflow2, family=binomial)
summary(smd)
smmd<- glmer(surv~ 1+ (1|yr) , data=pitflow2, family=binomial)
summary(smmd)

# plotting historical ftt dist ----
# from Jerry
histhist<- function(dat, year){
  att_all <- subset(dat, !yr %in% c(2011, year))$ftt
  att_dayr <- subset(dat, yr==year)$ftt
  
  # Set the break points for histogram
  breakvals = seq(1, 21, by=1)
  labelvals = seq(1, 20, by=1)
  
  # Define bins for each data point of histogram
  att_all.cut = cut(att_all, breaks=breakvals, labels=labelvals,right=FALSE)
  att_dayr.cut = cut(att_dayr, breaks=breakvals, labels=labelvals,right=FALSE)
  
  # populate bins delineated by breaks 
  att_all.freq = table(att_all.cut)
  att_dayr.freq = table(att_dayr.cut)
  
  # Calculate relative frequency
  att_all.relfreq = att_all.freq/nrow(as.data.frame(att_all))
  att_dayr.relfreq = att_dayr.freq/nrow(as.data.frame(att_dayr))

  # use barplot to precisely control what you get...
  ymax<- max(max(att_all.relfreq),max(att_dayr.relfreq))
  barplot(att_all.relfreq, col=rgb(0,1,0,1/4), xlim=(c(1,23)), ylim=c(0,ymax),
    xlab='Travel Time (Day)', ylab='Percent', main=paste('Travel Time,', year))
  barplot(att_dayr.relfreq, col=rgb(1,0,0,1/4), xlim=(c(1,23)), add=T)
  legend(15, ymax, c('Historical FTT','Observed FTT'), pch=15,
    col=c('lightgreen','pink'), bty='n')
  legend(15, ymax, c(' ',' '), pch=22, bty='n')
}

windows()
par(mfrow=c(3,5))
for(i in unique(pitflow2$yr)){
  histhist(pitflow2, i)
}








