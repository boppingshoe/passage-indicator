
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

# travel time and flow and stuff ----
load(file= paste0(wd, '/data compile/pit_flow.Rdata'))
pit_flow$ftt<- as.numeric(pit_flow$ftt)
pit_flow$ftm<- pit_flow$ftt*60 # travel time in min
pit_flow$flow<- ifelse(pit_flow$dis_lgs>100, '3high',
  ifelse(pit_flow$dis_lgs<60, '1low', '2med'))

hist(pit_flow$ftt, breaks=50)
hist(subset(pit_flow, srrt=='11H')$ftt, col='red', breaks=50)
hist(subset(pit_flow, srrt=='11W')$ftt, col='blue', breaks=20, add=TRUE)

abline(v= quantile(pit_flow$ftm, c(90, 95, 97.5)/100))
# 90%       95%       97.5% 
# 5520.992  8089.155 11507.052 
hist(pit_flow[pit_flow$ftm<4500,'ftm'], breaks=50)
hist(subset(pit_flow, srrt=='11H'&ftm<4500)$ftm, col='red', breaks=50, add=TRUE)
hist(subset(pit_flow, srrt=='11W'&ftm<4500)$ftm, col='blue', breaks=50, add=TRUE)

pf2<- subset(pit_flow, ftm>2000&ftm<3800)
round(table(pf2$tag_site)/ table(pit_flow$tag_site), 2)
summary(pf2$jday) ; summary(pit_flow$jday)
pf2d<- ifelse(pit_flow$ftm<2000, 0,
  ifelse(pit_flow$ftm>3800, 0, 1))
boxplot(pit_flow$jday~pf2d)
table(pit_flow$srrt, pf2d)

windows()
par(mfrow=c(2,2))
for(i in unique(pit_flow$yr)){
  hist(subset(pit_flow, yr==i)$ftt, breaks=500, main=i,
    xlim=c(0, max(pit_flow$ftt)), ylim=c(0, 120))
  print(summary(subset(pit_flow, yr==i)$ftt))
}
tbl1<- with(pit_flow, table(yr, pf2d))
tbl1[,2]/rowSums(tbl1)

tftt<- pit_flow$ftt^(-0.9)
fmdl1<- lm(tftt~ jday+ dis_lgs, data=pit_flow)
summary(fmdl1)
windows()
par(mfrow=c(2,2))
plot(fmdl1)

MASS::boxcox(fmdl1, lambda = seq(-1, -.8, 0.05)) # lambda=-0.9

require(lme4)
year<- 2017
fmmdl1<- lmer(log(ftt)~ jday+ dis_lgs+ (1|yr),
  data=subset(pit_flow, ftt<600 & yr!=year))
plot(fmmdl1)
summary(fmmdl1)

betty<- fixef(fmmdl1)
exp(betty[1]+ betty[2]*150+ betty[3]*80)


require(nlme)
year<- 2018
fmmdl1<- gls(log(ftt)~ jday+ dis_lgs, correlation= corAR1(), method='REML', data=subset(pit_flow, ftt<600 & yr!=year))
fmmdl2<- lme(log(ftt)~ jday+ dis_lgs, random= ~1|yr, method='REML', data=subset(pit_flow, ftt<600 & yr!=year))

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
exp(betty[1]+ betty[2]*140+ betty[3]*55) # low flow
exp(betty[1]+ betty[2]*140+ betty[3]*75+ betty[4]) # med flow
exp(betty[1]+ betty[2]*140+ betty[3]*110+ betty[5]) # high flow









