
# this is a Shiny web app for monitoring LGS adult Chinooka passage

library(shiny)
library(arm)
library(statmod)
library(HDInterval)
load(file= 'data/ad_dat.Rdata')
load(file= 'data/pitflow2.Rdata')
options(scipen=999) # keep plot from displaying scientific notation

# summary stats
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
# travel time function
tt_func<- function(dat, betties, year, strt, cutoff){

  # specify start and cutoff dates
  st<- as.numeric(format(as.Date(paste0(year, '-', strt)), format='%j'))
  en<- as.numeric(format(as.Date(paste0(year, '-', cutoff)), format='%j'))
  # subset data based on cutoff dates
  subdat<- subset(dat, yr==year & jday>=st & jday<=en)
  subdat[subdat$jday2>en|is.na(subdat$jday2), c('jday2','ftt')]<- NA
  
  # simulate for mortality
  subdat$surv<- rbinom(nrow(subdat), 1, 0.96)
  subdat2<- subset(subdat, surv==1)
  
  vmtx<- cbind(1, subset(subdat2,
    select=c('sjday','sdis_ihr','skm','sihr_temp','mig_his') ))
  vmtx$mig_his<- ifelse(vmtx$mig_his=='trans', 1, 0)
  # fttsim<- betties %*% t(vmtx) # for inv gauss
  fttsim<- 157/(betties %*% t(vmtx))

  # conversion rate
  jday2sim<- t(t(fttsim)+ subdat2$jday)
  conall<- apply(jday2sim, 1, function(x) sum(ceiling(x)<=en)/nrow(subdat) )
  obsconv<- mean(!is.na(subdat$jday2))
  # travel time
  fttobs<- getCI(subdat$ftt)
  # fttsim<- fttsim[jday2sim<=en & fttsim>0]
  fttsim[jday2sim>en | fttsim<0]<- NA
  # # summary table for ftt, n and obs conv
  # sumtab<- cbind(c("Min.", "2.5%", "Median","Mean", "97.5%", "Max."),
  #   round(fttobs, 3), round(getCI(fttsim),3))
  # sumtab<- rbind(sumtab,
  #   cbind('n =', nrow(subdat), ' '),
  #   cbind('Obs Conv =', round(obsconv,3), ' ') )
  # colnames(sumtab)<- c(' ','Observed', 'Predicted FTT')
  # summary table for historical
  hisdat<- subset(dat, jday>=st & jday<=en)
  hisdat[hisdat$jday2>en|is.na(hisdat$jday2), c('jday2','ftt')]<- NA
  # sumhis<- cbind(c("Min.", "2.5%", "Median","Mean", "97.5%", "Max."),
  #   round(fttobs, 3), round(getCI(hisdat$ftt),3))
  # colnames(sumhis)<- c(' ','Observed', 'Historical FTT')
  
  out<- list()
  out$n<- nrow(subdat)
  # out$sumtab<- sumtab
  out$ftt<- subdat$ftt
  out$fttsim<- fttsim
  out$conall<- conall
  out$obsconv<- obsconv
  out$en<- en
  out$hisdat<- hisdat
  # out$sumhis<- sumhis
  return(out)
}
# summarize cumulative conversion from 4/1, loop for every day in specified period
prep_it<- function(dat, betties, year, nsim, allDates){
  prep_out<- list()
  prep_out$conv_pre<- matrix(NA, nrow=3, ncol=length(allDates))
  prep_out$conv_prehdi<- matrix(NA, nrow=2, ncol=length(allDates))
  prep_out$conv_obs<- rep(NA, length(allDates))
  
  for(i in 1:length(allDates)){
    out<- tt_func(dat, betties, year, strt='04-01', cutoff=allDates[i])
    # conversion
    prep_out$conv_pre[,i]<- c(quantile(out$conall,0.15),
      median(out$conall), quantile(out$conall,0.85))
    prep_out$conv_prehdi[,i]<- hdi(out$conall,0.80)
    prep_out$conv_obs[i]<- out$obsconv
  }
  return(prep_out)
}
# conv summary plots function
plot_conv<- function(a, b, year, conv_obs, conv_pre, conv_prehdi, hdiconv){
  plot(a,0, xlim=c(a,b),
    ylim=c(min(min(conv_obs, na.rm=TRUE), min(conv_pre, na.rm=TRUE))-0.05, 1.2),
    main=year, xlab=NA, ylab='Conversion', ty='n')
  
  if(hdiconv==TRUE) {
    lines(seq(a, b, 'day'), conv_prehdi[1,], col='red')
    lines(seq(a, b, 'day'), conv_prehdi[2,], col='red')
  }
  
  lines(seq(a, b, 'day'), conv_pre[1,], lty=2, lwd=1)
  lines(seq(a, b, 'day'), conv_pre[2,], lty=1, lwd=2)
  lines(seq(a, b, 'day'), conv_pre[3,], lty=2, lwd=1)
  lines(seq(a, b, 'day'), conv_obs, lwd=3, col='coral')
  legend(a, 1.2, c('Observed','Predicted Median','70% Pred Interval'),
    lty=c(1,1,2), lwd=c(3,2,1), col=c('coral',1,1), bty='n')
}
# plot historical ftt dist
histhist<- function(hisdat, year, strt, cutoff){
  att_all <- subset(hisdat, !yr %in% c(2011, year))$ftt
  att_dayr <- subset(hisdat, yr==year)$ftt
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
    xlab='Day', ylab='Percent', main=paste('Distribution of Travel Time,', year))
  barplot(att_dayr.relfreq, col=rgb(1,0,0,1/4), xlim=(c(1,23)), add=T)
  legend(15, ymax, c('Observed','Historical'), pch=15,
    col=c('pink','lightgreen'), bty='n')
  legend(15, ymax, c(' ',' '), pch=22, bty='n')
}

# Define UI for application that...
# make sliders and buttons ----
ui <- fluidPage(
   
  # Application title
   titlePanel("LGS Passage Indicator"),
   
  # Sidebar with a slider input for which year to predict
  sidebarLayout(
    column(3,
      numericInput("year", 
        label = "Select a year to display:", 
        value = 2017, min = 2005, max = 2018, step = 1),
      numericInput("nsim", "Number of simulations:", 
        value = 300, min = 100, max = 500, step = 100),
      actionButton(inputId = "reload", label = "Simulate!"),
      sliderInput("strt", label="Start Date:", 
        min = as.Date("2017-04-01","%Y-%m-%d"),
        max = as.Date("2017-06-30","%Y-%m-%d"),
        value = as.Date("2017-05-10"), timeFormat="%m/%d", step = 1),
      sliderInput("cutoff", label="Cutoff Date:", 
        min = as.Date("2017-04-01","%Y-%m-%d"),
        max = as.Date("2017-06-30","%Y-%m-%d"),
        value = as.Date("2017-06-30"), timeFormat="%m/%d", step = 1),
      checkboxInput("hdiconv", "80% HDI for Predicted Conversion", FALSE),
      checkboxInput("historic", "Historical Data", FALSE)
    ),
  
      # Show a plot of expected passage
      mainPanel(
        fluidRow(
          column(6, plotOutput("passage_plot")),
          column(6, plotOutput("pittag_counts"))
        ),
        
        # fluidRow(
        #   column(5, tableOutput("travel_time")),
        #   column(7, plotOutput("ftt_hist"))
        # ),
        
        # fluidRow(
        #   column(6, tableOutput("his_tabl")),
        #   column(6, plotOutput("histhist"))
        # )
        fluidRow(
          column(6, plotOutput("histhist")),
          column(6, plotOutput("ftt_hist"))
        )
      )
   )
)

# Define server logic required to display plots and table ----

server <- function(input, output) {
  pi_out<- reactive({
    input$reload # reload simulation
    isolate({
      po<- list()
      po$da_yr<- input$year
      po$nsim<- input$nsim
      po$cutoff<- format(input$cutoff, format='%m-%d')
      return(po)
    })
  })

  # travel time model calculate here ----
  tt_betties<- reactive({
    da_yr<- pi_out()$da_yr
    n_sim<- pi_out()$nsim
    
    mmdl<- lmer(I(157/ftt)~ sjday+ sdis_ihr+ skm+ sihr_temp+ mig_his+ (1|yr),
      data=subset(pitflow2, !yr %in% c(2011, da_yr)) )
    # betties<- sim(mmdl, n.sims=n_sim)@fixef # simulation using arm package
    coefs<- fixef(mmdl)
    vcdf<- as.data.frame(VarCorr(mmdl))
    sdesp<- vcdf[vcdf$grp=='Residual', 'sdcor']
    betties<- as.matrix(cbind(rnorm(n_sim, mean=coefs[1], sd=sdesp),
      coefs[2], coefs[3], coefs[4], coefs[5], coefs[6]))
    # # or simulation based on inv gauss
    # mmdl<- glmer(ftt~ sjday+ sdis_ihr+ skm+ sihr_temp+ mig_his+ (1|yr),
    #   data=subset(pitflow2, !yr %in% c(2011, da_yr)),
    #   family=inverse.gaussian(link='identity') )
    # coefs<- fixef(mmdl)
    # vcdf<- as.data.frame(VarCorr(mmdl))
    # varesp<- vcdf[vcdf$grp=='Residual', 'vcov']
    # betties<- as.matrix(cbind(rinvgauss(n_sim, mean=coefs[1], dispersion=varesp),
    #   coefs[2], coefs[3], coefs[4], coefs[5], coefs[6]))
    return(betties)
  })
  
  outtie<- reactive({ # output from tt_func
    da_yr<- pi_out()$da_yr
    betties<- tt_betties()
    strt<- format(input$strt, format='%m-%d')
    cutoff<- format(input$cutoff, format='%m-%d')
    tt_func(pitflow2, betties, da_yr, strt, cutoff)
  })
  
  # Display a plot for expected counts ----

  ct_mdl<- reactive({
    da_yr<- pi_out()$da_yr
    n_sim<- pi_out()$nsim
    
    dat<- subset(ad_dat, obs_yr>=2009 & lgs>0)
    mdl1<- lmer(log(lgs)~ log(lmn_1)+ (log(lmn_1)|obs_yr), 
      data= subset(dat, lmn_1>100 & obs_yr!=da_yr))
    sim1<- sim(mdl1, n.sims=n_sim)
    mdlout<- list()
    mdlout$da_yr<- da_yr
    mdlout$n_sim<- n_sim
    mdlout$sim1<- sim1
    return(mdlout)
  })

  output$passage_plot <- renderPlot({
    da_yr<- ct_mdl()$da_yr
    n_sim<- ct_mdl()$n_sim

    coefs<- ct_mdl()$sim1@fixef
    # plotting fish counts
    plyr<- subset(ad_dat, 
      obs_yr==da_yr & as.numeric(format(obs_date, format='%j'))<=outtie()$en)
    with(subset(ad_dat, obs_yr==da_yr), 
      plot(obs_date, cumsum(lgs), pch= 20, main= da_yr, ylim= c(0, sum(lgs)+10000),
        xlim= as.Date(c(paste0(da_yr,"-04-15"), paste0(da_yr,"-06-30"))),
        xaxt='n', xlab= NA, ylab= 'Cumulative Counts', ty='n'))
    mrk<- seq(as.Date(paste0(da_yr,"-04-15")), as.Date(paste0(da_yr,"-07-01")), "week")
    axis(side= 1, at= mrk, labels= substr(mrk, 6,10))
    # plot predicted lgs counts as uncertainties
    for(i in 1:n_sim){
      with(subset(plyr, lmn_1>0),
        lines(obs_date, cumsum(exp(coefs[i, 1]+ coefs[i, 2]*log(lmn_1)) ),
          col='grey50'))
    }
    
    with(subset(plyr, lmn_1>0), # observed lmn counts
      lines(obs_date, cumsum(lmn_1), pch=20, lwd=3, col='blue'))
    with(subset(plyr, lmn_1>0), # observed lgs counts
      lines(obs_date, cumsum(lgs), pch=20, lwd=3, col='red'))
    with(subset(ad_dat, obs_yr==da_yr),
      legend(min(obs_date)+12, sum(lgs)+7500,
        c('LMN counts','LGS counts','Expected LGS counts'),
        lty=1, lwd=c(3,3,10), cex=1.2,
        col=c('blue','red','grey50'), bty='n'))
  })
  
  # Display a plot for PIT-tag conversion -----  
  
  output$pittag_counts <- renderPlot({
    da_yr<- pi_out()$da_yr
    n_sim<- pi_out()$nsim
    betties<- tt_betties()
    strt<- format(input$strt, format='%m-%d')
    cutoff<- format(input$cutoff, format='%m-%d')
    a<- as.Date(paste0(da_yr, '-', strt))
    b<-as.Date(paste0(da_yr, '-', cutoff))
    allDates <- format(seq(a, b, 'day'), format='%m-%d')
    
    prep_out<- prep_it(pitflow2, betties, da_yr, n_sim, allDates)
    plot_conv(a, b, da_yr, prep_out$conv_obs, prep_out$conv_pre, prep_out$conv_prehdi, input$hdiconv)
  })
  
  # Display travel time summary table and a histogram -----

  # output$travel_time <- renderTable({
  #   if(pi_out()$da_yr>2004) {outtie()$sumtab}
  # })
  
  output$ftt_hist <- renderPlot({
    da_yr<- pi_out()$da_yr
    if(da_yr>2004) {
      ftt<- outtie()$ftt # define travel time (obs)
      ftt[is.na(ftt)]<- -10

      if(input$historic==FALSE){
        fttsim<- outtie()$fttsim # ftt (simulated)
        fttsim[is.na(fttsim)]<- -10
        xmax<- max(max(ftt), max(fttsim))
        fshis<- hist(fttsim, breaks=seq(-10, ceiling(xmax), 1), plot=FALSE)
        ftthis<- hist(ftt, breaks=seq(-10, ceiling(xmax), 1), plot=FALSE)
        ymax<- max(max(fshis$density[-1]),max(ftthis$density[-1]))
        bt<- 'Predicted'
      } else{
        hisftt<- outtie()$hisdat$ftt # historical ftt
        hisftt[is.na(hisftt)]<- -10
        xmax<- max(max(ftt), max(hisftt))
        fshis<- hist(hisftt, breaks=seq(-10, ceiling(xmax), 1), plot=FALSE)
        ftthis<- hist(ftt, breaks=seq(-10, ceiling(xmax), 1), plot=FALSE)
        ymax<- max(max(fshis$density[-1]),max(ftthis$density[-1]))
        bt<- 'Historic'
      }
      plot(fshis, freq=FALSE, col=rgb(0,0,0,.8),
        xlab='Day', main= paste('Travel Time (IHR-LGR),', da_yr),
        xlim=c(0, xmax+1), ylim=c(0, ymax*1.2))
      hist(ftt, breaks=seq(-10, ceiling(xmax), 1), freq=FALSE,
        col=rgb(1,1,1,.5), add=TRUE)
      legend(xmax*.6, ymax*1.3,
        c( paste('n=', outtie()$n),
          paste('Observed=',round(sum(ftthis$density[-1]),3)),
          paste0(bt,'= ',round(sum(fshis$density[-1]),3)) ),
        pch=c(1,0,15), col=c('white',1,1), bty='n')
    }
  })
  
  # output$his_tabl <- renderTable({
  #   if(pi_out()$da_yr>2004) {outtie()$sumhis}
  # })
  
  output$histhist <- renderPlot({
    da_yr<- pi_out()$da_yr
    if(da_yr>2004) {
      strt<- format(input$strt, format='%m-%d')
      cutoff<- format(input$cutoff, format='%m-%d')
      hisdat<- outtie()$hisdat
      
      histhist(hisdat, da_yr, strt, cutoff)
    }
  })
}

# Run the application ----
shinyApp(ui = ui, server = server)



















