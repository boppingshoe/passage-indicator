
# this is a Shiny web app for monitoring LGS adult Chinooka passage

library(shiny)
library(arm)
library(statmod)
load(file= 'data/ad_dat.Rdata')
load(file= 'data/pitflow2.Rdata')
options(scipen=999) # keep plot from displaying scientific notation

# summary function
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
  fttsim<- betties %*% t(vmtx)

  # conversion rate
  jday2sim<- t(t(fttsim)+ subdat2$jday)
  conall<- apply(jday2sim, 1, function(x) sum(ceiling(x)<=en)/nrow(subdat) )
  obsconv<- mean(!is.na(subdat$jday2))
  # travel time
  fttobs<- getCI(subdat$ftt)
  fttsim<- fttsim[jday2sim<=en]
  # summary table for ftt, n and obs conv
  sumtab<- cbind(c("Min.", "2.5%", "Median","Mean", "97.5%", "Max."),
    round(fttobs, 3), round(getCI(fttsim),3))
  sumtab<- rbind(sumtab,
    cbind('n =', nrow(subdat), ' '),
    cbind('Obs Conv =', round(obsconv,3), ' ') )
  colnames(sumtab)<- c(' ','Observed', 'Predicted FTT')

  out<- list()
  out$n<- nrow(subdat)
  out$sumtab<- sumtab
  out$ftt<- subdat$ftt
  out$fttsim<- fttsim
  out$conall<- conall
  out$obsconv<- obsconv
  out$en<- en
  return(out)
}
# summary function
prep_it<- function(dat, betties, year, nsim, allDates){
  prep_out<- list()
  prep_out$conv_pre<- matrix(NA, nrow=3, ncol=length(allDates))
  prep_out$conv_obs<- rep(NA, length(allDates))
  
  for(i in 1:length(allDates)){
    out<- tt_func(dat, betties, year, strt='04-01', cutoff=allDates[i])
    # conversion
    # prep_out$conv_pre[,i]<- c(min(out$conall), mean(out$conall), max(out$conall))
    prep_out$conv_pre[,i]<- c(quantile(out$conall,0.25),
      median(out$conall), quantile(out$conall,0.75))
    prep_out$conv_obs[i]<- out$obsconv
  }
  return(prep_out)
}
# conv summary plots function
plot_conv<- function(a, b, year, conv_obs, conv_pre){
  plot(a,0, xlim=c(a,b),
    ylim=c(min(min(conv_obs, na.rm=TRUE), min(conv_pre, na.rm=TRUE))-0.05, 1.2),
    main=year, xlab=NA, ylab='Conversion', ty='n')
  lines(seq(a, b, 'day'), conv_pre[1,], lty=2, lwd=1)
  lines(seq(a, b, 'day'), conv_pre[2,], lty=1, lwd=2)
  lines(seq(a, b, 'day'), conv_pre[3,], lty=2, lwd=1)
  lines(seq(a, b, 'day'), conv_obs, lwd=3, col='coral')
  legend(a, 1.2, c('Observed','Predicted Median','50% Pred Interval'),
    lty=c(1,1,2), lwd=c(3,2,1), col=c('coral',1,1), bty='n')
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
        value = 100, min = 100, max = 500, step = 100),
      actionButton(inputId = "reload", label = "Simulate!"),
      sliderInput("strt", label="Start Date:", 
        min = as.Date("2017-04-01","%Y-%m-%d"),
        max = as.Date("2017-06-30","%Y-%m-%d"),
        value = as.Date("2017-05-10"), timeFormat="%m/%d", step = 1),
      sliderInput("cutoff", label="Cutoff Date:", 
        min = as.Date("2017-04-01","%Y-%m-%d"),
        max = as.Date("2017-06-30","%Y-%m-%d"),
        value = as.Date("2017-06-30"), timeFormat="%m/%d", step = 1),
      checkboxInput("colors", "Green", FALSE)
    ),
  
      # Show a plot of expected passage
      mainPanel(
        fluidRow(
          column(6, plotOutput("passage_plot")),
          column(6, plotOutput("pittag_counts"))
        ),
        
        fluidRow(
          column(5, tableOutput("travel_time")),
          column(7, plotOutput("ftt_hist"))
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
    
    mmdl<- glmer(ftt~ sjday+ sdis_ihr+ skm+ sihr_temp+ mig_his+ (1|yr),
      data=subset(pitflow2, !yr %in% c(2011, da_yr)),
      family=inverse.gaussian(link='identity') )
    # betties<- sim(mmdl, n.sims=n_sim)@fixef # simulation using arm package
    # or simulation based on inv gauss
    coefs<- fixef(mmdl)
    vcdf<- as.data.frame(VarCorr(mmdl))
    varesp<- vcdf[vcdf$grp=='Residual', 'vcov']
    betties<- as.matrix(cbind(rinvgauss(n_sim, mean=coefs[1], dispersion=varesp),
      coefs[2], coefs[3], coefs[4], coefs[5], coefs[6]))
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
  
  # Display a plot for PIT-tag counts -----  
  
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
    plot_conv(a, b, da_yr, prep_out$conv_obs, prep_out$conv_pre)
  })
  
  # Display travel time summary table and a histogram -----

  output$travel_time <- renderTable({
    if(pi_out()$da_yr>2004) {outtie()$sumtab}
  })
  
  output$ftt_hist <- renderPlot({
    if(pi_out()$da_yr>2004) {
      ftt<- outtie()$ftt[!is.na(outtie()$ftt)] # define travel time (obs)
      fttsim<- outtie()$fttsim # ftt (simulated)
      # define params for plotting
      da_yr<- pi_out()$da_yr
      fsimhis<- hist(fttsim, breaks=50, plot=FALSE)
      ftthis<- hist(ftt, breaks=50, plot=FALSE)
      xmax<- max(max(ftt), max(fttsim))
      ymax<- max(max(fsimhis$density),max(ftthis$density))
      
      plot(fsimhis, freq=FALSE, col=rgb(0,0,0,.9), xlab='Day',
        main= paste('Distribution of Travel Time,', da_yr),
        xlim=c(min(min(ftt), min(fttsim)), xmax),
        ylim=c(0, ymax*1.1))
      if(input$colors==FALSE){
        hist(ftt, breaks=50, freq=FALSE, col=rgb(1,1,1,.5), add=TRUE)
        legend(xmax*0.8, ymax, c('Observed','Predicted'),
          pch=c(0,15), col=c(1,1), bty='n')
      } else{
          colors<- ifelse(ftthis$breaks<median(ftt), rgb(0,0.255,0,.7), rgb(1,1,1,.5))
          hist(ftt, breaks=50, freq=FALSE, col=colors, add=TRUE)
          legend(xmax*0.8, ymax, c('Observed','< median','Predicted'),
            pch=c(0,15,15), col=c(1,rgb(0,0.255,0,.6),1), bty='n')
      }
    }
  })

}

# Run the application ----
shinyApp(ui = ui, server = server)



















