
# this is a Shiny web app for monitoring LGS adult Chinooka passage

library(shiny)
library(arm)
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
  ms<- apply(fttsim, 1, median)

  # travel time
  Predicted<- getCI(ms)
  Observed<- getCI(subdat$ftt)
  # conversion rate
  jday2sim<- t(t(fttsim)+ subdat2$jday)
  con_all<- apply(jday2sim, 1, function(x) sum(ceiling(x)<=en)/nrow(subdat) )
  ExpectedConv<- getCI(con_all)
  sumtab<- cbind(c("Min.", "2.5%", "Median","Mean", "97.5%", "Max."),
    round(Predicted, 3), round(Observed, 3), round(ExpectedConv, 3))
  sumtab<- rbind(sumtab,
    cbind('n =', nrow(subdat), 'Obs Conv =', round(mean(!is.na(subdat$jday2)),3)) )
  colnames(sumtab)<- c(' ', 'Predicted FTT', 'Observed FTT', 'Predicted Conv')
  # pit-tag counts
  ihr_ct<- table(subdat$jday)
  gra_ct<- with(subdat, tapply(jday2, jday, function(x) sum(!is.na(x))))
  gra_pre<- apply(jday2sim, 1, 
    function(x) tapply(x, subdat2$jday, function(x) sum(ceiling(x)<=en)) )
  pit_ct<- cbind(ihr_ct, gra_ct, gra_pre)
  
  out<- list()
  out$n<- nrow(subdat)
  out$sumtab<- sumtab
  out$ms<- ms
  out$con_all<- con_all
  out$pit_ct<- pit_ct
  out$en<- en
  return(out)
}
# summary function
prep_it<- function(dat, betties, year, nsim, allDates){
  prep_out<- list()
  prep_out$ftt_pre<- matrix(NA, nrow=3, ncol=length(allDates))
  # prep_out$ftt_pre<- matrix(NA, nrow=nsim, ncol=length(allDates))
  prep_out$ftt_obs<- rep(NA, length(allDates))
  prep_out$conv_pre<- matrix(NA, nrow=3, ncol=length(allDates))
  # prep_out$conv_pre<- matrix(NA, nrow=nsim, ncol=length(allDates))
  prep_out$conv_obs<- rep(NA, length(allDates))
  
  for(i in 1:length(allDates)){
    out<- tt_func(dat, betties, year, strt='04-01', cutoff=allDates[i])
    # travel time
    prep_out$ftt_pre[,i]<- as.numeric(out$sumtab[c(1,3,6),2])
    # prep_out$ftt_pre[,i]<- out$ms
    prep_out$ftt_obs[i]<- as.numeric(out$sumtab[3,3])
    # conversion
    prep_out$conv_pre[,i]<- as.numeric(out$sumtab[c(1,4,6),4])
    # prep_out$conv_pre[,i]<- out$con_all
    prep_out$conv_obs[i]<- as.numeric(out$sumtab[7,4])
  }
  return(prep_out)
}
# conv summary plots function
plot_conv<- function(a, b, year, conv_obs, conv_pre){
  # conv
  plot(a,0, xlim=c(a,b),
    ylim=c(min(min(conv_obs, na.rm=TRUE), min(conv_pre, na.rm=TRUE))-0.05, 1),
    main=year, xlab=NA, ylab='Conversion', ty='n')
  # apply(conv_pre, 1, function(x) lines(seq(a, b, 'day'), x, col='grey30'))
  # lines(seq(a, b, 'day'), apply(conv_pre, 2, mean), lty=2, col='grey80')
  lines(seq(a, b, 'day'), conv_pre[1,], lty=2, lwd=1)
  lines(seq(a, b, 'day'), conv_pre[2,], lty=1, lwd=2)
  lines(seq(a, b, 'day'), conv_pre[3,], lty=2, lwd=1)
  lines(seq(a, b, 'day'), conv_obs, lwd=3, col='coral')
}
# ftt summary plot function
plot_ftt<- function(a, b, year, ftt_obs, ftt_pre){
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
        value = 2017, min = 1991, max = 2018, step = 1),
      numericInput("nsim", "Number of simulations:", 
        value = 100, min = 100, max = 500, step = 100),
      actionButton(inputId = "reload", label = "Simulate!"),
      checkboxInput("dat91", "Use data since 1991", value = FALSE),
      checkboxInput("raneff", "Random year effects", value = TRUE),
      sliderInput("strt", label="Start Date:", 
        min = as.Date("2017-04-01","%Y-%m-%d"),
        max = as.Date("2017-06-30","%Y-%m-%d"),
        value = as.Date("2017-05-10"), timeFormat="%m/%d", step = 1),
      sliderInput("cutoff", label="Cutoff Date:", 
        min = as.Date("2017-04-01","%Y-%m-%d"),
        max = as.Date("2017-06-30","%Y-%m-%d"),
        value = as.Date("2017-06-30"), timeFormat="%m/%d", step = 1)
    ),
  
      # Show a plot of expected passage
      mainPanel(
        fluidRow(
          column(6, plotOutput("passage_plot")),
          column(6, plotOutput("pittag_counts"))
        ),
        
        fluidRow(
          column(6, tableOutput("travel_time")),
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
    
    mmdl<- glmer(ftt~ sjday+ sdis_ihr+ skm+ sihr_temp+ mig_his+ (1|yr),
      data=subset(pitflow2, !yr %in% c(2011, da_yr)),
      family=inverse.gaussian(link='identity') )
    betties<- sim(mmdl, n.sims=n_sim)@fixef # simulation
    return(betties)
  })
  
  outtie<- reactive({ # output from tt_func
    n_sim<- pi_out()$nsim
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
    
    ran<- input$raneff
    d91<- input$dat91
    # fitting a linear regression model and simulate predicted coefs
    if (d91==TRUE) {
      dat<- subset(ad_dat, lgs>0)
    } else {dat<- subset(ad_dat, obs_yr>=2009 & lgs>0)}
    
    if (ran==TRUE) {
      mdl1<- lmer(log(lgs)~ log(lmn_1)+ (log(lmn_1)|obs_yr), 
        data= subset(dat, lmn_1>100 & obs_yr!=da_yr))
    } else {
      mdl1<- lm(log(lgs)~ log(lmn_1), data= subset(dat, lmn_1>100 & obs_yr!=da_yr))
    }
    sim1<- sim(mdl1, n.sims=n_sim)
    mdlout<- list()
    mdlout$da_yr<- da_yr
    mdlout$n_sim<- n_sim
    mdlout$sim1<- sim1
    return(mdlout)
  })

  output$passage_plot <- renderPlot({
    ran<- input$raneff
    if (ran==TRUE) {
      coefs<- ct_mdl()$sim1@fixef
    } else {
      coefs<- ct_mdl()$sim1@coef
    }
    da_yr<- ct_mdl()$da_yr
    n_sim<- ct_mdl()$n_sim
    ran<- input$raneff

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
    n_sim<- pi_out()$n_sim
    betties<- tt_betties()
    strt<- format(input$strt, format='%m-%d')
    cutoff<- format(input$cutoff, format='%m-%d')
    a<- as.Date(paste0(da_yr, '-', strt))
    b<-as.Date(paste0(da_yr, '-', cutoff))
    allDates <- format(seq(a, b, 'day'), format='%m-%d')
    
    prep_out<- prep_it(pitflow2, betties, da_yr, n_sim, allDates)
    plot_conv(a, b, da_yr, prep_out$conv_obs, prep_out$conv_pre)
    legend(as.Date(paste0(da_yr,'-06-10')), 0.3, c('Observed','Predicted'),
      lty=1, lwd=c(3,2), col=c('coral',1), bty='n')
  })
  
  # Display travel time summary table and a histogram for distribution of median -----

  output$travel_time <- renderTable({
    if(pi_out()$da_yr>2004) {outtie()$sumtab}
  })
  
  output$ftt_hist <- renderPlot({
    if(pi_out()$da_yr>2004) {
      ms<- outtie()$ms
      obs_med<- as.numeric(outtie()$sumtab[3,3])
      hist(ms, breaks=50,
        xlim= c(min(min(outtie()$ms),obs_med)-1, max(max(outtie()$ms),obs_med)+1),
        main= paste('Distribution of Median Travel Time,', pi_out()$da_yr),
        xlab= NULL)
      abline(v=obs_med, col='red', lwd=2)
      text(obs_med, 3, labels = 'Observed Median')
    }
  })

}

# Run the application ----
shinyApp(ui = ui, server = server)



















