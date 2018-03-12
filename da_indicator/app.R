
# this is a Shiny web app for monitoring LGS adult Chinooka passage

library(shiny)
library(arm)
load(file= 'data/ad_dat.Rdata')
load(file= 'data/pitflow2.Rdata')
options(scipen=999) # keep plot from displaying scientific notation
# pitflow2$sjday<- scale(pitflow2$jday)
# pitflow2$sdis_ihr<- scale(pitflow2$dis_ihr)
# pitflow2$skm<- scale(pitflow2$km)
# pitflow2$sihr_temp<- scale(pitflow2$ihr_temp)
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
tt_func<- function(dat, year, strt, cutoff, nsim=100, use_median='t'){
  mmdl<- lmer(I(1/ftt)~ jday+ dis_ihr+ km+ ihr_temp+ mig_his+ (1|yr),
    data=subset(dat, !yr %in% c(2011, year)) )
  # mmdl<- glmer(ftt~ sjday+ sdis_ihr+ skm+ sihr_temp+ mig_his+ (1|yr),
  #   data=subset(dat, !yr %in% c(2011, year)),
  #   family=inverse.gaussian(link='identity') ) # slow!
  betties<- sim(mmdl, n.sims=nsim)@fixef
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
    select=c('jday','dis_ihr','km','ihr_temp','mig_his') ))
  vmtx$mig_his<- ifelse(vmtx$mig_his=='trans', 1, 0)
  fttsim<- 1/(betties %*% t(vmtx))
  # for glmer inverse guassian model
  # vmtx<- cbind(1, subset(subdat2,
  #   select=c('sjday','sdis_ihr','skm','sihr_temp','mig_his') ))
  # vmtx$mig_his<- ifelse(vmtx$mig_his=='trans', 1, 0)
  # fttsim<- betties %*% t(vmtx) 
  
  if(use_median=='t'){
    ms<- apply(fttsim, 1, median)
  } else {
    ms<- apply(fttsim, 1, mean)
  }
  
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
  out$sumtab<- sumtab
  out$ms<- ms
  out$pit_ct<- pit_ct
  return(out)
}

# Define UI for application that... ----
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
        value = as.Date("2017-04-01"), timeFormat="%m/%d", step = 1),
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

# Define server logic required to display plots and table

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

  # Display a plot for expected counts -----  

  output$passage_plot <- renderPlot({
    
    da_yr<- pi_out()$da_yr
    en<- as.numeric(format(as.Date(paste0(da_yr, '-', pi_out()$cutoff)), format='%j'))
    n_sim<- pi_out()$nsim
    ran<- input$raneff
    d91<- input$dat91
    if (d91==TRUE) {
      dat<- subset(ad_dat, lgs>0)
    } else {dat<- subset(ad_dat, obs_yr>=2009 & lgs>0)}
    
    # fitting a linear regression model and simulate predicted coefs
    if (ran==TRUE) {
      mdl1<- lmer(log(lgs)~ log(lmn_1)+ (log(lmn_1)|obs_yr), 
        data= subset(dat, lmn_1>100 & obs_yr!=da_yr))
    } else {
      mdl1<- lm(log(lgs)~ log(lmn_1), data= subset(dat, lmn_1>100 & obs_yr!=da_yr))
    }
    
    sim1<- sim(mdl1, n.sims=n_sim)

    # plotting fish counts
    plyr<- subset(ad_dat, 
      obs_yr==da_yr & as.numeric(format(obs_date, format='%j'))<=en)
    with(subset(ad_dat, obs_yr==da_yr), 
      plot(obs_date, cumsum(lgs), pch= 20, main= da_yr, ylim= c(0, sum(lgs)+10000),
        xlim= as.Date(c(paste0(da_yr,"-04-15"), paste0(da_yr,"-06-30"))),
        xaxt='n', xlab= NA, ylab= 'Cumulative Counts', ty='n'))
    mrk<- seq(as.Date(paste0(da_yr,"-04-15")), as.Date(paste0(da_yr,"-07-01")), "week")
    axis(side= 1, at= mrk, labels= substr(mrk, 6,10))
    # plot predicted lgs counts as uncertainties
    for(i in 1:n_sim){
      if (ran==TRUE) {
        with(subset(plyr, lmn_1>0),
          lines(obs_date, cumsum(exp(sim1@fixef[i, 1]+ sim1@fixef[i, 2]*log(lmn_1)) ),
            col='grey50'))
      } else{
        with(subset(plyr, lmn_1>0),
          lines(obs_date, cumsum(exp(sim1@coef[i, 1]+ sim1@coef[i, 2]*log(lmn_1)) ),
            col='grey50'))
      }
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
  
  outtie<- reactive({
    n_sim<- pi_out()$nsim
    da_yr<- pi_out()$da_yr
    strt<- format(input$strt, format='%m-%d')
    cutoff<- format(input$cutoff, format='%m-%d')
    tt_func(pitflow2, da_yr, strt, cutoff, n_sim)
  })

  output$pittag_counts <- renderPlot({
    da_yr<- pi_out()$da_yr
    dates<- as.Date(as.numeric(row.names(outtie()$pit_ct)),
      origin=as.Date(paste0(da_yr,"-01-01")) )
    plot(dates, cumsum(outtie()$pit_ct[,1]), pch=20, main=da_yr,
      ylim=c(0, sum(outtie()$pit_ct[,1])+100), xaxt='n',
      xlim= c(min(dates), max(dates)),
      xlab=NA, ylab='Cumulative PIT-tag Counts', ty='n')
    mrks<- seq(min(dates), max(dates), "week")
    axis(side= 1, at= mrks, labels= substr(mrks, 6,10))
    
    for(s in 1:pi_out()$nsim){
      lines(dates, cumsum(outtie()$pit_ct[,s+2]), col='grey50', lwd=3, lty=2)
    }
    # lines(dates, cumsum(outtie()$pit_ct[,1]), pch=20, lwd=2, col='blue')
    lines(dates, cumsum(outtie()$pit_ct[,2]), pch=20, lwd=2, col='coral')
    legend(min(dates), sum(outtie()$pit_ct[,1])+100,
      c('Expected GRA PIT-tag counts','GRA PIT-tag counts'),
      lty=c(2,1), lwd=c(3,2), cex=1.2,
      col=c('grey50','coral'), bty='n')
    
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



















