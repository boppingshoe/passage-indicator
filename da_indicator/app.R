
# this is a Shiny web app for monitoring LGS adult Chinooka passage

library(shiny)
library(arm)
load(file= 'data/ad_dat.Rdata')
load(file= 'data/pitflow2.Rdata')
options(scipen=999) # keep plot from displaying scientific notation

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

tt_func<- function(dat, year, strt, cutoff, nsim, use_median='t'){
  mmdl<- lmer(I(1/ftt)~ jday+ dis_ihr+ km+ ihr_temp+ (1|mig_his)+ (1|yr),
    data=subset(dat, !yr %in% c(2011, year)) )
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
  sub_yr<- subset(da_yr, jday>=st & jday<=en)
  Observed<- getCI(sub_yr$ftt)
  # Observed<- getCI(subset(da_yr, jday>=strt & jday<=cutoff)$ftt)
  sumtab<- cbind(c("Min.", "2.5%", "Median","Mean", "97.5%", "Max."),
    round(Predicted, 2), round(Observed, 2))
  sumtab<- rbind(sumtab, cbind('n =',nrow(sub_yr),' '))
  colnames(sumtab)<- c(' ','Predicted','Observed')
  outtie<- list()
  outtie$sumtab<- sumtab
  outtie$ms<- ms
  return(outtie)
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
        plotOutput("passage_plot"),
        
        fluidRow(
          column(4, tableOutput("travel_time")),
          column(8, plotOutput("ftt_hist"))
        )
      )
   )
)

# Define server logic required to display a plot for expected counts -----

server <- function(input, output) {
  pi_out<- reactive({
    input$reload # reload simulation
    isolate({
      po<- list()
      po$da_yr<- input$year
      po$nsim<- input$nsim
      return(po)
    })
  })
  
  output$passage_plot <- renderPlot({
    
    # n_sim<- input$nsim
    n_sim<- pi_out()$nsim
    ran<- input$raneff
    d91<- input$dat91
    if (d91==TRUE) {
      dat<- subset(ad_dat, lgs>0)
    } else {dat<- subset(ad_dat, obs_yr>=2009 & lgs>0)}
    da_yr<- pi_out()$da_yr
    
    # fitting a linear regression model and simulate predicted coefs
    if (ran==TRUE) {
      mdl1<- lmer(log(lgs)~ log(lmn_1)+ (log(lmn_1)|obs_yr), 
        data= subset(dat, lmn_1>100 & obs_yr!=da_yr))
    } else {
      mdl1<- lm(log(lgs)~ log(lmn_1), data= subset(dat, lmn_1>100 & obs_yr!=da_yr))
    }
    
    sim1<- sim(mdl1, n.sims=n_sim)

    # plotting fish counts
    with(ad_dat[ad_dat$obs_yr==da_yr,], 
      plot(obs_date, cumsum(lgs), pch= 20, main= da_yr, ylim= c(0, sum(lgs)+10000),
        xlim= as.Date(c(paste0(da_yr,"-04-15"), paste0(da_yr,"-06-30"))),
        xaxt='n', xlab= NA, ylab= 'Cumulative Counts at LGS', ty='n'))
    mrk<- seq(as.Date(paste0(da_yr,"-04-15")), as.Date(paste0(da_yr,"-07-01")), "week")
    axis(side= 1, at= mrk, labels= substr(mrk, 6,10))
    # plot predicted lgs counts as uncertainties
    for(i in 1:n_sim){
      if (ran==TRUE) {
        with(subset(ad_dat, obs_yr==da_yr & lmn_1>0),
          lines(obs_date, cumsum(exp(sim1@fixef[i, 1]+ sim1@fixef[i, 2]*log(lmn_1)) ),
            col='grey50'))
      } else{
        with(subset(ad_dat, obs_yr==da_yr & lmn_1>0),
          lines(obs_date, cumsum(exp(sim1@coef[i, 1]+ sim1@coef[i, 2]*log(lmn_1)) ),
            col='grey50'))
      }
    }

         with(subset(ad_dat, obs_yr==da_yr), # observed lmn counts
       lines(obs_date, cumsum(lmn_1), pch=20, lwd=3, col='blue'))
     with(subset(ad_dat, obs_yr==da_yr), # observed lgs counts
       lines(obs_date, cumsum(lgs), pch=20, lwd=3, col='red'))
     with(subset(ad_dat, obs_yr==da_yr),
       legend(min(obs_date)+12, sum(lgs)+7500,
         c('LMN counts','LGS counts','Expected LGS counts'),
         lty=1, lwd=c(3,3,10), cex=1.2,
         col=c('blue','red','grey50'), bty='n'))
   })
  
  outtie<- reactive({
    n_sim<- pi_out()$nsim
    da_yr<- pi_out()$da_yr
    strt<- format(input$strt, format='%m-%d')
    cutoff<- format(input$cutoff, format='%m-%d')
    outtie<- tt_func(pitflow2, da_yr, strt, cutoff, n_sim)
  })
  
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



















