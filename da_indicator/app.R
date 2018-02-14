#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(arm)
load(file= 'data/adop.Rdata')
load(file= 'data/ad_dat.Rdata')

# Define UI for application that... ----
ui <- fluidPage(
   
   # Application title
   titlePanel("LGS Passage Indicator"),
   
   # Sidebar with a slider input for which year to predict
  sidebarLayout(
    column(3,
      selectInput("year", 
        label = "Select a year to display", 
        choices = list("2009" = 2009, "2010" = 2010, "2011" = 2011, "2012" = 2012,
          "2013" = 2013, "2014" = 2014, "2015" = 2015, "2016" = 2016, "2017" = 2017),
        selected = 2009),
      numericInput("nsim", "Number of simulations", 
        value = 100, min = 100, max = 500, step = 100),
      checkboxInput("raneff", "Random year effects", value = TRUE)),

      # Show a plot of expected passage
      mainPanel(
         plotOutput("passage_plot")
      )
   )
)

# Define server logic required to display a plot for expected counts -----
server <- function(input, output) {
  
  output$passage_plot <- renderPlot({
    da_yr<- input$year
    n_sim<- input$nsim
    ran<- input$raneff
    # fitting a linear regression model and simulate predicted coefs
    if (ran==TRUE) {
      mdl1<- lmer(log(lgs)~ log(lmn_1)+ (1|obs_yr), data= subset(adop, lmn_1>100 & obs_yr!=da_yr))
    } else {
      mdl1<- lm(log(lgs)~ log(lmn_1), data= subset(adop, lmn_1>100, obs_yr!=da_yr))
    }
    
    sim1<- sim(mdl1, n.sims=n_sim)

    # plotting fish counts
    with(ad_dat[ad_dat$obs_yr==da_yr,], 
      plot(obs_date, cumsum(lgs), pch=20, main=da_yr, ylim=c(0, sum(lgs)+10000),
        xlab=NA, ylab='Cumulative Counts at LGS', ty='n'))
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
       legend(min(obs_date), sum(lgs),
         c('LMN counts','LGS counts','Expected LGS counts'), lty=1, lwd=c(3,3,10),
         col=c('blue','red','grey50'), bty='n'))
   })
}

# Run the application ----
shinyApp(ui = ui, server = server)

