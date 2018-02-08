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
   titlePanel("Benner Indicator"),
   
   # Sidebar with a slider input for which year to predict
  sidebarLayout(
    column(3,
      selectInput("year", 
        label = "Select a year to display", 
        choices = list("2009" = 2009, "2010" = 2010, "2011" = 2011, "2012" = 2012,
          "2013" = 2013, "2014" = 2014, "2015" = 2015, "2016" = 2016, "2017" = 2017),
        selected = 2009)),

      # Show a plot of expected passage
      mainPanel(
         plotOutput("passage_plot")
      )
   )
)

# Define server logic required to... -----
server <- function(input, output) {
  
  # yr_out<- reactive({
  #     yr_set<- input$year
  #     return(yr_set)
  # })

    
  output$passage_plot <- renderPlot({
    da_yr<- input$year
    #da_yr<- yr_out()$yr_set
    mdl1<- lm(log(lgs)~ log(lmn_1), data= subset(adop, lmn_1>100, obs_yr!=da_yr))
    sim1<- sim(mdl1, n.sims=500)
    colnames(sim1@coef)
    
    with(ad_dat[ad_dat$obs_yr==da_yr,], 
       plot(obs_date, cumsum(lgs), pch=20, main=da_yr, ylim=c(0, sum(lgs)+10000),
         xlab=NA, ylab='Cumulative Counts at LGS', ty='n'))
     for(i in 1:500){
       with(subset(ad_dat, obs_yr==da_yr & lmn_1>0),
         lines(obs_date, cumsum(exp(sim1@coef[i, 1]+ sim1@coef[i, 2]*log(lmn_1)) ),
           col='grey50'))
     }
     # with(ad_dat[ad_dat$obs_yr==da_yr,],
     #   lines(obs_date, cumsum(lmn_1), pch=20, lwd=2, col='blue'))
     with(ad_dat[ad_dat$obs_yr==da_yr,],
       lines(obs_date, cumsum(lgs), pch=20, lwd=2, col='red'))
   })
}

# Run the application ----
shinyApp(ui = ui, server = server)

