#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(EpiModel)
source("Shiny_SEIQHRFNetModules.R")
library(plotly)
library(tidyverse)
library(shinythemes)


# Define UI for application that draws a histogram
ui <- fluidPage(
    
    theme = shinytheme("slate"),
    
    # Application title
    titlePanel(h1("AI4GOOD - Network Simulator",
                  h4("Luis Chaves' work for AI4GOOD",
                     h5("In the AI4GOOD simulator project we are using statistical models
                        to simulate the spread of COVID19 in refugee camps so that we action can
                        be taken in these sites with minimal resources. In this simple app, the user
                        is able to simulate different very basic scenarios. We are running similar
                        but much more complicated models in the background to find the best solution
                        as soon as possible.")))),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fluidRow(column(12,
                            column(6,
                            h5("Population parameters"),
                            sliderInput("n",
                                        "Number of individuals:",
                                        min = 10,
                                        max = 200,
                                        value = 50),
                            sliderInput("density",
                                        h5("Population density",
                                           h6("Percentage of the population each individual is connected to")),
                                        min = 0,
                                        max = 1,
                                        value = 0.05),
                            sliderInput("days",
                                        "Number of days simulated",
                                        min = 30,
                                        max = 90, 
                                        value = 30),
                     ),
                     column(6,
                            h5("Disease parameters"),
                            sliderInput("fat.rate",
                                        "Fatality rate (Proportion of those requiring
                                        hospitalisation that die per day - on average)",
                                        min = 0,
                                        max = 1, 
                                        value = 1/50),
                            sliderInput("quar.rate",
                                        "Quarantining rate (Proportion of infected people that quarantine per day)",
                                        min = 0, 
                                        max = 1,
                                        value = 1/30)
                     )),
                     h5("This plot below can help you get an idea of what different number of
               individuals and density may look like."),
                     plotOutput("networkFeel"),
                     h5(""),
                     actionButton("StartSim",
                                  "Press to start simulation", width  = "100%"),
            )
            
        ), 
        
        # Show a plot of the generated distribution
        mainPanel(
            h3("Results"),
            plotlyOutput("resPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$networkFeel = renderPlot({
        
        plot(igraph::erdos.renyi.game(input$n,
                                      input$density*(input$n)*(input$n-1)/2,
                                      type = "gnm"),
             vertex.size = 2,
             vertex.label = NA)
    })
    
    observeEvent(input$StartSim,{
        n = isolate(input$n)
        density = isolate(input$density)
        
        progress = Progress$new()
        on.exit(progress$close())
        progress$set(message = "Initialising network", value = 0)
        nw = network.initialize(n = n, directed = FALSE)
        
        formation = ~edges
        target.stats = density*(n*(n-1)/2)
        coef.diss = dissolution_coefs(dissolution = ~offset(edges),
                                      duration = 7) 
        
        
        est <- netest(nw,
                      formation,
                      target.stats,
                      coef.diss,
                      edapprox = T,
                      verbose = F)
        
        
        param = param.net(act.rate.se = 10,
                          inf.prob.se = 0.02,
                          act.rate.si = 10,
                          inf.prob.si = 0.05,
                          act.rate.sq = 2.5,
                          inf.prob.sq = 0.02,
                          ei.rate = 1/10,
                          iq.rate = isolate(input$quar.rate),#1/30, #c(rep(1/30, 60), rep(15/30, 120)), # time varying works
                          ih.rate = 1/100,
                          qh.rate = 1/100,
                          hr.rate = 1/15,
                          qr.rate = 1/20,
                          hf.rate = isolate(input$fat.rate),
                          hf.rate.overcap = 1/25,
                          hosp.cap = 5,
                          hosp.tcoeff = 0.5,
                          vital = F
        ) 
        
        init = init.net(i.num = 3,
                        r.num = 0,
                        e.num = 0,
                        s.num = n - 3,
                        f.num = 0,
                        h.num = 0,
                        q.num = 0
        )
        
        control = control.net(
            nsims = 1, 
            nsteps = isolate(input$days),
            # delete.nodes = T,  this does not work for now
            ncores = 1,#cores,
            initialize.FUN = custom.initialize.net, # this bit is just so that I can extract time
            exposure.FUN = exposure,
            infect.FUN = infect,
            quarantine.FUN = quarantining,
            hospitalize.FUN = RequireHospitalization,
            recover.FUN = recover,
            fatality.FUN = fatality,
            recovery.FUN = NULL,
            infection.FUN = NULL,
            departures.FUN = departures.net,
            get_prev.FUN = custom.get_prev.net,
            skip_check = FALSE,
            depend = F
        )
        
        
        progress$set("Starting simulation", value = 0.3)
        sim = custom.netsim(est, param, init, control)
        res = as.data.frame(sim)
        progress$set("Plotting results", value = 0.9)
        #################################################################
        ##                    OUTPUT HANDLING BELOW                    ##
        #################################################################
        output$resPlot <- renderPlotly(ggplotly(
            res %>% 
                select(s.num, e.num, i.num, q.num,
                       h.num, r.num, f.num, num, time) %>%
                rename(Susceptible = s.num,Exposed = e.num, Infected = i.num,
                       Quarantined = q.num,  `Require hospitalization` = h.num,
                       Recovered = r.num, Fatalities = f.num, Total = num) %>% 
                group_by(time) %>%
                summarise_all(~mean(.)) %>% 
                pivot_longer(-time) %>%
                ggplot(., aes(x = time, y = value, color = name))+
                geom_line(size = 1)+
                scale_color_brewer(palette = "Set1")+
                ggtitle("Number of individuals per compartment")+
                ylab("Number of people")+
                xlab("Time (in days)")+
                labs("Compartment"))
        )
        
        
        progress$set("Done, thanks for your patience.", value = 1)
    })
}
# Run the application 
shinyApp(ui = ui, server = server)
