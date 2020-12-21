setwd("C://Users/dav649/Documents/GitHub/IMOS_Toolbox/Plankton")

suppressPackageStartupMessages({
  library(shiny)
  library(tidyverse)
  library(plotly)
})

data <- read_csv("Output/NRS_CombinedWaterQuality.csv") %>% pivot_longer(-c(NRScode:IMOSsampleCode)) %>% drop_na()

ui <- fluidPage(
  titlePanel("BGC paramters exploration tool"),
  
  sidebarLayout(
    sidebarPanel(
      # station selector
      selectInput(inputId = 'station', label = "Select a station", choices = unique(data$Station), selected = 'Port Hacking', multiple = TRUE),
      # Date selector
      dateRangeInput("date", "Select a date range", start = "2009-01-01", end = "2020-11-30", min = "2009-01-01", max = Sys.Date()),
      # select parameter
      selectInput(inputId = 'parameter', label = 'Select a parameter', choices = unique(data$name), selected = 'Silicate_umol_L', multiple = TRUE),
      selectInput(inputId = 'depth', label = 'Select a depth', choices = NULL),
      # Select whether to overlay smooth trend line
      checkboxInput(inputId = "smoother", label = strong("Overlay smooth trend line"), value = FALSE),
      
    ),
    mainPanel(
      tabsetPanel(
      tabPanel("Plot", plotlyOutput("plot")),
      tabPanel("Data table", DT::DTOutput('table'))
      )
    )
  )
)

# Define server function
server <- function(session, input, output) {
  # select depths
  observe({
    updateSelectInput(session, "depth", "Select a depth", 
                      choices = data[data$Station == input$station & data$name == input$parameter,]$SampleDepth_m)
  })
  
  # Subset data
  selected <- reactive({
    req(input$date)
    validate(need(!is.na(input$date[1]) & !is.na(input$date[2]), "Error: Please provide both a start and an end date."))
    validate(need(input$date[1] < input$date[2], "Error: Start date should be earlier than end date."))
    data %>%
      filter(
        Station %in% input$station,
        SampleDateLocal > as.POSIXct(input$date[1]) & SampleDateLocal < as.POSIXct(input$date[2]),
        name %in% input$parameter, 
        SampleDepth_m == input$depth
        ) %>%
      mutate(Station = as.factor(Station),
             name = as.factor(name))
  })
  
  # Create timeseries object the plotOutput function is expecting
  output$plot <- renderPlotly({
          p <- ggplot(selected()) + geom_line(aes(SampleDateLocal, value, colour = Station)) +
          labs(x = "Time") +
          theme_bw()
        if(nlevels(unique(selected()$name)) > 1){
          p <- p + facet_grid(name~., scales = "free") 
        }
        if(input$smoother){
          p <- p + geom_smooth(aes(SampleDateLocal, value, colour = Station))
        }
        ggplotly(p) 
  })
  
  # create table output
  output$table <- DT::renderDataTable(
    selected()
  )
}

# Create Shiny object
shinyApp(ui = ui, server = server)

