

suppressPackageStartupMessages({
  library(shiny)
  library(tidyverse)
  library(plotly)
  library(leaflet)
  library(shinydashboard)
})

data <- read_csv("Output/NRS_CombinedWaterQuality.csv") %>% pivot_longer(-c(NRScode:IMOSsampleCode)) %>% drop_na() 

mapdata <- data %>% group_by(NRScode, Station, Longitude, Latitude) %>% summarise(count = n())

server <- function(input, output) {
  output$plot <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%  
      addCircleMarkers(
        lng = mapdata$Longitude,
        lat = mapdata$Latitude, 
        radius = mapdata$count/100, 
        label = mapdata$Station, 
        weight = 2)
  })
}

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
    menuItem("maps", tabName = "maps", icon = icon("th"))
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "dashboard",
            h2("Some info"),
  fluidRow(
    valueBox(
      value = nrow(mapdata), 
      subtitle = "Total sampling trips", 
      icon = icon("anchor")
    ),
    valueBox(
      value = nrow(data), 
      subtitle = "Total samples analysed",
      icon = icon("microscope")
    ),
  )),
  tabItem(tabName = "maps",
  fluidRow(
    h2("Sample Locations"),
    leafletOutput("plot")
  ))
)
)

ui <- dashboardPage(
  skin = "green",
  header = dashboardHeader(title = "BGC data parameters"),
  sidebar = sidebar,
  body = body
)

shinyApp(ui, server)
