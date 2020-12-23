# Shiny dashboard practice - do we want to do something like this as a front page?

suppressPackageStartupMessages({
  library(sf)
  library(shiny)
  library(plotly)
  library(leaflet)
  library(purrr)
  library(lubridate)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(shinydashboard)
  library(tidyverse)
})

rawD <- "RawData"
outD <- "Output"

# Prepare the data sets to be shown in the dashboard
PCIdata <- read_csv(paste0(rawD,.Platform$file.sep,"PCISampCPR.csv"), na = c("(null)", NaN)) %>% 
  rename(route = ROUTE, region = REGION, Latitude = LATITUDE, Longitude = LONGITUDE, SampleDate = SAMPLEDATEUTC)
Zdata <- read_csv(paste0(rawD,.Platform$file.sep,"ZSampCPR.csv"), na = c("(null)", NaN)) %>% 
  rename(route = ROUTE, region = REGION, Latitude = LATITUDE, Longitude = LONGITUDE, SampleDate = SAMPLEDATEUTC)
Pdata <- read_csv(paste0(rawD,.Platform$file.sep,"PSampCPR.csv"), na = c("(null)", NaN)) %>% 
  rename(route = ROUTE, region = REGION, Latitude = LATITUDE, Longitude = LONGITUDE, SampleDate = SAMPLEDATEUTC)
miles <- read_csv(paste0(rawD,.Platform$file.sep,"mileageCPR.csv"), na = c("(null)", NaN)) %>%
  rename(tow = TOW, SampleDate = SAMPLEDATEUTC, miles = MILES) %>%
  mutate(year = year(SampleDate)) %>% group_by(year) %>% summarise(miles = sum(miles)) %>%
  arrange(year) %>% mutate(miles = cumsum(miles))

mapdata <- Zdata %>% select(c("region", "Longitude", "Latitude")) %>% # select all the distinct samples from the P and Z datasets
  merge(Pdata %>% select(c("region", "Longitude", "Latitude")), all = TRUE) 

maptext <- mapdata %>% group_by(region) %>% summarise(textno = n()) %>% # make the text markers for each region
  mutate(x = c(130, 160, 160, 160, 160,160, 160, 100),
         y = c(-42, -20, -30, -10, -35, -60, -43, -30))

# data for the PCI maps
PCImap <- PCIdata %>% mutate(Latitude = round(Latitude/0.5, 0)*0.5, 
                             Longitude = round(Longitude/0.5, 0)*0.5,
                             mon = month(SampleDate, label = TRUE),
                             monn = quarter(SampleDate)) %>%
  group_by(Latitude, Longitude, monn, region) %>% summarise(PCI = round(mean(PCI)/0.5,0)*0.5)

# parameters used in mapping   
aus <-  ne_countries(country = c('australia', 'papua new guinea', 'indonesia'), scale = "medium", returnclass = 'sf')   
cols <- c("white" ,"darkseagreen3", "chartreuse4","forestgreen")
colscale <- scale_colour_gradientn(name = '', colors = cols)
mind <- min(miles$year)
maxd <- max(miles$year)
n <- floor(sum(miles$miles)/21600)

server <- function(input, output) {
  output$plot <- renderPlotly({
   p <- ggplot(mapdata) + geom_sf(data = aus, fill = "cornsilk") + 
     geom_point(aes(x= Longitude, y = Latitude, color = region)) +
     geom_text(data = maptext, aes(x, y, label = paste0("n = ", textno), colour = region), fontface = "bold" ) +
     coord_sf() + labs(x="", y="") +
     scale_y_continuous(limits = c(NA, 0), expand = c(0,0)) +
     theme_bw() + theme(legend.position = "none")
  p
   ggplotly(p)       
  })
  output$plot1 <-  renderPlotly({
    miles %>% split(.$year) %>%
    accumulate(~bind_rows(.x, .y)) %>%
    set_names(mind:maxd) %>%
    bind_rows(.id = "frame") %>%
    plot_ly(x = ~year, y = ~miles) %>%
    add_lines(frame = ~frame, showlegend = FALSE) %>%
    layout(xaxis = list(title = ""),
           yaxis = list(title = "Nautical Miles Towed")) %>% 
    animation_slider(hide = T)
  })
  output$plot2 <-  renderPlotly({
      p <- ggplot(PCImap) + geom_sf(data = aus, fill = "cornsilk") + 
    geom_point(aes(x= Longitude, y = Latitude, color = PCI)) +
    colscale +
    facet_wrap(~monn) +
    coord_sf() + labs(x="", y="") +
    scale_y_continuous(limits = c(NA, 0), expand = c(0,0)) +
    theme_bw() + theme(panel.background = element_rect(fill = "aliceblue"),
                       legend.position = "none", 
                       strip.background = element_blank())
  ggplotly(p)
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
            h2("Fun CPR survey facts"),
  fluidRow(
    valueBox(
      value = nrow(PCIdata), 
      subtitle = "Total phytoplankton colour index samples", 
      icon = icon("microscope")
    ),
    valueBox(
      value = nrow(Pdata), 
      subtitle = "Total phytoplankton samples", 
      icon = icon("microscope")
    ),
    valueBox(
      value = nrow(Zdata), 
      subtitle = "Total zooplanktn samples",
      icon = icon("microscope")
    ),
  ),
  fluidRow(
    h2("Sampling Progress"),
    plotlyOutput("plot1")
  ),
  fluidRow(
    h2("Hit the play button to see the animation"), 
    h2(paste0("We have towed more than ", n, " times around the world - WOW !!"))
    )),
  tabItem(tabName = "maps",
  fluidRow(
    column(width = 8, height = 14, align="center",
    h2("Sample Locations"),
    plotlyOutput("plot"))
  ),
  fluidRow(
    column(width = 8, height = 14, align="center",
           h2("PCI through the seasons"),
           plotlyOutput("plot2"))
  ))
))

ui <- dashboardPage(
  skin = "green",
  header = dashboardHeader(title = "The AusCPR survey"),
  sidebar = sidebar,
  body = body
)

shinyApp(ui, server)

##############################################
P <- map_data("world", "australia") %>% filter(region == 'Australia') %>% plot_ly(x = ~long, y = ~lat) %>% add_markers(colors = "black")
P

# using plot geo, doesn't allow selection of Australian portion of map only
p <- PCImap %>% plot_ly(x = ~Longitude, y = ~Latitude, color = ~PCI, colors = c("#EDF8E9", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#005A32"),
                        frame = ~monn, type = 'scattergeo', mode = "markers") %>%
  layout(geo = list(scope="world"),
         xaxis = list(range = c(110, 160)),
         yaxis = list(range = c(-60, -10)))
p

# make PCImap an sf obeject so that it can be the same data set as aus
PCImapsf <- st_as_sf(PCImap, coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
test <- aus %>% st_join(PCImapsf)

p <- plot_ly(test) %>%
  add_trace(color = ~PCI, mode = 'markers', type = 'scatter')
p

p <-  PCImapsf %>% plot_geo(color = ~PCI, colors = c("#EDF8E9", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#005A32"),
                               frame = ~monn, ) 
p

%>%
  add_trace(aus, color = I("cornsilk"), stroke = I("black"), span = I(1))
p




