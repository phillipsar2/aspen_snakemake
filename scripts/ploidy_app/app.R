## shiny app ---

library(shiny)
library(plotly)
library(ggplot2)
library(sf)
library(raster)
library(ggspatial)
library(USAboundaries)
library(purrr)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)


# Load data
winner_num <- read.csv("flow_cyt_predictions.2026-08-07.csv")

north_america <- ne_countries(continent = "North America", returnclass = "sf")
states <- us_states()

meta <- read.csv("megametadata.2026-08-07.csv")
meta$seqname <- gsub(x = meta$RQC.Seq.Unit.Name, pattern = ".fastq.gz", replacement = "")
winner_meta <- merge(winner_num, meta, by.x = "sample", by.y = "seqname")

# --- Data setup (runs once when app starts) ---
north_america <- ne_countries(continent = "North America", returnclass = "sf")
states <- us_states()

events_sf <- winner_meta %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

ploidy_levels <- c("diploid", "triploid", "unknown")

# --- UI ---
ui <- fluidPage(
  titlePanel("Ploidy Map"),
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput(
        "ploidy_filter", "Show ploidy:",
        choices = ploidy_levels,
        selected = ploidy_levels
      ),
      width = 3
    ),
    mainPanel(
      plotlyOutput("ploidy_map", height = "700px"),
      width = 9
    )
  )
)

# --- Server ---
server <- function(input, output, session) {
  
  filtered_data <- reactive({
    req(input$ploidy_filter)
    events_sf %>% filter(ploidy_call %in% input$ploidy_filter)
  })
  
  output$ploidy_map <- renderPlotly({
    ploidy_map <- ggplot() + 
      geom_sf(data = north_america, size = 1, color = "black", fill = NA) +
      geom_sf(data = states, size = 1, fill = NA, color = "black") +
      geom_sf(
        data = filtered_data(), size = 3,
        aes(color = ploidy_call, shape = ploidy_call,
            text = paste0("Sample: ", Sample.Name, "<br>Ploidy: ", ploidy_call)),
        alpha = 0.7
      ) +
      theme_bw() +
      theme(panel.grid = element_blank(), 
            axis.line = element_blank()) + 
      ylim(c(10, 70)) +
      xlim(c(165, 57)) +
      scale_color_manual(
        values = c("diploid" = "blue", "triploid" = "red", "unknown" = "black"),
        drop = FALSE
      ) +
      scale_shape_manual(
        values = c("diploid" = 16, "triploid" = 17, "unknown" = 16),
        drop = FALSE
      ) +
      labs(color = "Ploidy") +
      guides(shape = "none")
    
    ggplotly(ploidy_map, tooltip = "text")
  })
}

shinyApp(ui, server)
