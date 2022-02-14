





library(shiny)
library(waiter)
library(shinythemes)
library(tidyverse)
library(DT)
library(ggsci)
library(cowplot)

side_width <- 6

# Define UI for application that draws a histogram
ui <- navbarPage(
  "DopaDB",
  selected = "SN/VTA Markers",
  theme = shinythemes::themeSelector(),
  navbarMenu(
    "Spatial Transcriptomics",
    "----",
    "Dopaminergic Neurons",
    tabPanel("Primary Markers"),
    tabPanel(
      "SN/VTA Markers",
      sidebarLayout(
        position = "left",
        sidebarPanel(
          width = side_width,
          h3(helpText("Click on a gene to view results...")),
          hr(),
          DT::dataTableOutput("SN_VTA_MARKERS"),
          hr(),
          br(),
          numericInput(
            "padj",
            label = "Filter Adj. P Value",
            value = 0.01,
            min = 0,
            max = 1
          ),
          hr(),
          br(),
          numericInput(
            "lfc",
            label = "Filter Log2 Fold-change",
            value = 0,
            min = NA,
            max = NA
          ),
          hr(),
          br(),
          numericInput(
            "count_threshold",
            label = "Count Display Threshold (Log2)",
            value = 0.5,
            min = NA,
            max = NA
          )
        ),
        mainPanel(
          use_waiter(),
          width = 12 - side_width,
          wellPanel(plotOutput("SN_VTA_VIOLIN_PLOT")),
          wellPanel(plotOutput("SN_VTA_SPATIAL_PLOT")),
          h4(helpText("Debug")),
          wellPanel(verbatimTextOutput("debug"))
        )
      )
    ),
    tabPanel("Ageing"),
    tabPanel("SNCA-OVX"),
    "----",
    "Other cell types",
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  waiter_show()
  
  # SN_VTA
  SN_VTA = reactiveValues(
    METADATA = read_csv("input/markers/sn_vs_vta_da_metadata.csv"),
    COUNTS = read_csv(
      "input/markers/sn_vs_vta_da_counts.csv",
      col_names = c("sample_cell_id", "Gene Symbol", "SCT_count")
    )
  )
  
  SN_VTA_RESULTS <- reactive({
    read_csv(
      "input/markers/sn_vs_vta_da.csv",
      col_names = c("Gene Symbol",
                    "Adjusted P Value",
                    "Log2 Fold-change"),
      skip = 1
    ) %>%
      mutate(across(c(`Adjusted P Value`, `Log2 Fold-change`), ~ signif(.x, 3))) %>%
      filter(abs(`Log2 Fold-change`) > input$lfc) %>%
      filter(`Adjusted P Value` < input$padj)
  })
  
  # SN_VTA TABLE
  output$SN_VTA_MARKERS <- DT::renderDataTable({
    SN_VTA_RESULTS()
  },
  selection = "single",
  server = FALSE)
  
  observeEvent(input$SN_VTA_MARKERS_rows_selected, {
    SN_VTA$PLOT_DATA <- SN_VTA$COUNTS %>%
      filter(`Gene Symbol` == SN_VTA_RESULTS()[input$SN_VTA_MARKERS_rows_selected, ]$`Gene Symbol`) %>%
      inner_join(SN_VTA$METADATA) %>%
      mutate(SCT_count = as.double(SCT_count))
    
    # SN_VTA VIOLIN PLOT
    output$SN_VTA_VIOLIN_PLOT <- renderPlot({
      SN_VTA$PLOT_DATA %>%
        ggplot(aes(x = region,
                   y = SCT_count,
                   fill = region)) +
        geom_violin() +
        geom_jitter(width = 0.1) +
        scale_fill_d3() +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        theme_cowplot() +
        theme(legend.position = "none") +
        labs(x = "Region",
             title = paste0(unique(SN_VTA_RESULTS()[input$SN_VTA_MARKERS_rows_selected, ]$`Gene Symbol`)))
    })
    
    # SN_VTA SPATIAL PLOT
    output$SN_VTA_SPATIAL_PLOT <- renderPlot({
      SN_VTA$PLOT_DATA %>%
        # filter(SCT_count > 0) %>%
        ggplot(aes(
          x = x,
          y = y,
          colour = SCT_count > input$count_threshold
        )) +
        geom_point() +
        facet_wrap(vars(mouse_id),
                   scales = "free") +
        scale_y_reverse() +
        scale_color_d3() +
        theme_cowplot() +
        theme(
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          legend.position = "top"
        ) +
        panel_border()
    })
  })
  
  output$debug <- renderPrint(input$SN_VTA_MARKERS_rows_selected)
  
  waiter_hide()
  
}

# Run the application
shinyApp(ui = ui, server = server)
