library(shiny)
library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(plotly)

url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature22402/MediaObjects/41586_2017_BFnature22402_MOESM3_ESM.xlsx"
destfile <- "zika_data.xlsx"
if(!file.exists(destfile)) curl::curl_download(url, destfile)

col_names <- c("Position", "Allele_Ancestral", "Allele_Derived", 
               "Freq_Ancestral", "Freq_Derived", "Minor_Count", 
               "Total_Count", "Codon_Pos", "Degeneracy", 
               "Codon_Change", "AA_Change", "Protein")

df_wide <- read_excel(destfile, skip = 2, col_names = col_names)

df_derived <- df_wide %>%
  mutate(across(c(Position, Freq_Derived), as.numeric),
         Protein = as.factor(trimws(Protein))) %>%
  rename(Freq = Freq_Derived)

my_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
               "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", 
               "#FC8D62", "#8DA0CB")

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { background-color: #f8f9fa; font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }
      .justification-box { 
        background-color: #ffffff; 
        padding: 20px; 
        border-radius: 8px; 
        border-left: 5px solid #2c3e50; 
        margin-bottom: 25px; 
        box-shadow: 0 2px 4px rgba(0,0,0,0.05);
      }
      h3 { color: #2c3e50; font-weight: bold; }
      .highlight { color: #e67e22; font-weight: bold; }
    "))
  ),
  
  titlePanel("Practical Session: Zika Virus Genomic Interactivity (Part C)"),
  hr(),
  
  tabsetPanel(
    tabPanel("Genomic Distribution",
             fluidRow(
               column(12,
                      div(class = "justification-box",
                          h3("Scatter Plot: htmlwidget Integration"),
                          HTML("<p>We utilized <b>plotly</b> to address the challenge of overplotting in high-density genomic data. 
                               By converting the static ggplot to an interactive widget, we allow the identification of 
                               exact genomic positions and frequencies through tooltips. 
                               This interactivity is essential for pinpointing specific outliers across the Zika genome without visual clutter.</p>")
                      ),
                      plotlyOutput("scatterPlot", height = "500px")
               )
             )
    ),
    
    tabPanel("Mutational Load",
             sidebarLayout(
               sidebarPanel(
                 h4("Interactive Controls"),
                 selectInput("prot_select", "Select Protein Subset:", 
                             choices = levels(df_derived$Protein), 
                             selected = levels(df_derived$Protein), 
                             multiple = TRUE),
                 hr(),
                 div(style="font-size: 12px; color: #7f8c8d;",
                     "Select proteins to filter the mutational count analysis.")
               ),
               mainPanel(
                 div(class = "justification-box",
                     h3("Bar Plot: Shiny Input Widgets"),
                     HTML("<p>We implemented <b>shiny input widgets</b> for the protein-level analysis. Since the virus consists of distinct 
                          functional regions, allowing the filtering of the dataset enables a direct comparison between 
                          different proteins. 
                          This dynamic filtering helps test hypotheses regarding conservation versus evolution in specific viral components.</p>")
                 ),
                 plotOutput("barPlot", click = "plot_click"),
                 hr(),
                 h4("Individual Observation Details (Click a bar to view):"),
                 tableOutput("click_info")
               )
             )
    ),
    
    tabPanel("Frequency Analysis",
             sidebarLayout(
               sidebarPanel(
                 h4("Data Filtering"),
                 sliderInput("freq_range", "Frequency Threshold:", 
                             min = 0, max = 1, value = c(0, 1), step = 0.01),
                 p("Adjust the slider to isolate rare vs. common variants.")
               ),
               mainPanel(
                 div(class = "justification-box",
                     h3("Histogram: Built-in Interactive Options"),
                     HTML("<p>The <b>sliderInput</b> provides a built-in interactive method to manipulate the data visualization in real-time. 
                          By filtering allele frequencies, we can effectively analyze the rare variant skew 
                          characteristic of rapid outbreaks. This interactivity allows for a more focused exploration of 
                          mutations that cross the 0.05 threshold of viability.</p>")
                 ),
                 plotOutput("histPlot")
               )
             )
    )
  )
)

server <- function(input, output, session) {
  
  # 1. Scatter Plot Logic
  output$scatterPlot <- renderPlotly({
    p <- ggplot(df_derived, aes(Position, Freq, color = Protein, 
                                text = paste("Position:", Position, "<br>Frequency:", Freq))) + 
      geom_point(alpha=0.7) +
      scale_y_log10() +
      scale_color_manual(values = my_colors) +
      theme_minimal() +
      labs(title="Zika Mutation Frequencies Across the Genome")
    ggplotly(p, tooltip = "text")
  })
  
  output$barPlot <- renderPlot({
    df_filt <- df_derived %>% filter(Protein %in% input$prot_select)
    ggplot(df_filt, aes(x = Protein, fill = Protein)) +
      geom_bar() +
      scale_fill_manual(values = my_colors) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Number of Mutations Identified per Protein", y = "Count")
  })
  
  selected_data <- reactiveVal(NULL)
  
  observeEvent(input$plot_click, {
    lvls <- levels(df_derived$Protein)
    clicked_index <- round(input$plot_click$x)
    
    if(clicked_index > 0 && clicked_index <= length(lvls)) {
      prot_name <- lvls[clicked_index]
      
      res <- df_derived %>%
        filter(Protein == prot_name) %>%
        select(Position, Protein, Freq, AA_Change) %>%
        head(15)
      
      selected_data(res)
    }
  })
  
  output$click_info <- renderTable({
    selected_data()
  })
  
  output$histPlot <- renderPlot({
    df_filt <- df_derived %>% filter(Freq >= input$freq_range[1] & Freq <= input$freq_range[2])
    
    ggplot(df_filt, aes(x = Freq)) +
      geom_histogram(bins = 30, aes(fill = after_stat(x))) +
      scale_x_log10() +
      scale_fill_gradientn(colors = c("red", "yellow" ,"green"),
                           values = scales::rescale(c(0, 0.05, 1))) +
      theme_classic() +
      labs(title = "Frequency Distribution of Variants", x = "Frequency (log10)", y = "Count")
  })
}

shinyApp(ui, server)
