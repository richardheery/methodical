#' Created a Shiny app for displaying CpG methylation and transcript expression correlations
#' 
#' @param methrix_object A methrix object
#' @param transcript_table A data.frame with transcript expression values where row names are the names of the transcripts and column names are sample names. 
#' @param correlation_table A table with the calculated correlations for methylation of CpG sites and transcript expression, 
#' such as returned by combine_transcription_cpg_correlations()
#' @param sample_groups An optional factor giving the names of the groups that samples belong to. 
#' @return A Shiny App which creates scatter plots for input CpG names and transcripts or else from selected rows in the table.
#' @export
make_shiny_app = function(methrix_object, transcript_table, correlation_table, sample_groups = NULL){
  
  # Find samples in common between methrix object and transcript table
  common_samples = intersect(row.names(colData(methrix_object)), names(transcript_table))
  
  # Define UI for application that draws a histogram
  ui <- shiny::fluidPage(
  
      # Application title
      shiny::titlePanel("CpG Methylation-Transcript Expression Correlations"),
  
      # Sidebar with a slider input for number of bins 
      shiny::sidebarLayout(
          shiny::sidebarPanel(
            shiny::textInput("cpg_name", "Name of CpG", value = "chr1:24411422"),
            shiny::textInput("transcript_name", "Name of Transcript", value = "ENST00000003583"),
            shiny::actionButton("button", "Create Plot"),
            width = 3
          ),
  
          # Show a plot of the generated distribution
          shiny::mainPanel(
            plotly::plotlyOutput("corr_plot", width = "666px", height = "500px"),
            width = 8
          )
        ),
    DT::dataTableOutput("corr_table")
    )

  # Define server logic required to draw a histogram
  server <- function(input, output) {
  
        
        pp = shiny::eventReactive(input$button, {
          
        if(!is.null(input$corr_table_rows_selected)){
          cpg_name = correlation_table$cpg_name[input$corr_table_rows_selected]
          transcript_name = correlation_table$transcript_id[input$corr_table_rows_selected]
        } else {
          cpg_name = input$cpg_name
          transcript_name = input$transcript_name
        }
        
        # Get values for CpG site as a data.frame
        cpg_region = resize(GRanges(cpg_name), width = 2)
        cpg_values = methodical::extract_cpg_values_from_methrix(methrix = methrix_object, genomic_regions = cpg_region, samples = common_samples)
    
        # Get transcript values as a vector
        transcript_values = unlist(transcript_table[transcript_name, ])
      # title = paste(cpg_name, "Methylation and", transcript_name, "Expression")    
      return(list(
        cor_plot = plotly::layout(plotly::ggplotly(methodical::methylation_feature_scatter_plot(cpg_values = cpg_values, cpg_name = cpg_name, feature = transcript_values, 
          sample_groups = sample_groups, use_names = T, title = NULL, add_correlation = F, method = "p", alpha = 0.75)), legend = list(bgcolor = "#E2E2E2", x = 100, y = 0.5)))
        )
        })
          
     output$corr_plot = plotly::renderPlotly(pp()$cor_plot)
     
     output$corr_table = DT::renderDT(
       DT::datatable(correlation_table, selection = "single", filter = "top"))
  }
  
  # Make the application and return it
  shiny_app = shiny::shinyApp(ui = ui, server = server)
  return(shiny_app)
}


