source("src/server/attach_model_outputs.R")
source("src/server/prediction_helpers.R")
source("src/server/attach_eda_outputs.R")

# Server function for the Shiny app
server <- function(input, output, session) {
  # Initialize reactive values
  pred_made <- reactiveVal(FALSE)
  output$show_model_plots <- reactive({ !pred_made() })
  outputOptions(output, "show_model_plots", suspendWhenHidden = FALSE)
  
  # Attach visualizations to output
  attach_eda_outputs(output)
  attach_model_outputs(output)
  
  # Initialize prediction outputs
  clear_prediction_outputs(output)
  
  # Store model parameters for prediction
  model_params <- list(
    coxboost_model = coxboost_model,
    selected_features = selected_features,
    train_center = attr(X_train_selected_std, "scaled:center"),
    train_scale = attr(X_train_selected_std, "scaled:scale"),
    cutoff = median(risk_scores_test),
    risk_groups = risk_groups,
    risk_scores_test = risk_scores_test,
    lookup = lookup
  )
  
  # Handle prediction button click
  observeEvent(input$predict, {
    req(input$gene_csv)
    pred_made(TRUE)
    updateTabsetPanel(session, "main_tabs", selected = "Prediction")
    
    withProgress(message = 'Processing data...', {
      tryCatch({
        # Process the prediction with all parameters in a single object
        pred_results <- process_prediction(
          input$gene_csv$datapath, 
          model_params
        )
        
        # Attach the prediction outputs to the UI
        attach_prediction_outputs(
          output, 
          pred_results, 
          risk_groups,
          risk_scores_test, 
          c_index
        )
        
      }, error = function(e) {
        clear_prediction_outputs(output, e$message)
      })
    })
  })
  
  # Reset to model plots if file re-uploaded
  observeEvent(input$gene_csv, {
    pred_made(FALSE)
    clear_prediction_outputs(output)
    updateTabsetPanel(session, "main_tabs", selected = "Model")
  })
}