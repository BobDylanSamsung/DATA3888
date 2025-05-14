server <- function(input, output, session) {
  pred_result <- reactiveValues()
  # Flag for whether we are showing the model or prediction
  pred_made <- reactiveVal(FALSE)
  output$show_model_plots <- reactive({ !pred_made() })
  outputOptions(output, "show_model_plots", suspendWhenHidden = FALSE)
  
  attach_render_outputs(
    output,
    cvfit      = cvfit,
    best_lambda= best_lambda,
    gene_coefs = gene_coefs,
    fit_km     = fit_km,
    train_df   = train_df,
    train_lp   = train_lp,
    pred_result= pred_result,
    gse_eda    = gse_eda,
    summary_per_sample_eda = summary_per_sample_eda,
    summary_long_eda       = summary_long_eda,
    pca_all                = pca_all,
    top100_genes_eda       = top100_genes_eda,
    annotation_col_eda     = annotation_col_eda
  )
  
  
  observeEvent(input$predict, {
    req(input$gene_csv)
    pred_made(TRUE)
    updateTabsetPanel(session, "main_tabs", selected = "Prediction")
    df <- read.csv(input$gene_csv$datapath, header = FALSE, stringsAsFactors = FALSE)
    
    # Remove CSV header row
    df <- df[-1, , drop = FALSE]
    genes <- df[[1]]
    values <- as.numeric(df[[2]])
    names(values) <- genes
    
    newdata <- as.data.frame(t(values[selected_genes]))
    colnames(newdata) <- selected_genes
    
    # Predict
    # PREDICTION
    pred_result$risk_score <- predict(coxph_model, newdata = newdata, type = "lp")
    pred_result$risk_group <- ifelse(pred_result$risk_score > cutoff, "High", "Low")
    pred_result$hr <- exp(pred_result$risk_score)
    pred_result$fit_new <- survfit(coxph_model, newdata = newdata)
    pred_result$newdata <- newdata 
  })
  
  # Reset to model plots if file re-uploaded
  observeEvent(input$gene_csv, pred_made(FALSE))
}