# Define wrapper function at the top of your file
withDefaultSpinner <- function(ui_element) {
  shinycssloaders::withSpinner(
    ui_element,
    type = 8,
    color = "#2E9FDF",
    size = 1
  )
}

model_tab <- tabItem(
  tabName = "model",
  # -----------------------------------------------------------------
  # 1. Summary value boxes (2x2 grid) and Log-Rank Test results side by side
  # -----------------------------------------------------------------
  fluidRow(
    # Left column - Value boxes in 2x2 grid
    column(width = 6,
           box(
             title = "Key Metrics",
             status = "primary",
             solidHeader = TRUE,
             width = NULL,
             div(
               fluidRow(
                 column(width = 6, withDefaultSpinner(valueBoxOutput("opt_steps_box", width = NULL))),   # optimal boosting steps
                 column(width = 6, withDefaultSpinner(valueBoxOutput("c_index_box", width = NULL)))      # test-set C-index
               ),
               fluidRow(
                 column(width = 6, withDefaultSpinner(valueBoxOutput("n_features_box", width = NULL))),  # number of selected genes
                 column(width = 6, withDefaultSpinner(valueBoxOutput("ibs_box", width = NULL)))          # integrated brier score
               )
             )
           )
    ),
    # Right column - Log-rank test results
    column(width = 6,
           box(
             title = "Log-Rank Test Results",
             status = "info",
             solidHeader = TRUE,
             width = NULL,
             withDefaultSpinner(tableOutput("log_rank_results")),
             p("The log-rank test compares survival distributions between risk groups. 
                A significant p-value (<0.05) indicates that the model effectively separates patients 
                into groups with different survival outcomes.")
           )
    )
  ),
  
  # -----------------------------------------------------------------
  # 2. Survival Analysis Visualizations
  # -----------------------------------------------------------------
  fluidRow(
    column(width = 6,
           box(
             title = "Kaplan-Meier Survival Curve",
             status = "primary",
             solidHeader = TRUE,
             width = NULL,
             withDefaultSpinner(plotOutput("km_curve", height = "400px")),
             p("This plot shows survival probability over time for each risk group. 
                Separation between curves indicates how well the model stratifies patients 
                into distinct prognostic groups. Wider separation suggests better risk discrimination.")
           )
    ),
    column(width = 6,
           box(
             title = "Cumulative Hazard Plot",
             status = "warning",
             solidHeader = TRUE,
             width = NULL,
             withDefaultSpinner(plotOutput("chaz_curve", height = "400px")),
             p("The cumulative hazard plot shows the accumulated risk of the event over time 
                for each risk group. Steeper slopes indicate higher rates of events, providing 
                insight into how risk accumulates differently between groups.")
           )
    )
  ),
  
  # -----------------------------------------------------------------
  # 3. Model performance
  # -----------------------------------------------------------------
  fluidRow(
    # ROC Curve - left side
    column(width = 6,
           box(
             title = "ROC Curve Analysis",
             status = "primary",
             solidHeader = TRUE,
             width = NULL,
             withDefaultSpinner(plotOutput("roc_curve_plot", height = "400px")),
             p("The ROC curve shows the tradeoff between sensitivity and specificity at 
                different classification thresholds. The Area Under the Curve (AUC) measures 
                overall discriminative ability, with values closer to 1 indicating better 
                performance.")
           )
    ),
    # Time-dependent AUC - right side
    column(width = 6,
           box(
             title = "Time Dependent AUC",
             status = "info",
             solidHeader = TRUE,
             width = NULL,
             withDefaultSpinner(plotOutput("time_auc_plot", height = "400px")),
             p("This plot shows how the model's discriminative ability (AUC) changes over time. 
                It helps identify when the model performs best and whether its predictive power 
                remains stable or deteriorates with longer follow-up.")
           )
    )
  ),
  fluidRow(
    box(
      title = "Brier Score Plot",
      status = "info",
      solidHeader = TRUE,
      width = 6,
      withDefaultSpinner(plotOutput("brier_curve")),
      p("The Brier score measures prediction accuracy at different time points, with lower 
         values indicating better calibration. This plot helps assess how well the predicted 
         probabilities match actual outcomes over the follow-up period. Values below 0.25 
         generally indicate good predictive accuracy.")
    ),
    box(
      title = "Feature Importance",
      status = "success",
      solidHeader = TRUE,
      width = 6,
      withDefaultSpinner(plotOutput("feature_importance_plot", height = "500px")),
      p("This plot shows the relative importance of each gene in the final model. 
       Features are sorted by their absolute effect size. Blue bars indicate genes associated with 
       improved survival (protective), while red bars show genes associated with worse survival (risk). 
       Longer bars represent stronger effects on patient prognosis.")
    )
  ),
  fluidRow(
    box(
      title = "Model Calibration",
      status = "primary",
      solidHeader = TRUE,
      width = 12,
      withDefaultSpinner(plotOutput("calibration_plot", height = "450px")),
      p("This calibration plot assesses how well the predicted survival probabilities match the observed outcomes. 
       Each point represents a group of patients with similar predicted risk. 
       Points close to the diagonal line indicate good calibration. 
       Above the line means the model underestimates risk, below means it overestimates risk.")
    )
  )
)