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
             p("The plot shows the Kaplan–Meier survival curves for high and low risk groups predicted by the model (blue = low risk, yellow = high risk). 
             The x-axis is follow-up time (months), and the y-axis is survival probability. The shaded areas around the curves represent the 95% confidence intervals. 
             The survival curve for the high risk group drops faster and ends lower, suggesting a worse prognosis. The log-rank test gives a p-value of 0.097, 
             which does not reach the 0.05 level of statistical significance. This means that, under the current sample size and model, 
             there is a trend of difference between the two groups, but it is not strong enough. Increasing the sample size may help better distinguish between groups.
")
           )
    ),
    column(width = 6,
           box(
             title = "Cumulative Hazard Plot",
             status = "warning",
             solidHeader = TRUE,
             width = NULL,
             withDefaultSpinner(plotOutput("chaz_curve", height = "400px")),
             p("This cumulative hazard plot shows how events accumulate over time for the high and low risk groups. 
             The x-axis represents follow-up time (months), and the y-axis shows cumulative hazard. The high-risk group’s curve rises steeply early on, 
             with a much higher cumulative hazard than the low-risk group. In contrast, the low-risk group increases more slowly, 
             and the gap between the two groups grows over time. This plot helps visualize how events accumulate across different risk levels, 
             and shows that the model can effectively separate high- and low-risk patients: those in the high-risk group experience events earlier and more frequently, 
             demonstrating the model’s good risk stratification ability.
")
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
             p("This plot shows the ROC curve for the model’s prediction of 12-month mortality risk. The x-axis represents the false positive rate 
               (1 – specificity), and the y-axis shows the true positive rate (sensitivity). The diagonal dashed line represents the baseline of random prediction. 
               The area under the curve (AUC) is 0.769 (95% CI: 0.503–1), indicating that the model has a moderate ability to distinguish between outcomes. ")
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
             p("This plot shows the time-dependent AUC curve and its confidence intervals for the model’s survival prediction on the test set at different 
             follow-up time points (6, 12, 18, and 24 months). The AUC is approximately 0.62 at 6 months, increases to around 0.77 at 12 months, 
             peaks at about 0.82 at 18 months, and slightly drops to 0.79 at 24 months. As follow-up time increases, the model’s discriminative 
             ability improves gradually, reaches its best performance in the middle period, and then slightly declines. The wide confidence intervals, 
             especially at early time points, suggest estimation uncertainty due to limited sample size. This plot helps identify when the model performs 
             best and assess its stability over time.
")
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
      p("This plot shows the time-dependent Brier Score of the model on the test set at different follow-up time points. 
      The x-axis represents follow-up months, and the y-axis is the Brier Score. The dashed line at 0.25 serves as a reference threshold. 
      The Brier Score starts around 0.13, increases gradually to about 0.22 at 30 months, then steadily decreases to around 0.10 by 48 months, 
      and remains below 0.25 throughout. This indicates that the predicted survival probabilities match the actual outcomes well, 
      with better calibration and more stable prediction accuracy in the later follow-up periods.
")
    ),
    box(
      title = "Feature Importance",
      status = "success",
      solidHeader = TRUE,
      width = 6,
      withDefaultSpinner(plotOutput("feature_importance_plot", height = "500px")),
      p("This plot ranks gene importance based on the absolute values of their coefficients in the CoxBoost model. The x-axis shows the regression coefficients, 
        which reflect each gene's impact on survival outcomes. Positive values (yellow) indicate that higher expression of the gene is associated with worse prognosis, 
        while negative values (blue) suggest a protective effect. The most influential gene is AGRN with a coefficient of about +0.23, 
        indicating that its high expression significantly increases event risk. On the protective side, PRDM16 and 
        CAMTA1 show strong negative effects, suggesting their high expression is linked to better survival. 
        NMNAT1 also contributes protective value. The impact of other genes is weaker. 
        Overall, these four genes are the main drivers of the model’s risk stratification.")
    )
  ),
  fluidRow(
    box(
      title = "Model Calibration",
      status = "primary",
      solidHeader = TRUE,
      width = 12,
      withDefaultSpinner(plotOutput("calibration_plot", height = "450px")),
      p("This calibration plot divides the predicted 6-month survival probabilities in the test set into two groups, sorted in ascending order, 
      with about 10 patients per group (point size reflects sample size). The x-axis shows the model’s average predicted survival for each group, 
      and the y-axis shows the actual observed survival. Error bars represent 95% confidence intervals. 
      Perfect calibration would align with the diagonal dashed line. 
      The predicted points for both groups are close to the diagonal line in the plot, indicating that the model's 6-month survival predictions for patients are relatively accurate.
      The low-risk group lies slightly below the line (predicted ≈ 0.80, observed ≈ 0.75), suggesting a slight overestimation of about 5%. 
      Overall, the error bars are short and all intersect the diagonal, indicating that the model has very good calibration for predicting 6-month survival 
      across different risk levels.
")
    )
  )
)