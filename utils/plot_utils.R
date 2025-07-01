positivity_check_histograms <- function(df, treatment, confounders, bins = 30) {
  # ensure treatment is a factor for consistent coloring
  df[[treatment]] <- as.factor(df[[treatment]])
  
  for (conf in confounders) {
    # start the plot
    p <- ggplot(df, aes_string(x = conf, fill = treatment))
    
    # add the appropriate histogram/bar layer
    if (is.numeric(df[[conf]])) {
      p <- p +
        geom_histogram(
          position = "identity",
          alpha    = 0.5,
          bins     = bins,
          color    = "black",
          linewidth = 0.2
        )
    } else {
      p <- p +
        geom_bar(
          position = "identity",
          stat     = "count",
          alpha    = 0.5,
          color    = "black",
          linewidth = 0.2
        )
    }
    
    # labels & styling
    p <- p +
      labs(
        title = sprintf("Positivity Check for %s", conf),
        x     = conf,
        y     = "Count",
        fill  = "treatment"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "top",
        plot.title     = element_text(face = "bold", hjust = 0.5)
      )
    
    print(p)
  }
}

plot_ps_hist <- function(df,
                         ps_vec,
                         treatment_col,
                         bins = 30,
                         plot_title = "Propensity Score Distribution by Treatment") {
  stopifnot(length(ps_vec) == nrow(df))
  
  # 1) add PS column
  df$propensity_score <- ps_vec
  
  # 2) coerce treatment → factor, PS → numeric
  df[[treatment_col]]   <- as.factor(df[[treatment_col]])
  df$propensity_score   <- as.numeric(df$propensity_score)
  
  # 3) drop bad PS
  plot_df <- df[is.finite(df$propensity_score), , drop = FALSE]
  
  # 4) plot
  p <- ggplot(plot_df, aes(x = propensity_score, fill = !!rlang::sym(treatment_col))) +
    geom_histogram(
      position  = "identity",
      alpha     = 0.5,
      bins      = bins,
      color     = "black",
      linewidth = 0.2,
      na.rm     = TRUE
    ) +
    labs(
      title = plot_title,
      x     = "Propensity Score",
      y     = "Count",
      fill  = "treatment"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title     = element_text(face = "bold", hjust = 0.5)
    )
  
  print(p)
}

plot_iptw_weights <- function(df,
                              weights,
                              treatment_col,
                              bins       = 30,
                              plot_title = "IPTW Weights by Treatment") {
  stopifnot(length(weights) == nrow(df))
  
  # 1) add weight column
  df$iptw_weight <- as.numeric(weights)
  
  # 2) ensure treatment is a factor
  df[[treatment_col]] <- as.factor(df[[treatment_col]])
  
  # 3) drop non‐finite weights
  plot_df <- df[is.finite(df$iptw_weight), , drop = FALSE]
  
  # 4) plot
  p <- ggplot(plot_df, aes_string(x = "iptw_weight", fill = treatment_col)) +
    geom_histogram(
      position  = "identity",
      alpha     = 0.5,
      bins      = bins,
      color     = "black",
      linewidth = 0.2,
      na.rm     = TRUE
    ) +
    labs(
      title = plot_title,
      x     = "IPTW Weight",
      y     = "Count",
      fill  = "treatment"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title     = element_text(face = "bold", hjust = 0.5)
    )
  
  print(p)
}

plot_predicted_risk_hist <- function(df,
                                     pred_vec,
                                     group_col,          # <- NEW
                                     bins       = 30,
                                     plot_title = "Predicted outcome probability") {
  stopifnot(length(pred_vec) == nrow(df))
  
  plot_df <- df %>%
    mutate(
      pred_prob = as.numeric(pred_vec),
      !!group_col := factor(.data[[group_col]])
    )
  
  p <- ggplot(plot_df, aes(x = pred_prob, fill = !!rlang::sym(group_col))) +
    geom_histogram(
      position  = "identity",
      alpha     = 0.5,
      bins      = bins,
      colour    = "black",
      linewidth = 0.2
    ) +
    labs(
      title = plot_title,
      x     = "Predicted P(Y = 1)",
      y     = "Count",
      fill  = group_col                # legend title now “day_28_flg”
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title      = element_text(face = "bold", hjust = 0.5)
    )
  
  print(p)
  invisible(p)
}

