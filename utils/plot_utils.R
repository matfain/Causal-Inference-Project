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


# ---------------------------------------------------------------------------
# Plot helpers For CATE Notebooks 
# ---------------------------------------------------------------------------
plot_cate_vs_modifier <- function(df, cate_vec, modifier, learner_name, modifier_labels, ylim_range = c(-0.2, 0.2)) {
  
  plot_df <- df %>% mutate(tau = cate_vec, mod = .data[[modifier]])
  
  if (modifier == "sofa_first_cat") {
    plot_df$mod <- factor(plot_df$mod, levels = c("0-2.5", "2.5-5", "5-7.5", "7.5-10", "10+"))
  }
  
  x_label <- modifier_labels[[modifier]] %||% modifier
  
  if (is.numeric(plot_df$mod)) {
    p <- ggplot(plot_df, aes(mod, tau)) +
      geom_point(alpha = .2) +
      geom_smooth(method = "loess", se = FALSE, colour = "red") +
      ylim(ylim_range[1], ylim_range[2]) +
      labs(x = x_label, y = "Estimated CATE")
  } else {
    p <- ggplot(plot_df, aes(mod, tau, fill = mod)) +
      geom_boxplot(outlier.alpha = .15) +
      ylim(ylim_range[1], ylim_range[2]) +
      guides(fill = "none") +
      labs(x = x_label, y = "Estimated CATE")
  }
  print(p); invisible(p)
}

plot_cate_grid <- function(grid_df, modifier, learner_name) {
  
  if (is.numeric(grid_df$mod_value)) {
    p <- ggplot(grid_df, aes(mod_value, cate_mean)) +
      geom_line() + geom_point() +
      labs(title = sprintf("Mean CATE vs %s (%s)", modifier, learner_name),
           x = modifier, y = "Mean CATE")
  } else {
    p <- ggplot(grid_df, aes(mod_value, cate_mean)) +
      geom_col(fill = "steelblue") +
      labs(title = sprintf("Mean CATE by %s (%s)", modifier, learner_name),
           x = modifier, y = "Mean CATE")
  }
  print(p); invisible(p)
}



# [Previous functions remain the same...]

# ---------------------------------------------------------------------------
# NEW: Additional CATE Analysis Functions
# ---------------------------------------------------------------------------

# 1. Compare CATE distributions across all learners
compare_cate_distributions <- function(s_results, t_results) {
  # Combine all CATE estimates
  all_cates <- bind_rows(
    tibble(tau = s_results[["S-rf"]]$tau, learner = "S-learner (RF)"),
    tibble(tau = s_results[["S-xgb"]]$tau, learner = "S-learner (XGB)"),
    tibble(tau = t_results[["rf"]]$tau, learner = "T-learner (RF)"),
    tibble(tau = t_results[["xgb"]]$tau, learner = "T-learner (XGB)")
  )
  
  # Remove non-finite values
  all_cates <- all_cates %>% filter(is.finite(tau))
  
  # Calculate reasonable x-axis limits (trim extreme outliers)
  tau_range <- quantile(all_cates$tau, c(0.01, 0.99), na.rm = TRUE)
  x_limits <- c(tau_range[1] - 0.05, tau_range[2] + 0.05)
  
  # Density plot
  p <- ggplot(all_cates, aes(x = tau, fill = learner)) +
    geom_density(alpha = 0.4, adjust = 1.5) +  # Smoother density
    geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
    labs(
         x = "Estimated CATE", 
         y = "Density") +
    theme_minimal(base_size = 12) +
    xlim(x_limits) +
    theme(legend.position = "bottom")
  
  print(p)
  invisible(p)
}


# 2. Plot CATE by pairs of modifiers
plot_cate_by_modifier_pairs <- function(df, tau_vec, mod1, mod2, learner_name, modifier_labels = list()) {
  df_plot <- df %>%
    mutate(tau = tau_vec)
  
  # Get labels
  lab1 <- modifier_labels[[mod1]] %||% mod1
  lab2 <- modifier_labels[[mod2]] %||% mod2
  
  # Check variable types
  is_num1 <- is.numeric(df_plot[[mod1]])
  is_num2 <- is.numeric(df_plot[[mod2]])
  
  # For continuous × categorical
  if (is_num1 && !is_num2) {
    # Recode factor levels for legend
    if (mod2 %in% names(modifier_labels)) {
      df_plot[[mod2]] <- factor(df_plot[[mod2]], 
                                levels = unique(df_plot[[mod2]]),
                                labels = paste(lab2, "=", unique(df_plot[[mod2]])))
    }
    
    p <- ggplot(df_plot, aes_string(x = mod1, y = "tau", color = mod2)) +
      geom_point(alpha = 0.3) +
      geom_smooth(method = "loess", se = FALSE) +
      labs(title = paste("CATE by", lab1, "and", lab2, "-", learner_name),
           x = lab1,
           y = "Estimated CATE",
           color = lab2) +
      theme_minimal()
  }
  # For categorical × continuous  
  else if (!is_num1 && is_num2) {
    # Recode factor levels for legend
    if (mod1 %in% names(modifier_labels)) {
      df_plot[[mod1]] <- factor(df_plot[[mod1]], 
                                levels = unique(df_plot[[mod1]]),
                                labels = paste(lab1, "=", unique(df_plot[[mod1]])))
    }
    
    p <- ggplot(df_plot, aes_string(x = mod2, y = "tau", color = mod1)) +
      geom_point(alpha = 0.3) +
      geom_smooth(method = "loess", se = FALSE) +
      labs(title = paste("CATE by", lab2, "and", lab1, "-", learner_name),
           x = lab2,
           y = "Estimated CATE",
           color = lab1) +
      theme_minimal()
  }
  # For categorical × categorical
  else if (!is_num1 && !is_num2) {
    # Create labels for x-axis if service_unit or sofa_first_cat
    if (mod1 == "service_unit") {
      df_plot[[mod1]] <- factor(df_plot[[mod1]], 
                                levels = c("MICU", "SICU"),
                                labels = c("Medical ICU", "Surgical ICU"))
    }
    if (mod2 %in% names(modifier_labels) && mod2 != "service_unit" && mod2 != "sofa_first_cat") {
      df_plot[[mod2]] <- factor(df_plot[[mod2]], 
                                levels = unique(df_plot[[mod2]]),
                                labels = paste(lab2, "=", unique(df_plot[[mod2]])))
    }
    
    # Create summary with sample sizes
    summary_df <- df_plot %>%
      group_by(!!sym(mod1), !!sym(mod2)) %>%
      summarize(mean_tau = mean(tau), 
                se_tau = sd(tau)/sqrt(n()), 
                n = n(),
                .groups = 'drop')
    
    p <- ggplot(summary_df, aes_string(x = mod1, y = "mean_tau", fill = mod2)) +
      geom_col(position = "dodge") +
      geom_errorbar(aes(ymin = mean_tau - se_tau, ymax = mean_tau + se_tau),
                    position = position_dodge(0.9), width = 0.2) +
      geom_text(aes(label = paste0("n=", n)), 
                position = position_dodge(0.9), 
                vjust = -0.5,
                size = 3) +
      labs(title = paste("Mean CATE by", lab1, "and", lab2, "-", learner_name),
           x = lab1,
           y = "Mean CATE",
           fill = lab2) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))
  }
  # For continuous × continuous
  else {
    # Create 2D heatmap
    df_binned <- df_plot %>%
      mutate(
        mod1_bin = cut(!!sym(mod1), breaks = 10),
        mod2_bin = cut(!!sym(mod2), breaks = 10)
      ) %>%
      group_by(mod1_bin, mod2_bin) %>%
      summarize(mean_cate = mean(tau, na.rm = TRUE), .groups = 'drop')
    
    p <- ggplot(df_binned, aes(x = mod1_bin, y = mod2_bin, fill = mean_cate)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      labs(title = paste("CATE Heatmap:", lab1, "vs", lab2, "-", learner_name),
           x = lab1,
           y = lab2,
           fill = "Mean CATE") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  print(p)
  invisible(p)
}

# 3. Analyze extreme CATE patients
analyze_cate_extremes <- function(df, tau_vec, learner_name, n_extreme = 50, modifier_labels = list()) {
  df_cate <- df %>%
    mutate(
      tau = tau_vec,
      patient_id = row_number(),
      tau_percentile = percent_rank(tau)
    )
  
  # Get top and bottom responders
  top_responders <- df_cate %>%
    arrange(desc(tau)) %>%
    slice_head(n = n_extreme) %>%
    mutate(group = "Top beneficiaries")
  
  bottom_responders <- df_cate %>%
    arrange(tau) %>%
    slice_head(n = n_extreme) %>%
    mutate(group = "Potential harm")
  
  extreme_patients <- bind_rows(top_responders, bottom_responders)
  
  # Compare key characteristics
  vars_to_compare <- c("age", "sapsi_first", "sofa_first_cat", "service_unit")
  
  plots <- list()
  for (var in vars_to_compare) {
    # Get label
    var_label <- modifier_labels[[var]] %||% var
    
    if (is.numeric(extreme_patients[[var]])) {
      p <- ggplot(extreme_patients, aes_string(x = "group", y = var, fill = "group")) +
        geom_boxplot() +
        labs(title = paste(var_label, "by CATE group -", learner_name),
             y = var_label) +
        theme_minimal()
    } else {
      p <- extreme_patients %>%
        count(group, !!sym(var)) %>%
        ggplot(aes_string(x = var, y = "n", fill = "group")) +
        geom_col(position = "dodge") +
        labs(title = paste(var_label, "distribution by CATE group -", learner_name),
             x = var_label,
             y = "Count") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    plots[[var]] <- p
    print(p)
  }
  
  # Summary statistics table
  summary_stats <- extreme_patients %>%
    group_by(group) %>%
    summarize(
      mean_tau = mean(tau),
      mean_age = mean(age),
      mean_sapsi = mean(sapsi_first),
      prop_micu = mean(service_unit == "MICU"),
      prop_chf = mean(chf_flg == 1),
      .groups = 'drop'
    )
  
  print(kable(summary_stats, caption = paste("Extreme CATE Groups -", learner_name)))
  
  invisible(list(plots = plots, extreme_patients = extreme_patients, summary = summary_stats))
}
