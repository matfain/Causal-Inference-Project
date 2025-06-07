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
        fill  = treatment
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "top",
        plot.title     = element_text(face = "bold", hjust = 0.5)
      )
    
    print(p)
  }
}