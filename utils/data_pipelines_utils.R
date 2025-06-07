# Impute missing values in a data.frame

impute_df <- function(df, method = c("mean", "median")) {
  method <- match.arg(method)
  
  for (col in names(df)) {
    x <- df[[col]]
    
    # only act if there are NAs
    if (anyNA(x)) {
      
      # numeric columns → mean or median
      if (is.numeric(x)) {
        if (method == "mean") {
          imp <- mean(x, na.rm = TRUE)
        } else {
          imp <- median(x, na.rm = TRUE)
        }
        
        # character or factor → most frequent (mode)
      } else if (is.factor(x) || is.character(x)) {
        tbl <- table(x)
        imp <- names(tbl)[which.max(tbl)]
        
        # if original was factor, ensure imp is a factor level
        if (is.factor(x) && ! imp %in% levels(x)) {
          levels(x) <- c(levels(x), imp)
        }
        
        # other types (logical, dates, etc.) → skip
      } else {
        next
      }
      
      # replace NAs
      x[is.na(x)] <- imp
      df[[col]] <- x
    }
  }
  
  df
}