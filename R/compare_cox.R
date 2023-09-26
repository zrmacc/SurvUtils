# Purpose: Compare Cox models via cross-validation.
# Updated: 2023-09-26


#' Compare Cox Models
#' 
#' Compares Cox models based on two sets of predictors with respect to 
#' their C-statistics via K-fold cross validation.
#' 
#' @param status Event status.
#' @param time Event time.
#' @param x1 Model matrix for the first model.
#' @param x2 Model matrix for the second model.
#' @param folds Cross validation folds.
#' @param simple_cstat If TRUE uses the unweighted C-stat, otherwise
#'   uses the inverse probability of censoring weighted C-stat.
#' @export 
CompreCoxCstat <- function(
  status,
  time,
  x1,
  x2,
  folds = 10,
  simple_cstat = TRUE
) {
  
  # Data frames.
  df1 <- data.frame(time, status, x1)
  df2 <- data.frame(time, status, x2)
  
  # Cross validation folds.
  fold_id <- seq_len(length(time)) %% folds
  df1_split <- split(df1, fold_id) 
  df2_split <- split(df2, fold_id)
  
  fold_results <- lapply(seq_len(folds), function(i) {

    # Evaluation data.
    eval1 <- df1_split[[i]]
    eval2 <- df2_split[[i]]
    
    # Training data.
    train1 <- do.call(rbind, df1_split[-i])
    train2 <- do.call(rbind, df2_split[-i])
    
    # Model fits.
    fit1 <- survival::coxph(
      survival::Surv(time, status) ~ ., data = train1)
    fit2 <- survival::coxph(
      survival::Surv(time, status) ~ ., data = train2)
    
    # Predict.
    eval1$risk <- stats::predict(object = fit1, newdata = eval1)
    eval2$risk <- stats::predict(object = fit2, newdata = eval2)
    
    # Cstat.
    if (simple_cstat) {
      cstat1 <- Cstat(data = eval1)
      cstat2 <- Cstat(data = eval2)
    } else {
      cstat1 <- WeightedCstat(train_data = train1, test_data = eval1)
      cstat2 <- WeightedCstat(train_data = train2, test_data = eval2)
    }
    
    # Output.
    out <- data.frame(fold = i, cstat1 = cstat1, cstat2 = cstat2)
    return(out)
  })
  
  # Overall results.
  out <- do.call(rbind, fold_results)
  out$diff <- out$cstat1 - out$cstat2
  out$ratio <- out$cstat1 / out$cstat2
  return(out)
}
