library(tidyverse)
library(pwr)


#' Calculate Proportion of Correct Classification (PCC)
#'
#' Computes the Proportion of Correct Classification (PCC) based on class label counts.
#'
#' The PCC is calculated as the sum of squared class proportions. It represents the probability 
#' that a randomly drawn instance from the dataset belongs to the most frequent class.
#'
#' @param class_labels_tbl A table containing the counts of class labels.
#'
#' @return A numeric value representing the Proportion of Correct Classification (PCC).
#' @examples
#' class_labels <- table(c("A", "A", "B", "A", "B", "B", "B"))
#' calculate_pcc(class_labels)
#' # Output: 0.4285714
#'
#' @export
calculate_pcc <- function(class_labels_tbl) {
  class_proportions <- prop.table(class_labels_tbl)
  
  pcc <- sum(class_proportions^2)
  
  return(pcc)
}


#' Calculate Cohen's h Effect Size
#'
#' Computes Cohen's h, which is a standardized effect size used to measure the difference 
#' between two proportions.
#'
#' Cohen's h is calculated as:
#' \deqn{h = 2 \times (\asin(\sqrt{p1}) - \asin(\sqrt{p2}))}
#' where `p1` and `p2` are proportions.
#'
#' This measure is useful in power analysis and hypothesis testing involving proportions.
#'
#' @param p1 Numeric. A proportion (must be between 0 and 1).
#' @param p2 Numeric. Another proportion (must be between 0 and 1).
#'
#' @return A numeric value representing Cohen's h.
#' @examples
#' calculate_cohens_h(0.6, 0.4)
#' # Output: 0.413
#'
#' calculate_cohens_h(0.8, 0.5)
#' # Output: 0.684
#'
#' @export
calculate_cohens_h <- function(p1, p2) {
  if (p1 < 0 || p1 > 1 || p2 < 0 || p2 > 1) {
    stop("Proportions p1 and p2 must be between 0 and 1.")
  }
  
  h <- 2 * (asin(sqrt(p1)) - asin(sqrt(p2)))
  
  return(h)
}


#' Calculate Sample Size Across Accuracy Levels
#'
#' Computes the required sample size (`n`) for a given set of desired accuracies based on category counts.  
#' It considers class distribution, computes a multiplier for minimum target accuracy,  
#' and calculates the required sample size using Cohen's h and power analysis.
#'
#' @param name Character. A name for the variable (default: "Variable").
#' @param cat_counts Numeric vector. The counts of each category in the dataset (length 2 or 3).
#' @param accs Numeric vector. Desired accuracy levels (e.g., 0.75, 0.85).
#' @param alpha Numeric. Significance level for power analysis (default: 0.05).
#' @param power Numeric. Statistical power for sample size calculation (default: 0.80).
#'
#' @return A data frame with:
#'   \item{name}{Variable name}
#'   \item{cat_1, cat_2}{Counts for first two categories}
#'   \item{distribution}{Proportion of the majority class}
#'   \item{multiplier}{Adjustment factor based on class imbalance}
#'   \item{pcc}{Proportion of correct classifications (baseline)}
#'   \item{min_target_acc}{Minimum achievable accuracy based on class balance}
#'   \item{desired_accs}{List of desired accuracies, including minimum target accuracy}
#'   \item{n}{Required sample size for each accuracy level}
#'
#' If class imbalance is too extreme, the function returns an error message.
#'
#' @examples
#' calculate_n_across_accs(cat_counts = c(50, 100), accs = c(0.75, 0.85))
#'
#' @importFrom pwr pwr.p.test
#' @export
calculate_n_across_accs <- function(name='Variable', cat_counts, accs, alpha=0.05, power=0.80) {
    
    if (length(cat_counts) == 2) {
        cat_1 <- cat_counts[1]
        cat_2 <- cat_counts[2]
        distribution <- pmax(cat_1, cat_2) / (cat_1 + cat_2)
        
    } else if (length(cat_counts) == 3) {
        cat_1 <- cat_counts[1]
        cat_2 <- cat_counts[2]
        cat_3 <- cat_counts[3]
        distribution <- pmax(cat_1, cat_2, cat_3) / (cat_1 + cat_2 + cat_3)
    }
    
    multiplier <- ifelse(distribution >= 0.90, 1.15, 
                         ifelse(distribution >= 0.80, 1.20, 
                                ifelse(distribution >= 0.50, 1.25, NA)))
    
    pcc <- calculate_pcc(cat_counts)
    min_target_acc <- pcc * multiplier
    
    if (min_target_acc >= 1) {
        message(paste0('Outcome is too imbalanced; minimum target accuracy is too high: ', round(min_target_acc * 100, 2), '%'))
    } else {
        desired_accs <- c(min_target_acc, accs)
        
        res <- tryCatch({
            hs <- sapply(desired_accs, function(acc) calculate_cohens_h(acc, pcc))
            ns <- sapply(hs, function(h) pwr.p.test(h = h, sig.level = alpha, power = power, alternative = "greater")$n * 10)
            ns <- ceiling(ns)  # Round up
            
            data.frame(
              name = name,
              cat_1 = cat_1,
              cat_2 = cat_2,
              distribution = distribution,
              multiplier = multiplier,
              pcc = pcc,
              min_target_acc = min_target_acc,
              desired_accs = desired_accs,
              n = ns
            )
        }, error = function(e) {
            
            if (min_target_accuracy > 0.90) {
                message("Consider lowering your range of desired accuracies. The minimum target acc is too high: ", min_target_acc)
                data.frame(
                  name = name,
                  cat_1 = cat_1,
                  cat_2 = cat_2,
                  distribution = distribution,
                  multiplier = multiplier,
                  pcc = pcc,
                  min_target_acc = min_target_acc,
                  desired_accs = desired_accs
                )
            }  else {
                message(e)
            }   
        })
        
        return(res)
    }
}
                         

#' Check Post-Hoc Statistical Power of a Model Compared to Baseline
#'
#' This function calculates the post-hoc power of a model's observed accuracy compared to a 
#' predefined baseline Proportion of Correct Classifications (PCC). If the power is at least 80%, 
#' it prints a message indicating that the model performs significantly better than the baseline at 
#' the specified family-wise error rate (FWER).
#'
#' @param pcc Numeric. The baseline proportion of correct classifications (e.g., from a majority class baseline).
#' @param observed_acc Numeric. The observed accuracy of the model being evaluated.
#'
#' @return Prints a message if post-hoc power is at least 80\%.
#' 
#' @details 
#' This function uses Cohen's h to measure the effect size between the observed accuracy 
#' and the baseline PCC. It then uses a power calculation (`pwr.p.test`) assuming a one-sided test.
#' If the post-hoc power is greater than or equal to 0.80, it prints a message indicating
#' the model's performance improvement percentage over the baseline.
#'
#' @note This function assumes the presence of external variables:
#' \itemize{
#'   \item \code{df_res}: A data frame from which the sample size (number of rows) is derived.
#'   \item \code{fwer_alpha}: The significance level used in power analysis.
#'   \item \code{alpha}: Displayed in the printed message to indicate the significance threshold.
#' }
#'
#' @seealso \code{\link[pwr]{pwr.p.test}}, \code{\link{calculate_cohens_h}}
#'
#' @examples
#' \dontrun{
#' df_res <- data.frame(pred = c(1, 0, 1), true = c(1, 0, 0))
#' fwer_alpha <- 0.05
#' alpha <- 0.05
#' check_posthoc_power(pcc = 0.65, observed_acc = 0.78)
#' }                        
check_posthoc_power <- function(pcc, observed_acc) {
    # Check power if at least 80%
    h <- calculate_cohens_h(observed_acc, pcc)
    power <- pwr.p.test(h = h, sig.level = fwer_alpha, alternative = "greater", n = nrow(df_res))$power
    
    if (power>=0.80) {
        multiplier  <- round(observed_acc/pcc*100, 2)
        print(paste0('Model performs ', multiplier, '% better than baseline PCC at controlled FWER alpha = ', alpha, ' and at least the standard 80% power.')) 
    }  
}