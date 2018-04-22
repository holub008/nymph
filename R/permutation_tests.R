#' S4 permutation test class
#' @exportClass prmt
setClass("prmt", representation(permutation_results = "data.frame",
                                original_permutation_results = "list"),
         validity = function(object) {
           errors <- c()
           if (any(is.na(colnames(object@permutation_results))) || # sneakily relies on short circuit evaluation
               is.null(colnames(object@permutation_results)) ||
               any(colnames(object@permutation_results) == '')) {
             error <- "One or more bootstrapped statistic name(s) missing"
             errors <- c(errors, error)
           }
           
           if (any(is.na(names(object@original_permutation_results))) ||
               is.null(names(object@original_permutation_results)) ||
               any(colnames(object@permutation_results) != names(object@original_permutation_results))) {
             error <- "One or more original permutation statistic name(s) missing"
             errors <- c(errors, error)
           }
           
           if (length(errors) == 0) TRUE else errors
         })

#' S4 permutation test summary  class
#' @exportClass summary.prmt
setClass("summary.prmt", representation(summary = "data.frame",
                                        permutations = "integer",
                                        level = "numeric",
                                        alternative = "character"),
         validity = function(object) {
           errors <- c()
           
           required_summary_columns <- c("statistic", "ci_lower", "ci_upper", "actual", "p_value")
           if (is.null(colnames(object@summary)) ||
               length(setdiff(required_summary_columns, colnames(object@summary))) != 0 &&
               length(required_summary_columns) != length(colnames(object@summary))) {
             errors <- c(errors,
                         paste0("Summary column names do not match expectation of ", paste(required_summary_columns, collapse = ',')))
           }
           
           if (any(object@summary$p_value < 0 | object@summary$p_value > 1)) {
             errors <- c(errors,
                         "Summary contains a p_value that does not represent a valid probability")
           }
           
           if (!(object@level <= 1 && object@level >= 0)) {
             errors <- c(errors,
                         "Significance level must fall in [0,1]")
           }
           
           valid_alternatives <- c('less', 'greater', 'two.sided')
           if (!(object@alternative %in% valid_alternatives)) {
             errors <- c(errors,
                         paste0("Alternative must be one of ", paste(valid_alternatives, collapse = ',')))
           }
           
           if (length(errors) == 0) TRUE else errors
         })

#' Visualize observed statistics against null distributions
#'
#' @param object an prmt object
#' @param nrows the number of rows in the resulting graphic (each entry corresponds to a single statistic)
#' @param ncols the number of columns in the resulting graphic
#'
#' @author kholub
#'
#' @examples
#' data(mpg)
#' sample1 <- mpg %>% filter(class == 'suv')
#' sample2 <- mpg %>% filter(class == 'compact')
#' mcr_results <- mcrd_test(sample1, sample2, median = median(hwy),
#'                                           mean = mean(hwy))
#' plot(mcr_results)
#'
#' @export plot.prmt
plot.prmt <- function(x, ..., nrows = 1, ncols = NULL) {
  ncols <- if (is.null(ncols)) ncol(x@permutation_results) else ncols
  
  previous_state <- par(mfrow = c(nrows, ncols))
  tryCatch({
    for (stat_name in colnames(x@permutation_results)) {
      stat_permuted <- unlist(x@permutation_results[[stat_name]])
      stat_observed <- unlist(x@original_permutation_results[[stat_name]])
      
      x_limits = c(min(stat_permuted, stat_observed),
                   max(stat_permuted, stat_observed))
      
      hist(stat_permuted, main = stat_name,
           xlab = "Statistic under permutation",
           xlim = x_limits)
      abline(v = stat_observed, col = 'blue')
    }
  },
  finally = { par(previous_state) })
}

setMethod("plot", signature("prmt"), plot.prmt)

#' Compute empirical confidence intervals for statistics under permutation testing
#'
#' @param object an object of class prmt
#' @param parm a vector of statistic names to produce intervals for
#' @param level the width of confidence intervals
#'
#' @author kholub
#' 
#' @examples
#' data(mpg)
#' sample1 <- mpg %>% filter(class == 'suv')
#' sample2 <- mpg %>% filter(class == 'compact')
#' mcr_results <- mcrd_test(sample1, sample2, median = median(hwy),
#'                                            mean = mean(hwy))
#' confint(mcr_results, parm = c("mean", "median"))
#'
#' @export confint.prmt
confint.prmt <- function(object, parm, level = .95, ...) {
  stopifnot(level >= 0 && level <= 1)
  
  summary_df <- as.data.frame(summary(object, level = level))
  
  cis<- summary_df[summary_df$statistic %in% parm, c('ci_lower', 'ci_upper')]
  percentiles <- as.character(100 * c((1 - level) / 2, level + (1 - level) / 2))
  colnames(cis) <- sapply(percentiles, function(x){ paste0(x,'%') })
  rownames(cis) <- parm
  
  as.matrix(cis)
}

setMethod("confint", signature("prmt"), confint.prmt)

#' Print a permutation test
#' 
#' @author kholub
#'
#' @export print.prmt
print.prmt <- function(x, ...) { 
  cat("Monte Carlo randomization test of", as.character(nrow(x@permutation_results)), 
      " trials on the following statistics:", paste0(colnames(x@permutation_results), collapse = ','))
}

setMethod("print", signature(x = "prmt"), print.prmt)
setMethod("show", signature(object = "prmt"), function(object) { print.prmt(object) })

#' Summarize the results of a permutation test
#'
#' Yields p values & CIs against null (empirical) distribution of statistics.
#' For two sided hypotheses, reported p-value is the minimum one sided p-value with a Bonferroni correction.
#'
#' @param object a prmt object to be summarized
#' @param level the width of empirical confidence intervals
#' @param alternative the alternative used for tests of significance. one of 'less', 'greater', 'two.sided'
#'
#' @author kholub
#'
#' @examples
#' data(mpg)
#' sample1 <- mpg %>% filter(class == 'suv')
#' sample2 <- mpg %>% filter(class == 'compact')
#' mcr_results <- mcrd_test(sample1, sample2, mean = mean(hwy))
#' summary(mcr_results)
#'
#' @importFrom dplyr bind_rows
#'
#' @export summary.prmt
summary.prmt <- function(object, level = .05, alternative = 'two.sided', ...) {
  stopifnot(level <= 1 && level >= 0)
  
  total_summary <- data.frame()
  for (statistic_name in colnames(object@permutation_results)) {
    permuted_differences <- unlist(object@permutation_results[[statistic_name]])
    observed_difference <- object@original_permutation_results[[statistic_name]]
    cdf <- ecdf(permuted_differences)
    
    lower_tail <- cdf(observed_difference)
    if (alternative == 'less') {
      empirical_p <- lower_tail
      empirical_ci <- c(-Inf, quantile(permuted_differences, level))
    }
    else if (alternative == 'greater') {
      empirical_p <- 1 - lower_tail
      empirical_ci <- c(quantile(permuted_differences, 1 - level), Inf)
    }
    else if (alternative == 'two.sided') {
      # by selecting the minimum, we are implicitly performing 2 tests, so we apply a Bonferroni correction
      empirical_p <- 2 * min(1 - lower_tail, lower_tail)
      empirical_ci <- quantile(permuted_differences, c(level / 2, 1 - level / 2))
    }
    else {
      stop('Unrecognized "alternative" argument. Must be one of "less", "greater", "two.sided"')
    }
    total_summary <- rbind(total_summary, list(statistic = statistic_name,
                                               ci_lower = empirical_ci[1],
                                               ci_upper = empirical_ci[2],
                                               actual = observed_difference,
                                               p_value = empirical_p),
                           stringsAsFactors = FALSE)
  }
  
  rownames(total_summary) <- NULL
  
  new('summary.prmt', 
      summary = total_summary,
      permutations = nrow(object@permutation_results),
      level = level,
      alternative = alternative)
}

setMethod("summary", signature("prmt"), summary.prmt)

#' Access permutation test summary data
#'
#' @author kholub
#'
#' @export as.data.frame.summary.prmt
as.data.frame.summary.prmt <- function(x, ...) { 
  x@summary
}

setMethod("as.data.frame", signature("summary.prmt"), as.data.frame.summary.prmt)

#' Print a permutation test summary
#' 
#' @author kholub
#'
#' @export print.summary.prmt
print.summary.prmt <- function(x, ...) { 
  cat("Permutation test of", x@permutations, "permutations against alternative of", x@alternative, "at significance", as.character(x@level), '\n')
  print(x@summary, row.names = FALSE)
}

setMethod("print", signature(x = "summary.prmt"), print.summary.prmt)
setMethod("show", signature(object = "summary.prmt"), function(object) { print.summary.prmt(object) })