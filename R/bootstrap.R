#' S4 bootstrap class
#' @exportClass bootstrap
setClass("bootstrap", representation(bootstrapped_statistics = "data.frame"),
         validity = function(object) {
           errors <- c()
           if (any(is.na(colnames(object@bootstrapped_statistics))) || # sneakily relies on short circuit evaluation
               is.null(colnames(object@bootstrapped_statistics)) ||
               any(colnames(object@bootstrapped_statistics) == '')) {
              error <- "One or more bootstrapped statistic name(s) missing"
              errors <- c(errors, error)
           }
           
           if (length(errors) == 0) TRUE else errors
         })

#' S4 summary.bootstrap class
#' @exportClass summary.bootstrap
setClass("summary.bootstrap", representation(summary = "data.frame", 
                                             n_samples = "integer",
                                             interval_level = "numeric"),
         validity = function(object) {
           errors <- c()
           
           required_colnames <- c('statistic_name', 'mean', 'standard_error', 'ci_lower_bound', 'ci_upper_bound')
           summary_columns <- colnames(object@summary)
           if (length(intersect(required_colnames, summary_columns)) != length(summary_columns)) {
             error <- paste0("Bootstrap summary must include the following columns: ", paste(required_colnames, collapse = ','))
             errors <- c(errors, error)
           }
           
           if (object@interval_level >= .5 || object@interval_level < 1) {
             error <- paste0("Interval level must be inside range [.5, 1)")
             errors <- c(errors, error)
           }
           
           if (length(errors) == 0) TRUE else errors
         })

#' Perform boostrap resampling
#'
#' Perform bootstrap resampling for arbitrary statistics using tidy semantics
#' Note that the global mc.cores can be used to control level of parallelism in computing statistics on resamples.
#'
#' @param data data frame or matrix containing observations to be bootstrapped
#' @param ... statistics to compute 
#' @param trials number of resamplings to perform
#'
#' @author kholub
#' @examples
#' boot_results <- bootstrap(iris, sepal_ratio = sum(Sepal.Length) / sum(Sepal.Width),
#'                                 petal_ratio = sum(Petal.Length) / sum(Petal.Width))
#'
#' @importFrom dplyr sample_n
#' @importFrom parallel mclapply
#'
#' @export bootstrap
bootstrap <- function(data, ...,
                      trials = 1e2) {
  # note that the call object includes a sentinel first value - we ignore it
  statistics_call <- substitute(list(...))[-1]
  statistic_names <- names(statistics_call)
  
  if (length(statistic_names) < 1) {
    stop('Must specify one or more statistics to bootstrap')
  }
  else if (is.null(statistic_names) || any(statistic_names == '')) {
    stop('One or more statistic(s) specified without a name - all statistics must have a name')
  }
  
  results <- mclapply(1:trials, function(trial_ix) {
    resample <- sample_n(data, nrow(data), replace = TRUE)
    results_row <- sapply(unname(statistics_call), function(x){
      eval(x, resample, parent.frame())
    })
    
    results_row
  })
  
  # better alternatives, but keep it base R
  results <- data.frame(do.call(rbind, results))
  colnames(results) <- statistic_names
  
  new('bootstrap',
      bootstrapped_statistics = results)
}

#' Compute confidence intervals for a bootstrapped statistic
#'
#'
#' @export confint.bootstrap
confint.bootstrap <- function(object, parm, level = .95, ...) {
  summary_df <- as.data.frame(summary(object, empirical_interval_level = level))
  
  cis<- summary_df[summary_df$statistic_name %in% parm, c('ci_lower_bound', 'ci_lower_bound')]
  colnames(cis) <- as.character(c((1 - level) / 2, level + (1 - level) / 2))
  
  as.matrix(cis)
}

print.bootstrap <- function(x, ...) {
  
}

#' Summarize a bootstrap result
#'
#' @param object a bootstrap result to be summarized
#' @param empirical_ci_level the lower quantile used to compute empirical intervals on the bootstrapped statistic
#'
#' @author kholub
#' @examples 
#' boot_results <- bootstrap(iris, sepal_ratio = sum(Sepal.Length) / sum(Sepal.Width),
#'                                 petal_ratio = sum(Petal.Length) / sum(Petal.Width))
#' summary(boot_results)
#'
#' @importFrom dplyr bind_rows
#'
#' @export summary.bootstrap
summary.bootstrap <- function(object, level = .95, ...) {
  stopifnot(level >= .5 && level < 1)
  
  total_summary <- data.frame()
  for (statistic_name in colnames(object@bootstrapped_statistics)) {
    bootstrapped_statistic <- unlist(object@bootstrapped_statistics[[statistic_name]])
    statistic_mean = mean(bootstrapped_statistic)
    statistic_se = sd(bootstrapped_statistic)
    
    empirical_ci <- quantile(bootstrapped_statistic, c((1 - level) / 2, level + (1 - level) / 2))
    
    total_summary <- bind_rows(total_summary, list(statistic_name =  statistic_name,
                                                   mean = statistic_mean,
                                                   standard_error = statistic_se,
                                                   ci_lower_bound = empirical_ci[1],
                                                   ci_upper_bound = empirical_ci[2]
    ))
  }
  
  new('summary.bootstrap',
      summary = total_summary,
      n_samples = nrow(object@bootstrapped_statistics),
      interval_level = level)
}

setMethod("summary", signature("bootstrap"), summary.bootstrap)

#' Produce a dataframe of bootstrap summarization
#'
#' @author kholub
#'
#' @export as.data.frame.summary.bootstrap
as.data.frame.summary.bootstrap <- function(x, ...) {
  x@summary
}

#' Print a bootstrap result summary
#'
#' @param x an object of class summary.boostrap to be printed
#'
#' @author kholub
#'
#' @export print.summary.bootstrap
print.summary.bootstrap <- function(x, ...) {
  cat("Number of boostrap samples:", x@n_samples, '\n')
  cat("Empirical interval level:", x@interval_level, '\n\n')
  
  print(x@summary, row.names = FALSE)
}

setMethod("print", signature(x = "summary.bootstrap"), print.summary.bootstrap)
setMethod("show", signature(object = "summary.bootstrap"), function(object) { print.summary.bootstrap(object) })