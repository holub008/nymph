#' S4 permutation test power class
#' @exportClass prm_power
setClass("prm_power", representation(experiments = "list",
                                     sample_sizes = "data.frame"), # allow for multiple groups in the case of imbalanced designs
         validity = function(object) {
           errors <- c()
           
           if (length(object@experiments) < 1) {
             error <- "Must simulate at least one experiment to determine power"
           }
           
           if ( is(object@experiments[[1]], 'prmt')) {
           }
           
           if (!statistic_names_valid(distinct(object@experiments[[1]]$statistic))) {
             error <- "One or more statistic name(s) missing from permutation test results"
             errors <- c(errors, error)
           }
           
           if (nrow(sample_sizes) < 1) {
             error <- "Must provide valid sample sizes"
             errors <- c(errors, error)
           }
           
           if (length(errors) == 0) TRUE else errors
         })

#' S4 randomization test power class
#' @exportClass mcr_power
setClass("mcr_power", contains = "prm_power")

#' S4 permutation test power summary class
#' @exportClass summary.prm_power
setClass("summary.prm_power", representation(statistic_summary = "data.frame",
                                             sample_sizes = "data.frame", # allow for multiple groups in the case of imbalanced designs
                                             alternative = 'character',
                                             alpha = 'numeric',
                                             experiments = 'integer'), 
         validity = function(object) {
           errors <- c()
           
           required_summary_columns <- c('statistic', 'power', 'average_effect', 'alternative', 'alpha')
           if ( length(setdiff(required_simulation_columns, colnames(object@statistic_summary))) != 0 &&
                length(required_simulation_columns) != length(colnames(object@statistic_summary))) {
             error <- 'Missing required statistic summary columns'
             errors <- c(errors, error)
           }
           
           if (!statistic_names_valid(object@statistic_summary$statistic)) {
             error <- "One or more statistic name(s) missing from power summary"
             errors <- c(errors, error)
           }
           
           if (nrow(sample_sizes) < 1) {
             error <- "Must provide valid sample sizes"
             errors <- c(errors, error)
           }
           
           if (length(errors) == 0) TRUE else errors
         })

#' Perform a power analysis of a randomization test
#'
#' @author kholub
#'
#' @export mcr_power
mcr_power <- function(generate_data, statistics, group_var, alternative,
                      statistic_names = names(statistics),
                      test_trials = 1e2,
                      power_trials = 1e2) {
  stopifnot(power_trials > 0)
  stopifnot(test_trials > 0)
  
  # since mcr_test can run in parallel, always run serially to avoid nested parallelization
  experiments <- lapply(1:power_trials, function(trial_ix) {
    mcr_test(generate_data(), statistics, group_var, statistic_names, test_trials)
  })
  
  group_sizes <- as.data.frame(table(generate_data()[,group_var]))
  name(group_sizes) <- c(group_var, 'size')
   
  new('mcr_power',
      simulation_results = simulations_df,
      sample_sizes = group_sizes,
      alternative = alternative)
}

#' Perform a power analysis of a randomization test
#' Note this is a convinience wrapper for \code{\link{mcr_power}} in the case of two sample differences
#'
#' @author kholub
#'
#' @export mcrd_power
mcrd_power <- function(generate_sample1, generate_sample2, alternative, ...,
                       statistic_names = names(statistics),
                       test_trials = 1e2,
                       power_trials = 1e2) {
  stopifnot(power_trials > 0)
  stopifnot(test_trials > 0)
  
  # since mcr_test can run in parallel, always run serially here to avoid nested parallelization
  experiments <- lapply(1:power_trials, function(trial_ix) {
    mcrd_test(generate_sample1(), generate_sample2, ..., group_var, statistic_names, test_trials)
  })
  
  group_sizes <- data.frame(sample = c('1', '2'), 
                            size = c(nrow(generate_sample1()), nrow(generate_sample2())))
  
  new('mcr_power',
      simulation_results = simulations_df,
      sample_sizes = group_sizes,
      alternative = alternative)
}

#' Visualize the power analysis of a randomization test
#' If alpha is left unspecified, the distribution of p-values is visualized. If alpha is specified, power is plotted as a function of alpha
#'
#' @author kholub
#'
#' @export plot.prmt
plot.prm_power <- function(x, ..., statistic, alternative, alpha = seq(.01, .99, by = .01)) {
  stopifnot(statistic %in% distinct(object@experiments[[1]]$statistic))
  stopifnot(is.null(alpha) || is.na(alpha) || length(alpha) > 0)
  
  if (is.null(alpha) || is.na(alpha)) {
    p <- sapply(x@experiments, function(experiment){
      s <- as.data.frame(summary(experiment, alternative = alternative))
      s[s$statistic == statistic, 'p_value']
    })
    
    hist(p, 
         main = paste0('p-values for', statistic),
         xlab = 'p-value',
         freq = FALSE)
  }
  else {
    power <- sapply(alpha, function(sig) {
      s <- as.data.frame(summary(x, alternative = alternative, alpha = sig))
      s[s$statistic == statistic, 'power']
    })
    
    plot(alpha, power, 
         main = paste0('Power analysis of ', statistic),
         xlab = 'alpha',
         ylab = 'power')
  }
}

setMethod("plot", signature("prm_power"), plot.prmt)

#' Print a permutation test power analysis
#' 
#' @author kholub
#'
#' @export print.prm_power
print.prm_power <- function(x, ...) { 
  cat("Power analysis of a permutation test of", as.character(length(x@experiments)), 
      " experiments on the following statistics:", paste0(distinct(x@experiments[[1]]$statistic), collapse = ','))
}

setMethod("print", signature(x = "prm_power"), print.prm_power)
setMethod("show", signature(object = "prm_power"), function(object) { print.prm_power(object) })

#' Summarize a power analysis of a randomization test
#'
#' @author kholub
#'
#'@export summary.prm_power
summary.prm_power <- function(object, alpha, alternative, ...) {

  simulation_list <- lapply(object@experiments, function(experiment) {
    prm_summary <- summary(experiment, alternative = alternative)
    
    as.data.frame(prm_summary[, c('p_value', 'actual', 'statistic')])
  })
  
  simulation_results <- as.data.frame(do.call(rbind, simulation_results))
  
  summary_list <- lapply(distint(simulation_results$statistics), function(statistic) {
    statistic_simulation <- simulation_results[simulation_results$statistic == statistic,]
    
    list(statistic = statistic,
         power = sum(statistic_simulation$p_value > alpha),
         average_effect = mean(statistic_simulation$actual))
  })
  
  statistic_summary <- as.data.frame(do.call(rbind, summary_list))
  
  new('summary.prm_power',
      statistic_summary = statistic_summary,
      group_sizes = object@group_sizes,
      alternative = alternative,
      alpha = alpha,
      experiments = nrow(simulation_results))
}

setMethod("summary", signature("prm_power"), summary.prm_power)

#' Access a permutation test power analysis summarization
#'
#' @author kholub
#'
#' @export as.data.frame.summary.prm_power
as.data.frame.summary.prm_power <- function(x, ...) {
  x@statistic_summary
}

setMethod("as.data.frame", signature("summary.prm_power"), as.data.frame.summary.prm_power)

#' Print a power analysis summary of a permutation test
#'
#' @author kholub
#'
#'@export print.summary.prm_power
print.summary.prm_power <- function(x, ...) {
  cat("Power analysis of", x@experiments, "experiments with alternative of", x@alternative, "at significance", as.character(x@alpha), '\n')
  cat("Group sizes:\n")
  print(x@sample_sizes, row.names = FALSE)
  print(x@statistic_summary, row.names = FALSE)
}

setMethod("print", signature(x = "summary.prm_power"), print.summary.prm_power)
setMethod("show", signature(object = "summary.prm_power"), function(object) { print.summary.prm_power(object) })
