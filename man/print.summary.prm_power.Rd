% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power.R
\name{print.summary.prm_power}
\alias{print.summary.prm_power}
\title{Print a power analysis summary of a permutation test}
\usage{
print.summary.prm_power(x, ...)
}
\arguments{
\item{x}{a summary.prm_power object to be printed}
}
\description{
Print a power analysis summary of a permutation test
}
\examples{
gen_data <- function(){ data.frame(x = c(rnorm(50), rnorm(50, 1)), class = rep(c('a', 'b'), each = 50)) }
mcrp <- mcr_power(gen_data, c(mean = function(df){ mean(df[df$class == 'a', 'x']) - mean(df[df$class == 'b', 'x']) }, 
                              median = function(df){ median(df[df$class == 'a', 'x']) - median(df[df$class == 'b', 'x']) }),
    'class')
s <- summary(mcrp, .05, 'two.sided')
print(s)

}
\author{
kholub
}
