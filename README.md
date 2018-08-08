# nymph
Randomization tests for nonparametric inference

## Installation
```R
devtools::install_git('https://github.com/holub008/nymph')
```

## Purpose
Randomization methods provide a powerful toolset for performing inference on a wide range of statistics with few distributional constraints on the data. I believe these methods are underutilized, partially as a result of poor library support; existing R packages (e.g. [coin](https://cran.r-project.org/web/packages/coin/index.html) & [perm](https://cran.r-project.org/web/packages/perm/index.html)) provide a narrow range of functionality or obfuscate the underlying simplicity of the methods from the user. nymph aims to provide the practitioner a simple & generic interface to this class of methods.

## Examples

### Inference
Using the iris dataset, we'd like to investigate if the ratio of length to width of Virginica flower petals are different from that of Virginica sepals. We use the _mcrd_ function to compute ratio differences between petals and sepals.
```R
library(dplyr)

petal_measurements <- iris %>% 
    filter(Species == 'virginica') %>% 
    mutate(len = Petal.Length,
           width = Petal.Width) %>%
    select(len, width)
sepal_measurements <- iris %>% 
    filter(Species == 'virginica') %>% 
    mutate(len = Sepal.Length,
           width = Sepal.Width) %>%
    select(len, width)

set.seed(55414)
mcrt <- mcrd_test(petal_measurements, sepal_measurements, 
                  lw_ratio = mean(len / width),
                  length_proportion = mean(len / (len + width)))
summary(mcrt)
```
with results:
```
Permutation test of 1000 permutations against alternative of two.sided at significance 0.95 
         statistic    ci_lower   ci_upper     actual p_value
          lw_ratio -0.16575289 0.17387946 0.55020960       0
 length_proportion -0.01267495 0.01315527 0.04388741       0
```
Quite convincing that a difference exists! For a visualization:
```R
plot(mcrt)
```

### Power Analysis
Here we perform a power analysis of a contrived experiment - a two treatment ('a' & 'b') experiment with a minimally impactful effect size of 1 against otherwise standard normal populations.
```R
gen_data_s1 <- function(){ data.frame(x = rnorm(50), class = 'a') }
gen_data_s2 <- function(){ data.frame(x = rnorm(50, 1), class = 'b') }
mcrp <- mcrd_power(gen_data_s1, gen_data_s2, mean = mean(x), median = median(x), 
                    test_trials = 1e2)
summary(mcrp)
```
With result:
```
Power analysis of 100 experiments with alternative of two.sided at significance 0.05 
Group sizes:
 sample size
      1   50
      2   50
 statistic power average_effect
      mean     1     -0.9829345
    median  0.99     -0.9878372
```

We can also visualize the distribution of p-values (from repeated simulation of the experiment):
```R
plot(mcrp, statistic = 'median', alternative = 'two.sided', 
     alpha = NULL)
```

Or we can visualize how power varies across desired inferential FPRs:
```R
plot(mcrp, statistic = 'median', alternative = 'two.sided', 
     alpha = seq(.01, .2, by = .01))
```

## Implementation
  * No dependencies
    * Portability & ease of install
    * Functionality should be transparent to the user
    * Suggests the [parallel](http://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf) package, which comes standard with R >= 2.14, for accelerated computations
  * [S4 object system](https://stat.ethz.ch/R-manual/R-devel/library/methods/html/Introduction.html)
    * Ensure consistency & validity across a range of tests
    * Fewer redundant objects for the user to comprehend
  * [Non-standard evaluation](http://developer.r-project.org/nonstandard-eval.pdf) wrappers for common use cases
    * Reduce boilerplate for interactive use
