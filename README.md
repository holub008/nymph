# nymph
Randomization tests for nonparametric inference

## Installation
Hopefully this project makes it to CRAN. Until then:

```R
devtools::install_git('https://github.com/holub008/nymph')
```

## Purpose
Randomization methods provide a powerful toolset for performing inference on a wide range of statistics with few distributional constraints on the data. I believe these methods are underutilized, partially as a result of poor library support; existing R libraries (e.g. [coin](https://cran.r-project.org/web/packages/coin/index.html) & [perm](https://cran.r-project.org/web/packages/perm/index.html)) provide a narrow range of functionality or obfuscate the underlying simplicity of the methods from the user. nymph aims to provide the practitioner a simple & generic interface to this class of methods.

## Examples
WIP

## Implementation
  * No dependencies
    * Portability & ease of install
    * Functionality should be transparent to the user
    * Suggests the [parallel](http://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf), which comes standard with R >= 2.14, for accelerated computations
  * [S4 object system](https://stat.ethz.ch/R-manual/R-devel/library/methods/html/Introduction.html)
    * Ensure consistency & validity across a range of different tests
    * Fewer redundant objects for the user to comprehend
  * [Non-standard evaluation](http://developer.r-project.org/nonstandard-eval.pdf) wrappers for common use cases
    * Reduce boilerplate for interactive use
