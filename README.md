# multiStockassessment

The `multiStockassessment` package extends the
[`stockassessment`](https://github.com/fishfollower/SAM) to allow linked
stocks (See e.g. Albertsen, C. M., Nielsen, A., and Thygesen, U. H.
(2018) Connecting single-stock assessment models through correlated
survival. *ICES Journal of Marine Science*, 75(1), 235-244. doi:
[10.1093/icesjms/fsx114](https://dx.doi.org/10.1093/icesjms/fsx114)).
The package can be installed with
`devtools`:

``` r
devtools::install_github("calbertsen/multi_SAM", subdir = "multiStockassessment")
```

## Short example

To use the `multiStockassessment` package, the individual stocks must
first be fitted with the `stockassessment` package and combined to a
`samset` using the `c(...)` function.

### Preparing data using stockassessment

The code below loads the North Sea cod data from the stockassessment
package, fits it, simulates two new data sets and fit those
individually.

``` r
library(stockassessment)

data(nscodData)
data(nscodConf)
data(nscodParameters)

set.seed(9876)
fit <- sam.fit(nscodData, nscodConf, nscodParameters,sim.condRE = FALSE)

datSim <- simulate(fit,nsim=2)
fitSim <- do.call("c",lapply(datSim,function(x)sam.fit(x,nscodConf, nscodParameters)))
```

### Creating correlation structure

``` r
library(multiStockassessment)
```

After the simulated data have been fitted and combined by the `c`
function, the `multiStockassessment` package can be used to fit a model
with correlated survival.

The `suggestCorStructure` function can be used to set up the correlation
matrix. The `nAgeClose` argument can be used to construct the banded
structures used in the paper. The result is a logical symmetric matrix
indicating whether a correlation should be fixed at zero. The dimension
of the matrix is the number of stocks times the number of ages in the
data.

For instance, `suggestCorStructure(fitSim,nAgeClose=2)` creates a matrix
where ages less than 2 appart are correlated between the stocks.
Likewise,

``` r
cs <- suggestCorStructure(fitSim,nAgeClose=2)
cs
```

    ##        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11]
    ##  [1,]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE  TRUE  TRUE  TRUE
    ##  [2,]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE  TRUE  TRUE
    ##  [3,]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE  TRUE
    ##  [4,]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE
    ##  [5,]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE
    ##  [6,]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE
    ##  [7,] FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ##  [8,] FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ##  [9,]  TRUE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [10,]  TRUE  TRUE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [11,]  TRUE  TRUE  TRUE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [12,]  TRUE  TRUE  TRUE  TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE
    ##       [,12]
    ##  [1,]  TRUE
    ##  [2,]  TRUE
    ##  [3,]  TRUE
    ##  [4,]  TRUE
    ##  [5,] FALSE
    ##  [6,] FALSE
    ##  [7,]  TRUE
    ##  [8,]  TRUE
    ##  [9,]  TRUE
    ## [10,]  TRUE
    ## [11,]  TRUE
    ## [12,]  TRUE

The `suggestCorStructure` function has several options to configure the
correlation structure, and the result can be modified by hand for
complete freedom.

### Fitting the model

The model is fitted by the `multisam.fit` function. The function
requires a `samset` and a correlaiton structure. By default, the
correlation structure given models the partial correlations between
cohorts. This can be changed to regular correlations by setting
`usePartialCors =
    FALSE`.

``` r
obj <- multisam.fit(fitSim,cs)
```

``` r
obj
```

    ## Multi-SAM model with 2 stocks: log likelihood is -259.3252. Convergence OK.

### Investigating the result

Most plotting and table functions from the `stockassessment` package
have a corresponding version in `multiStockassessment`.

The fitted correlation structure can be plotted by

``` r
corplot(obj)
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Other plotting functions implemented from the `stockassessment` package
are

``` r
fbarplot(obj)
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
ssbplot(obj)
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
ssbplot(obj)
```

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
tsbplot(obj)
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
recplot(obj)
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
catchplot(obj)
```

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
srplot(obj)
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
fitplot(obj,1)
```

![](README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Note that the `fitplot` function above requires the stock to be
specified.

``` r
fitplot(obj,2)
```

![](README_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

### Comparing models

Models can be compared with the `modeltable` function. For instance, the
model above can be compared to individual assessments by

``` r
cs2 <- suggestCorStructure(fitSim,nAgeClose=0)
obj2 <- multisam.fit(fitSim,cs2)
```

``` r
modeltable(c(obj,obj2))
```

    ##       log(L) #par      AIC Pval( M1 -> M2 )
    ## M1 -259.3252   84 686.6503               NA
    ## M2 -265.0501   68 666.1002        0.7809027
