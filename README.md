
<!-- README.md is generated from README.Rmd. Please edit that file -->

# didunit

<!-- badges: start -->
<!-- badges: end -->

Calculates differences-in-differences treatment effects for every
treated unit. The more granular estimator allows for diagnostics,
aggregation on custom variables, and imposes composition balance even in
unbalanced panels. It extends the functionality of the
[did](https://bcallaway11.github.io/did/) package in applications where
within treatment-group heterogeneity is high and can generalize to
non-binary treatments to produce dose-response functions.

## Installation

You can install the development version of didunit like so:

You can install the development version of didwrappers

``` r
#install.packages("devtools")
#devtools::install_github("ransiw/didunit", build_vignettes = TRUE)
```

## Example

This is a basic example generating code with a sample dataset:

``` r
library(didunit)
```

Simulate sample data

``` r
simdata = sim_data()
head(simdata[,c("unit","time","treatg","dosage","y")])
#>   unit time treatg dosage          y
#> 1    1    5     10      1 -0.5304577
#> 2    1    6     10      1  2.6677010
#> 3    1    7     10      1  4.7607068
#> 4    1    8     10      1  5.2621771
#> 5    1    9     10      1  7.2237843
#> 6    1   10     10      1 19.6710919
```

Run the first-step `att_it()` function to get estimates for each treated
unit at each time and tabulate a sample of these.

``` r
attobject = att_it(yname = "y", tname = "time", gname = "treatg", idname ="unit", data = simdata)
attdf = attit_table(attobject)
head(attdf)
#>   id group  t         att       lci       uci
#> 1  1    10  6  1.46091865 -1.106683  4.000371
#> 2  1    10  7  0.02580734 -2.288994  2.650115
#> 3  1    10  8 -1.00526502 -3.675418  1.265061
#> 4  1    10  9 -0.15021294 -3.141658  2.538498
#> 5  1    10 10 10.07057382  7.505619 12.508717
#> 6  1    10 11  7.59592035  5.025431 10.105523
```

Now aggregate all post-treatment effects of a treated unit.

``` r
agtobject = aggite(attobject,type="unit")
aggite_table(agtobject)
#>    egt   att.egt   lci.egt  uci.egt
#> 1    1 10.268125  7.962492 12.57376
#> 2    2 10.860571  9.007725 12.71342
#> 3    3 10.459844  8.160576 12.75911
#> 4    4 11.244004  9.402217 13.08579
#> 5    5 11.549739  9.238233 13.86124
#> 6    6  9.453484  7.596821 11.31015
#> 7    7  9.251990  6.947623 11.55636
#> 8    8  9.892762  8.041271 11.74425
#> 9    9  9.564710  7.243194 11.88623
#> 10  10 12.168149 10.316889 14.01941
```

Sample a data frame with different dosage amounts and aggregate to the
unit level.

``` r
attobject = att_it(yname = "y", tname = "time", gname = "treatg", idname ="unit", customnames = "dosage", data = sim_data(dosage = rep(c(1,2),each=5)))
agtobject = aggite(attobject,type="unit")
aggite_table(agtobject)
#>    egt   att.egt   lci.egt  uci.egt
#> 1    1  9.445631  7.337315 11.55395
#> 2    2 11.248110  8.555919 13.94030
#> 3    3  8.263400  6.170324 10.35648
#> 4    4 10.157019  7.468782 12.84526
#> 5    5 10.619653  8.530252 12.70905
#> 6    6 19.591504 16.915919 22.26709
#> 7    7 21.582226 19.476365 23.68809
#> 8    8 21.340354 18.657675 24.02303
#> 9    9 20.016970 17.919413 22.11453
#> 10  10 20.142061 17.468207 22.81592
```

Aggregate to the dosage level by specifying dosage as the `type`.

``` r
agtobject = aggite(attobject,type="dosage")
aggite_table(agtobject)
#>   egt   att.egt   lci.egt  uci.egt
#> 1   1  9.867205  8.815234 10.91918
#> 2   2 20.564065 19.462440 21.66569
```
