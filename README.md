# Functions for Survival Analysis

Zachary R. McCaw <br>
Updated: 2024-08-04




```r
suppressPackageStartupMessages({
  library(dplyr)
  library(SurvUtils)
})
```

# Installation


```r
devtools::install_github(repo = "zrmacc/SurvUtils")
```

# Data Generation

Generates survival data with exponential event times and censoring. Optionally, the subject-specific event rate may depend on a set of covariates and/or a gamma-frailty.


```r
data <- SurvUtils::GenData(
  base_event_rate = 1.0,
  censoring_rate = 0.25,
  n = 100,
  tau = 4.0
)
head(data)
```

```
##   idx       time status
## 1   1 1.87244564      0
## 2   2 1.39954366      1
## 3   3 1.62289630      0
## 4   4 1.44909639      1
## 5   5 0.39637045      1
## 6   6 0.07584898      1
```

# Estimation

## One Sample

### Kaplan-Meier

* Tabulates the cumulative hazard and survival functions, along with variance estimates and confidence intervals.


```r
km_tab <- SurvUtils::TabulateKM(data)
head(km_tab)
```

```
## # A tibble: 6 × 13
##      time censor events   nar    haz cum_haz cum_haz_var cum_haz_lower
##     <dbl>  <dbl>  <dbl> <dbl>  <dbl>   <dbl>       <dbl>         <dbl>
## 1 0            0      0   100 0       0         0              0      
## 2 0.00196      0      1   100 0.01    0.01      0.0001         0.00141
## 3 0.00416      0      1    99 0.0101  0.0201    0.000202       0.00503
## 4 0.0266       1      0    98 0       0.0201    0.000202       0.00503
## 5 0.0321       0      1    97 0.0103  0.0304    0.000308       0.00981
## 6 0.0411       0      1    96 0.0104  0.0408    0.000417       0.0153 
## # ℹ 5 more variables: cum_haz_upper <dbl>, surv <dbl>, surv_var <dbl>,
## #   surv_lower <dbl>, surv_upper <dbl>
```

### Event Rate, Percentile, Restricted Mean Survival

* Calculate the event rate at a point in time.


```r
# Rate.
SurvUtils::OneSampleRates(data, tau = 1.0)
```

```
##   tau      rate         se     lower     upper
## 1   1 0.3061041 0.04771619 0.2160165 0.4006978
```


```r
# Percentile: median.
SurvUtils::OneSamplePercentiles(data, p = 0.5)
```

```
##   prob      time     lower     upper
## 1  0.5 0.5514958 0.4060927 0.8082641
```


```r
# RMST.
SurvUtils::OneSampleRMST(data, tau = 1.0)
```

```
##   tau       auc         se     lower     upper
## 1   1 0.5896377 0.03643433 0.5182277 0.6610477
```

## Two Sample

### Generate Data

```r
data0 <- SurvUtils::GenData(
  base_event_rate = 1.0,
  censoring_rate = 0.25,
  n = 100,
  tau = 4.0
)
data0$arm <- 0

data1 <- SurvUtils::GenData(
  base_event_rate = 0.5,
  censoring_rate = 0.25,
  n = 100,
  tau = 4.0
)
data1$arm <- 1
data <- rbind(data0, data1)
```

### Compare Rates

```r
SurvUtils::CompareRates(data, tau = 1.0)
```

```
## Marginal Statistics:
##   arm tau  rate     se
## 1   0   1 0.320 0.0487
## 2   1   1 0.683 0.0489
## 
## 
## Contrasts:
##   stat   est    se lower upper        p
## 1   rd 0.362 0.069 0.227 0.498 1.51e-07
## 2   rr 2.130 0.358 1.530 2.960 6.79e-06
## 3   or 4.560 1.450 2.450 8.510 1.77e-06
```

### Compare RMSTs

```r
SurvUtils::CompareRMSTs(data, tau = 1.0)
```

```
## Marginal Statistics:
##   tau   auc     se lower upper arm
## 1   1 0.593 0.0373 0.520 0.666   0
## 2   1 0.807 0.0332 0.742 0.872   1
## 
## 
## Contrasts:
##   stat   est    se lower upper        p
## 1   rd 0.214 0.050 0.116 0.312 1.92e-05
## 2   rr 1.360 0.102 1.170 1.580 4.33e-05
```

### Compare Cox Models
Compare the predictive performance of Cox models based on different sets of covariates with respect to their c-statistics on held-out data via k-fold cross validation.

```r
# Simulate data.
n <- 1000
x1 <- rnorm(n)
x2 <- rnorm(n)
data <- SurvUtils::GenData(
  covariates = cbind(x1, x2),
  beta_event = c(1.0, -1.0),
  simple = FALSE
)

# Evaluate.
eval <- CompreCoxCstat(
  status = data$status,
  time = data$time,
  x1 = data %>% dplyr::select(x1, x2),
  x2 = data %>% dplyr::select(x1)
)

head(round(eval, digits = 3))
```

```
##   fold cstat1 cstat2  diff ratio
## 1    1  0.766  0.674 0.092 1.137
## 2    2  0.732  0.641 0.091 1.142
## 3    3  0.795  0.759 0.037 1.048
## 4    4  0.793  0.710 0.084 1.118
## 5    5  0.829  0.700 0.129 1.184
## 6    6  0.715  0.651 0.064 1.099
```

# Inference

For a tutorial on influence functions and the perturbation bootstrap, see [this vignette](https://github.com/zrmacc/SurvUtils/blob/master/vignettes/perturbation_bootstrap.pdf).

# Plotting


```r
# Generate data.
arm1 <- SurvUtils::GenData(base_event_rate = 0.8)
arm1$arm <- 1
arm0 <- SurvUtils::GenData(base_event_rate = 1.0)
arm0$arm <- 0
data <- rbind(arm1, arm0)
```

## One Sample

### Standard Kaplan-Meier


```r
x_breaks <- seq(from = 0.0, to = 4.0, by = 0.50)
data0 <- data %>% dplyr::filter(arm == 0)
fit0 <- Temporal::FitParaSurv(data0)  # Optional parametric fit. 
q_km <- SurvUtils::PlotOneSampleKM(data0, fit = fit0, x_breaks = x_breaks, x_max = 4)
q_nar <- SurvUtils::PlotOneSampleNARs(data0, x_breaks = x_breaks, x_max = 4)
cowplot::plot_grid(
  plotlist = list(q_km, q_nar),
  align = "v",
  axis = "l",
  ncol = 1,
  rel_heights = c(3, 1)
)
```

<img src="README_files/figure-html/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

### AUC


```r
x_breaks <- seq(from = 0.0, to = 4.0, by = 0.50)
data0 <- data %>% dplyr::filter(arm == 0)
q_auc <- SurvUtils::PlotOneSampleAUC(data0, x_breaks = x_breaks, x_max = 4, tau = 3)
q_nar <- SurvUtils::PlotOneSampleNARs(data0, x_breaks = x_breaks, x_max = 4)
cowplot::plot_grid(
  plotlist = list(q_auc, q_nar),
  align = "v",
  axis = "l",
  ncol = 1,
  rel_heights = c(3, 1)
)
```

<img src="README_files/figure-html/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

## Two Sample


```r
x_breaks <- seq(from = 0.0, to = 4.0, by = 0.50)
contrast <- Temporal::CompParaSurv(data)  # Optional parametric fit. 
q_km <- SurvUtils::PlotTwoSampleKM(data, contrast = contrast, x_breaks = x_breaks, x_max = 4)
q_nar <- SurvUtils::PlotTwoSampleNARs(data, x_breaks = x_breaks, x_max = 4)
cowplot::plot_grid(
  plotlist = list(q_km, q_nar),
  align = "v",
  axis = "l",
  ncol = 1,
  rel_heights = c(3, 1)
)
```

<img src="README_files/figure-html/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />
