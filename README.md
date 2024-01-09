# Functions for Survival Analysis

Zachary R. McCaw <br>
Updated: 2024-01-09


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
##   idx covariates true_event_rate frailty event_time censor_time       time
## 1   1          1               1       1 0.41809511    3.330890 0.41809511
## 2   2          1               1       1 0.09365577    4.263774 0.09365577
## 3   3          1               1       1 0.07381663    5.782096 0.07381663
## 4   4          1               1       1 0.87577447    4.178681 0.87577447
## 5   5          1               1       1 3.87965972    3.258151 3.25815086
## 6   6          1               1       1 0.03210211    3.647888 0.03210211
##   status
## 1      1
## 2      1
## 3      1
## 4      1
## 5      0
## 6      1
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
## 2 0.00556      0      1   100 0.01    0.01      0.0001         0.00141
## 3 0.0248       0      1    99 0.0101  0.0201    0.000202       0.00503
## 4 0.0294       1      0    98 0       0.0201    0.000202       0.00503
## 5 0.0321       0      1    97 0.0103  0.0304    0.000308       0.00981
## 6 0.0355       0      1    96 0.0104  0.0408    0.000417       0.0153 
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
##   tau     rate        se     lower     upper
## 1   1 0.402005 0.0508188 0.3023984 0.4993984
```


```r
# Percentile: median.
SurvUtils::OneSamplePercentiles(data, p = 0.5)
```

```
##   prob      time     lower     upper
## 1  0.5 0.6788454 0.5055086 0.9872552
```


```r
# RMST.
SurvUtils::OneSampleRMST(data, tau = 1.0)
```

```
##   tau       auc         se     lower     upper
## 1   1 0.6219704 0.03884637 0.5458329 0.6981079
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
## 1   0   1 0.345 0.0511
## 2   1   1 0.602 0.0535
## 
## 
## Contrasts:
##   stat   est    se lower upper        p
## 1   rd 0.257 0.074 0.112 0.402 0.000510
## 2   rr 1.740 0.301 1.240 2.450 0.001260
## 3   or 2.870 0.912 1.540 5.350 0.000897
```

### Compare RMSTs

```r
SurvUtils::CompareRMSTs(data, tau = 1.0)
```

```
## Marginal Statistics:
##   tau   auc     se lower upper arm
## 1   1 0.653 0.0361 0.582 0.724   0
## 2   1 0.803 0.0304 0.743 0.862   1
## 
## 
## Contrasts:
##   stat  est     se  lower upper       p
## 1   rd 0.15 0.0472 0.0573 0.242 0.00151
## 2   rr 1.23 0.0824 1.0800 1.400 0.00206
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
  beta_event = c(1.0, -1.0)
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
## 1    1  0.774  0.703 0.071 1.101
## 2    2  0.790  0.699 0.091 1.130
## 3    3  0.821  0.692 0.128 1.185
## 4    4  0.760  0.641 0.119 1.186
## 5    5  0.803  0.667 0.136 1.204
## 6    6  0.822  0.678 0.144 1.213
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
