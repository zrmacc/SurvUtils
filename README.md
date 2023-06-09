# Functions for Survival Analysis

Zachary R. McCaw <br>
Updated: 2023-06-09


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
## 1   1          1               1       1 1.25882533   1.7383592 1.25882533
## 2   2          1               1       1 3.77287448   0.2994819 0.29948189
## 3   3          1               1       1 0.86277663   9.9509909 0.86277663
## 4   4          1               1       1 0.05351559   7.3693073 0.05351559
## 5   5          1               1       1 0.33578540   1.1028062 0.33578540
## 6   6          1               1       1 0.74077214  16.8707338 0.74077214
##   status
## 1      1
## 2      0
## 3      1
## 4      1
## 5      1
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
## 2 0.00619      1      0   100 0       0         0              0      
## 3 0.0103       0      1    99 0.0101  0.0101    0.000102       0.00142
## 4 0.0132       0      1    98 0.0102  0.0203    0.000206       0.00508
## 5 0.0317       1      0    97 0       0.0203    0.000206       0.00508
## 6 0.0331       0      1    96 0.0104  0.0307    0.000315       0.00991
## # ℹ 5 more variables: cum_haz_upper <dbl>, surv <dbl>, surv_var <dbl>,
## #   surv_lower <dbl>, surv_upper <dbl>
```

#### Influence function

The influence function expansion of the Kaplan-Meier estimator is:

$$
\sqrt{n}\big\{ \hat{S}(t) - S(t) \big\} = \frac{1}{\sqrt{n}}\sum_{i=1}^{n}\psi_{i}(t) + o_{p}(1)
$$
where:

$$
\psi_{i}(t) = -S(t) \int_{0}^{t}\frac{dM_{i}(u)}{n^{-1} Y(u)}
$$
and $dM_{i}$ is the counting-process martingale:

$$
dM_{i}(u) = dN_{i}(u) - Y_{i}(u)dH(u)
$$


```r
n <- 1000
data <- GenData(n = n)
tau <- 1.0

# Influence function calculation.
influence <- KMInfluence(data, tau = tau)

# Estimated variance of the KM estimator at time tau.
inf_var <- mean(influence^2) / n
inf_se <- sqrt(inf_var)
print(round(inf_se, digits = 4))
```

```
## [1] 0.0163
```

### Event Rate, Percentile, Restricted Mean Survival

* Calculate the event rate at a point in time.


```r
# Rate.
SurvUtils::OneSampleRates(data, tau = 1.0)
```

```
##   tau      rate         se     lower     upper
## 1   1 0.3559764 0.01628601 0.3241441 0.3879014
```


```r
# Percentile: median.
SurvUtils::OneSamplePercentiles(data, p = 0.5)
```

```
##   prob     time     lower     upper
## 1  0.5 0.678371 0.6172791 0.7609683
```


```r
# RMST.
SurvUtils::OneSampleRMST(data, tau = 1.0)
```

```
##   tau       auc         se     lower     upper
## 1   1 0.6248738 0.01173439 0.6018748 0.6478728
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
## 1   0   1 0.366 0.0508
## 2   1   1 0.635 0.0508
## 
## 
## Contrasts:
##   stat   est     se lower upper        p
## 1   rd 0.268 0.0718 0.128 0.409 0.000186
## 2   rr 1.730 0.2770 1.270 2.370 0.000592
## 3   or 3.010 0.9300 1.640 5.510 0.000377
```


```r
SurvUtils::CompareRMSTs(data, tau = 1.0)
```

```
## Marginal Statistics:
##   tau   auc     se lower upper arm
## 1   1 0.582 0.0369  0.51 0.654   0
## 2   1 0.760 0.0303  0.70 0.819   1
## 
## 
## Contrasts:
##   stat   est     se  lower upper        p
## 1   rd 0.178 0.0477 0.0845 0.272 0.000190
## 2   rr 1.310 0.0978 1.1300 1.510 0.000361
```

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
