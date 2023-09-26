# Functions for Survival Analysis

Zachary R. McCaw <br>
Updated: 2023-09-26




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
## 1   1          1               1       1 0.15032178   0.8281659 0.15032178
## 2   2          1               1       1 0.08965564  10.5189404 0.08965564
## 3   3          1               1       1 0.14340693   2.5481145 0.14340693
## 4   4          1               1       1 0.71691154   7.6195282 0.71691154
## 5   5          1               1       1 3.02119978   5.6989281 3.02119978
## 6   6          1               1       1 0.55638190   3.0264381 0.55638190
##   status
## 1      1
## 2      1
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
##     time censor events   nar    haz cum_haz cum_haz_var cum_haz_lower
##    <dbl>  <dbl>  <dbl> <dbl>  <dbl>   <dbl>       <dbl>         <dbl>
## 1 0           0      0   100 0       0         0              0      
## 2 0.0253      0      1   100 0.01    0.01      0.0001         0.00141
## 3 0.0454      0      1    99 0.0101  0.0201    0.000202       0.00503
## 4 0.0512      0      1    98 0.0102  0.0303    0.000306       0.00977
## 5 0.0859      0      1    97 0.0103  0.0406    0.000412       0.0152 
## 6 0.0897      0      1    96 0.0104  0.0510    0.000521       0.0212 
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
## 1   1 0.3591066 0.01630011 0.3272344 0.3910474
```


```r
# Percentile: median.
SurvUtils::OneSamplePercentiles(data, p = 0.5)
```

```
##   prob     time     lower     upper
## 1  0.5 0.671038 0.6238443 0.7326507
```


```r
# RMST.
SurvUtils::OneSampleRMST(data, tau = 1.0)
```

```
##   tau       auc         se     lower     upper
## 1   1 0.6249941 0.01167324 0.6021149 0.6478732
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
## 1   0   1 0.401 0.0522
## 2   1   1 0.558 0.0527
## 
## 
## Contrasts:
##   stat   est     se  lower upper      p
## 1   rd 0.158 0.0742 0.0125 0.303 0.0333
## 2   rr 1.390 0.2240 1.0200 1.910 0.0389
## 3   or 1.890 0.5770 1.0400 3.440 0.0364
```

### Compare RMSTs

```r
SurvUtils::CompareRMSTs(data, tau = 1.0)
```

```
## Marginal Statistics:
##   tau   auc     se lower upper arm
## 1   1 0.625 0.0377 0.551 0.699   0
## 2   1 0.729 0.0352 0.660 0.798   1
## 
## 
## Contrasts:
##   stat   est     se   lower upper      p
## 1   rd 0.104 0.0516 0.00254 0.205 0.0445
## 2   rr 1.170 0.0901 1.00000 1.360 0.0471
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
## 1    1  0.772  0.641 0.130 1.203
## 2    2  0.796  0.707 0.089 1.126
## 3    3  0.826  0.748 0.079 1.105
## 4    4  0.776  0.640 0.136 1.213
## 5    5  0.723  0.602 0.120 1.200
## 6    6  0.803  0.744 0.059 1.079
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

<img src="README_files/figure-html/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

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

<img src="README_files/figure-html/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

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

<img src="README_files/figure-html/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />
