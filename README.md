# Functions for Survival Analysis

Zachary R. McCaw <br>
Updated: 2022-10-09


```r
suppressPackageStartupMessages({
  library(dplyr)
  library(SurvUtils)
})
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
##   idx covariates true_event_rate frailty event_time censor_time      time
## 1   1          1               1       1  1.9405879   6.1577900 1.9405879
## 2   2          1               1       1  0.4403402   4.4846209 0.4403402
## 3   3          1               1       1  0.7509104   0.5108565 0.5108565
## 4   4          1               1       1  0.1005392   1.0944711 0.1005392
## 5   5          1               1       1  0.3138042   0.6715436 0.3138042
## 6   6          1               1       1  0.6166296  12.2708344 0.6166296
##   status
## 1      1
## 2      1
## 3      0
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
##     time censor events   nar    haz cum_haz cum_haz_var cum_haz_…¹ cum_h…²  surv
##    <dbl>  <dbl>  <dbl> <dbl>  <dbl>   <dbl>       <dbl>      <dbl>   <dbl> <dbl>
## 1 0           0      0   100 0       0         0           0        0      1    
## 2 0.0442      0      1   100 0.01    0.01      0.0001      0.00141  0.0710 0.99 
## 3 0.0583      1      0    99 0       0.01      0.0001      0.00141  0.0710 0.99 
## 4 0.0693      0      1    98 0.0102  0.0202    0.000204    0.00505  0.0808 0.980
## 5 0.0723      1      0    97 0       0.0202    0.000204    0.00505  0.0808 0.980
## 6 0.0907      1      0    96 0       0.0202    0.000204    0.00505  0.0808 0.980
## # … with 3 more variables: surv_var <dbl>, surv_lower <dbl>, surv_upper <dbl>,
## #   and abbreviated variable names ¹​cum_haz_lower, ²​cum_haz_upper
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
## [1] 0.0165
```

### Event Rate, Percentile, Restricted Mean Survival

* Calculate the event rate at a point in time.


```r
# Rate.
SurvUtils::OneSampleRates(data, tau = 1.0)
```

```
##   tau      rate         se    lower    upper
## 1   1 0.3889673 0.01654321 0.356504 0.421273
```


```r
# Percentile: median.
SurvUtils::OneSamplePercentiles(data, p = 0.5)
```

```
##   prob     time     lower     upper
## 1  0.5 0.724318 0.6666131 0.8084555
```


```r
# RMST.
SurvUtils::OneSampleRMST(data, tau = 1.0)
```

```
##   tau       auc         se     lower    upper
## 1   1 0.6469971 0.01167413 0.6241162 0.669878
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
## 1   0   1 0.431 0.0523
## 2   1   1 0.579 0.0522
## 
## 
## Contrasts:
##   stat   est     se   lower upper      p
## 1   rd 0.148 0.0739 0.00325 0.293 0.0451
## 2   rr 1.340 0.2030 0.99900 1.810 0.0507
## 3   or 1.820 0.5490 1.00000 3.280 0.0483
```


```r
SurvUtils::CompareRMSTs(data, tau = 1.0)
```

```
## Marginal Statistics:
##   tau   auc     se lower upper arm
## 1   1 0.688 0.0352 0.619 0.757   0
## 2   1 0.727 0.0331 0.662 0.791   1
## 
## 
## Contrasts:
##   stat    est     se   lower upper     p
## 1   rd 0.0382 0.0483 -0.0565 0.133 0.430
## 2   rr 1.0600 0.0723  0.9230 1.210 0.431
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
fit0 <- Temporal::FitParaSurv(data0)
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

<img src="README_files/figure-html/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

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

<img src="README_files/figure-html/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

## Two Sample


```r
x_breaks <- seq(from = 0.0, to = 4.0, by = 0.50)
contrast <- Temporal::CompParaSurv(data)
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

<img src="README_files/figure-html/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />
