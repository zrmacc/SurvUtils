# Functions for Survival Analysis

Zachary R. McCaw <br>
Updated: 2022-03-19



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
## 1   1          1               1       1  0.2689345 10.79424136 0.26893452
## 2   2          1               1       1  0.2141906  0.65922102 0.21419058
## 3   3          1               1       1  1.1131212  0.42156700 0.42156700
## 4   4          1               1       1  0.2849961  0.09927574 0.09927574
## 5   5          1               1       1  0.6185080  4.85548160 0.61850796
## 6   6          1               1       1  1.3814466  5.62214240 1.38144665
##   status
## 1      1
## 2      1
## 3      0
## 4      0
## 5      1
## 6      1
```

# Tabulate Kaplan-Meier

* Tabulates the cumulative hazard and survival functions, along with variance estimates and confidence intervals.


```r
km_tab <- SurvUtils::TabulateKM(data)
head(km_tab)
```

```
## # A tibble: 6 × 13
##      time censor events   nar    haz cum_haz var_cum_haz cum_haz_lower
##     <dbl>  <dbl>  <dbl> <dbl>  <dbl>   <dbl>       <dbl>         <dbl>
## 1 0            0      0   100 0       0         0              0      
## 2 0.00144      1      0   100 0       0         0              0      
## 3 0.00954      0      1    99 0.0101  0.0101    0.000102       0.00142
## 4 0.0199       0      1    98 0.0102  0.0203    0.000206       0.00508
## 5 0.0217       0      1    97 0.0103  0.0306    0.000312       0.00987
## 6 0.0219       0      1    96 0.0104  0.0410    0.000421       0.0154 
## # … with 5 more variables: cum_haz_upper <dbl>, surv <dbl>, var_surv <dbl>,
## #   surv_lower <dbl>, surv_upper <dbl>
```

