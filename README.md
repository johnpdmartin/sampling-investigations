# sampling-investigations

This directory contains results from sampling and simulations using R programming

1. Quantile regression model coefficient variance estimation: Tested and adapted sample quantile CI estimators to the task of quantile regression model coefficient variance estimation. Needed to use max half CIs to increase coverage and worked out how to adapt standard regression formula for slope variance (unchanged) and intercept variance (obvious change). Nominal 95% coverage for homoscedastic iid cases.

2. Sample quantile estimation: Suitably transforming the reference frame for quantile variance estimation can improve smoothed quantile estimating function variance estimates of sample quantile standard errors and get good agreement with bootstrap estimates of the linear programming (unsmoothed) quantile estimating function solution. The approach has potential for quantile regression where variances estimates of Y-bX are calculated rather than X as calculated in the first paper "empirical variance distribution ..."

3. Median estimation: Suitably transforming the reference frame for jackknife variance estimation can improve jackknife variance estimates of median confidence intervals and get good agreement with bootstrap estimates.  

4. Synthetic Research: Investigation of suitable A/B sampling sizes based on public data about Zulily customer base. This report is a synthetic construction of a customer base using public data describing high level numbers for the size of the customer base and the relative behaviour (purchases per year) of new and old customers. The intent was to see if internal structure of large customner bases could impact on confidence intervals of A/B samplimng estimates. As demonstrated stratified sampling can decrease the confidence interval of estimates. 

the end
