# sampling-investigations

This directory contains results from sampling and simulations using R programming

1. THe Reimann Zeta function is noisy, complex and infinite. Using a conjugate pair version of the function results in a smooth function with FM modulation behaviour, independent of the distance from the critical line. The smooth envelope function will help allow better estimates of the asymptotic behaviour of the Reimann Zeta function.
 
2. Quantile regression model coefficient variance estimation, homoscedastic iid case with collinearity: Tested and implemented covariance matrix terms for model intercept coefficient CI bivariate case. Collinearity between the explanatory variables, results in variance inflation of slope CIs of collinear variables. The variable inflation factor is present in intercept CI estimator but is finely balanced by numerator covariance term to leave intercept CI relatively insensitive to collinearity. Improved bivariate performance is found to occur compared to earlier paper examples. Old r code also had mean(x) which needs to be mean(c(x)) to ensure mean is vector calculation.
 
3. Quantile regression model coefficient variance estimation, linear expanding horn heteroscedastic id case: Tested and adapted quantile regression model coefficient CI estimators to linear expanding horn heteroscedastic id case. Used auxiliary regression modelling of the heteroscedasticity in the quantile regression model residuals. The slope CI estimators adapted with no change for univariate case, again using max half CIs. The terms in the standard regression formula for intercept variance required use of median(x,w_auxreg) for constant term, weighted medians for slope terms and use of mean half CI for adjusted slope variance (in the intercept formula). Nominal 95% coverage for n=1000 cases, weaker for smaller samples. Covariance term correction needs further work. The minimum point estimate CI occurs for the auxiliary regression weighted median and is robust to sample size.

4. Quantile regression model coefficient variance estimation, homoscedastic iid case: Tested and adapted sample quantile CI estimators to the task of quantile regression model coefficient variance estimation. Needed to use max half CIs to increase coverage and worked out how to adapt standard regression formula for slope variance (unchanged) and intercept variance (obvious change). Nominal 95% coverage for homoscedastic iid cases.

5. Sample quantile estimation: Suitably transforming the reference frame for quantile variance estimation can improve smoothed quantile estimating function variance estimates of sample quantile standard errors and get good agreement with bootstrap estimates of the linear programming (unsmoothed) quantile estimating function solution. The approach has potential for quantile regression where variances estimates of Y-bX are calculated rather than X as calculated in the first paper "empirical variance distribution ..."

6. Median estimation: Suitably transforming the reference frame for jackknife variance estimation can improve jackknife variance estimates of median confidence intervals and get good agreement with bootstrap estimates.  

7. Synthetic Research: Investigation of suitable A/B sampling sizes based on public data about Zulily customer base. This report is a synthetic construction of a customer base using public data describing high level numbers for the size of the customer base and the relative behaviour (purchases per year) of new and old customers. The intent was to see if internal structure of large customner bases could impact on confidence intervals of A/B samplimng estimates. As demonstrated stratified sampling can decrease the confidence interval of estimates. 

the end
