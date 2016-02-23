library(Runuran);library(quantreg);library(Bolstad);library(scatterplot3d)
# START MANUAL SETTING 1
# set heteroscedasticity and population distribution slope and noise of interest
sampsize <- 1000;replicates <- 1000;p_con <- 1;multi <- 0.5;lambda <- 15; multi2 <- 10
adj_het <- "adjHet" # "adjHet"   "none"
method_df <- "perc_scale" # "perc_scale" # "org_scale" #  
# initialise some important arrays
x_samp_org <- matrix(rep(1,sampsize*replicates),nrow=replicates,ncol=sampsize);x_samp <- x_samp_org
x_samp2_org <- x_samp_org; x_samp2 <- x_samp2_org
y_samp_org <- matrix(rep(1,sampsize*replicates),nrow=replicates,ncol=sampsize);y_samp <- y_samp_org
med_pt_est_x <- matrix(rep(1,replicates),nrow=replicates,ncol=1);med_pt_est_x2 <- med_pt_est_x
med_pt_est_x_wgt <- matrix(rep(1,replicates),nrow=replicates,ncol=1);med_pt_est_x2_wgt <- med_pt_est_x_wgt
output <- matrix(rep(1,40*15),nrow=15,ncol=40)
#filename <- "oneterm_Hetero_x_lin_1000reps_p15_dfadj_in_perc_scale_aux_reg" #  _perc_scale" # _org_scale" #
filename <- "v3runif_intrighttwoterms_m360_0_xlin_p15_dfadj_in_perc_scale_aux_reg" #  _perc_scale" # _org_scale" #
# filename <- "intassessv3engelproxy_intrightwmm3vers_trimadjmedHetero_x_lin_1000reps_p15_dfadj_in_perc_scale_aux_reg" #  _perc_scale" # _org_scale" #

# setting up original datasets and any transformation for quantile regression
set.seed(592)  
for (reps_setup in 1:replicates) {
#0#-1500#-1800#-5000
    x_samp_org[reps_setup,] <- runif(sampsize,-0,360)+(-360)#exp(rnorm(sampsize,7.5,.5)^1)+(-5000)#runif(sampsize,40,400)#runif(sampsize,-180,180)#exp(rnorm(sampsize,0,2))#runif(sampsize,-10,10)#seq(0.025,.975,length.out=20)
    x_samp2_org[reps_setup,] <- runif(sampsize,-10,10)+(0)#(x_samp_org[reps_setup,])^2
  y_samp_org[reps_setup,] <- multi*(x_samp_org[reps_setup,])+multi2*(x_samp2_org[reps_setup,])+
    (3+0/20+1*(x_samp_org[reps_setup,]-(-360))/20)^p_con*rnorm(sampsize,0,lambda) #rnorm(sampsize,0,lambda) #
  
  x_samp[reps_setup,] <- x_samp_org[reps_setup,] # log(x_samp_org[reps_setup,]) #  
  x_samp2[reps_setup,] <- x_samp2_org[reps_setup,] # log(x_samp2_org[reps_setup,]) #  
  y_samp[reps_setup,] <- y_samp_org[reps_setup,] # log(y_samp_org[reps_setup,]) #  
  
  med_est_x <- rq((x_samp[reps_setup,])~(1),tau=.5,weights=rep(1,sampsize))
  med_pt_est_x[reps_setup] <- med_est_x$coefficients[1]
  
    med_est_x2 <- rq((x_samp2[reps_setup,])~(1),tau=.5,weights=rep(1,sampsize))
    med_pt_est_x2[reps_setup,] <- med_est_x2$coefficients[1]
  
}

# setting quantile values of interest
taus_set <- c(0.025,.05,.1,.2,.25,.3,.4,.5,.6,.7,.75,.8,.9,.95,.975)
# taus_set <- c(0.5)
#taus_set <- c(.5,.6,.7,.75,.8,.9,.95,.975)

# calculating proxy population parameters for quantile regression based on whole of repeated samples

x_large <- c(x_samp)
x2_large <- c(x_samp2)
y_large <- c(y_samp)


bmk_model <- rq(y_large~x_large+x2_large,tau=taus_set,method="fn")
print(bmk_model) 
m_mk <- bmk_model$coefficient[seq(1,3*length(taus_set)-2,length.out=length(taus_set))]
l_mk <- bmk_model$coefficient[seq(2,3*length(taus_set)-1,length.out=length(taus_set))]
l2_mk <- bmk_model$coefficient[seq(3,3*length(taus_set),length.out=length(taus_set))]
print("estimated popn regression slopes and intercept")
print(taus_set)
print(round(m_mk,3))
print(round(l_mk,3))
print(round(l2_mk,3))



# graphing an example of the scatterplot to be analysed in quantile regression
for ( i in 1:1) {
x_sample <- x_samp[i,]
x_sample2 <- x_samp2[i,]
y_sample <- y_samp[i,]
n_obs <- length(x_sample)
t_adj <- abs(qt(.025,n_obs))

model_fit_rq <- rq((y_sample)~(x_sample+x_sample2),tau=.5)
model_fit_rq
rq_beta0 <- model_fit_rq$coefficients[1]
rq_beta1 <- model_fit_rq$coefficients[2]
rq_beta2 <- model_fit_rq$coefficients[3]
summary(model_fit_rq,se="boot",R=600)

s3d <- scatterplot3d(x=x_sample,y=x_sample2,z=y_sample,angle=-20)
# plot(y=y_sample,x=x_sample)#,ylim=c(-2,2),xlim=c(-2,2))#,log="xy")
s3d$points(z=(rq_beta0+rq_beta1*x_sample+rq_beta2*x_sample2),x=x_sample,y=x_sample2,pch=20,col="red")
s3d$points(z=fitted.values(rq((y_sample)~(x_sample+x_sample2),tau=.1)),x=x_sample,y=x_sample2,col="green",pch=20)
s3d$points(z=fitted.values(rq((y_sample)~(x_sample+x_sample2),tau=.9)),x=x_sample,y=x_sample2,col="blue",pch=20)
# points(y=fitted.values(rq((y_sample)~(x_sample+x_sample2),tau=.2)),x=x_sample,col="green",pch=20)
# points(y=fitted.values(rq((y_sample)~(x_sample+x_sample2),tau=.8)),x=x_sample,col="green",pch=20)
# points(y=fitted.values(rq((y_sample)~(x_sample+x_sample2),tau=.025)),x=x_sample,col="green",pch=20)
# points(y=fitted.values(rq((y_sample)~(x_sample+x_sample2),tau=.975)),x=x_sample,col="green",pch=20)
# points(y=fitted.values(rq(y_sample~1,tau=.5)),x=fitted.values(rq(x_sample~1,tau=.5)),col="magenta",pch=20)

# lmfit <- fitted.values(lm(y_sample~x_sample))
#   plot(y=y_sample-lmfit,x=x_sample,ylim=c(-6,6),main="lm - CI shown is one std")#,log="xy")
# abline(lm(y_sample-lmfit~x_sample),col="red")
# points(y=mean(y_sample-lmfit),x=mean(x_sample),col="magenta")
# lm_par <- summary(lm(y_sample-lmfit~x_sample))$coefficients[1:4]
# x_sort <- sort(x_sample);fit_sort <- lmfit[order(x_sample)]
# y_lm <- lm_par[1]+lm_par[2]*x_sort
# y_high_slope_int <- mean(y_sample-lmfit)-(lm_par[2]+lm_par[4])*mean(x_sample)
# y_high_slope_lm <- y_high_slope_int+(lm_par[2]+lm_par[4])*x_sort
# lines(y=y_high_slope_lm,x=x_sort,col="red",lty=2)
# y_high_int_lm <- lm_par[1]+lm_par[3]+0*lm_par[2]*x_sort
# lines(y=y_high_int_lm,x=x_sort,col="red",lty=2)
# y_high_comb_lm <- sqrt((y_high_int_lm-y_lm)^2+(y_high_slope_lm-y_lm)^2)+y_lm
# lines(y=y_high_comb_lm,x=x_sort,col="black",lty=2)
# y_low_slope_int <- mean(y_sample-lmfit)-(lm_par[2]-lm_par[4])*mean(x_sample)
# y_low_slope_lm <- y_low_slope_int+(lm_par[2]-lm_par[4])*x_sort
# lines(y=y_low_slope_lm,x=x_sort,col="red",lty=2)
# y_low_int_lm <- lm_par[1]-lm_par[3]+0*lm_par[2]*x_sort
# lines(y=y_low_int_lm,x=x_sort,col="red",lty=2)
# y_low_comb_lm <- -sqrt((y_low_int_lm-y_lm)^2+(y_low_slope_lm-y_lm)^2)+y_lm
# lines(y=y_low_comb_lm,x=x_sort,col="black",lty=2)
# 
# model_fit_rq <- rq(y_sample~x_sample,tau=0.5)
# rqfit <- fitted.values(model_fit_rq)
# plot(y=y_sample-rqfit,x=x_sample,ylim=c(-6,6),main="rq - CI shown is one std")#,log="xy")
# intrqfit <- rq(y_sample-rqfit~x_sample,tau=.5)
# intrqfithigh <- rq(y_sample-rqfit~x_sample,tau=.56)
# intrqfitlow <- rq(y_sample-rqfit~x_sample,tau=.44)
# 
# abline(h=(intrqfit$coefficients[1]+intrqfit$coefficients[2]*median(x_sample))/t_adj,col="green")
# abline(h=(intrqfithigh$coefficients[1]+intrqfithigh$coefficients[2]*median(x_sample))/t_adj,col="green",lty=2)
# abline(h=(intrqfitlow$coefficients[1]+intrqfitlow$coefficients[2]*median(x_sample))/t_adj,col="green",lty=2)
# points(y=rq(y_sample-rqfit~1,tau=.5)$coefficients,x=rq(x_sample~1,tau=.5)$coefficients,col="magenta")
# rq_par <- summary(rq(y_sample-rqfit~x_sample),se="boot")$coefficients[1:4]
# abline(h=median(x_sample)*rq_par[4],col="red",lty=3)
# abline(h=sqrt(rq_par[3]^2+(median(x_sample)*rq_par[4])^2),col="red",lty=3)
# abline(h=-median(x_sample)*rq_par[4],col="red",lty=3)
# abline(h=-sqrt(rq_par[3]^2+(median(x_sample)*rq_par[4])^2),col="red",lty=3)
# x_sort <- sort(x_sample)
# y_rq <- rq_par[1]+rq_par[2]*x_sort
# y_high_slope_int <- median(y_sample-rqfit)-(rq_par[2]+rq_par[4])*median(x_sample)
# y_high_slope_rq <- y_high_slope_int+(rq_par[2]+rq_par[4])*x_sort
# lines(y=y_high_slope_rq,x=x_sort,col="red",lty=3)
# y_high_int_rq <- rq_par[1]+rq_par[3]+0*rq_par[2]*x_sort
# lines(y=y_high_int_rq,x=x_sort,col="red",lty=3)
# y_high_comb_rq <- sqrt((y_high_int_rq-y_rq)^2+(y_high_slope_rq-y_rq)^2)+y_rq
# lines(y=y_high_comb_rq,x=x_sort,col="black",lty=3)
# y_low_slope_int <- median(y_sample-rqfit)-(rq_par[2]-rq_par[4])*median(x_sample)
# y_low_slope_rq <- y_low_slope_int+(rq_par[2]-rq_par[4])*x_sort
# lines(y=y_low_slope_rq,x=x_sort,col="red",lty=3)
# y_low_int_rq <- rq_par[1]-rq_par[3]+0*rq_par[2]*x_sort
# lines(y=y_low_int_rq,x=x_sort,col="red",lty=3)
# y_low_comb_rq <- -sqrt((y_low_int_rq-y_rq)^2+(y_low_slope_rq-y_rq)^2)+y_rq
# lines(y=y_low_comb_rq,x=x_sort,col="black",lty=3)
# yres <- resid(model_fit_rq)
# y2 <- abs(yres)
# res_model <- rq(y2~poly(x_sample,6),tau=0.5)
# var_est <- fitted.values(res_model)
# q_maxes <- rq(yres~1,tau=c(.025,.975))$coefficients[1:2]
# 
# trimdatax <- subset(x_sample,!(yres > q_maxes[2] | yres < q_maxes[1] | var_est == 0 )) 
# trimdatay <- subset(y2,!(yres > q_maxes[2] | yres < q_maxes[1] | var_est == 0 )) 
# trimdatavar_est <- subset(abs(var_est),!(yres > q_maxes[2] | yres < q_maxes[1] | var_est == 0  )) 
# var_adj <- var_est/(suppressMessages(sintegral(trimdatax,trimdatay)$value)/
#                       suppressMessages(sintegral(trimdatax,trimdatay/trimdatavar_est)$value))
# pts_adj <- ifelse(!(var_est ==0 | abs(yres/var_adj) > max(abs(yres))  ), 1, 0)
# #   print(paste(reps,sum(pts_adj)))
# res_adj <- ifelse(pts_adj,yres/var_adj,yres)
# wgt_adj <- ifelse(pts_adj,1/var_adj,1)
# wgt_adj <- ifelse(wgt_adj > 0,wgt_adj,0)
# 
# res_adjfin <- res_adj
# 
# #     plot(y=y2,x=x_sample)
# #      points(y=var_est,x=x_sample,col="red")
# #     plot(y=yres,x=x_sample,ylim=c(-max(y2,abs(res_adj)),max(y2,abs(res_adj))))
# #     # points(y=res_adj,x=x_sample,col="red")
# #     points(y=res_adjfin,x=x_sample,col="green")
# #      plot(y=res_adjfin,x=x_sample)
# #      points(y=fitted.values(rq((res_adjfin)~(x_sample),tau=.5)),x=x_sample,col="red",pch=20)
# #      points(y=fitted.values(rq((res_adjfin)~(x_sample),tau=.1)),x=x_sample,col="green",pch=20)
# #      points(y=fitted.values(rq((res_adjfin)~(x_sample),tau=.9)),x=x_sample,col="green",pch=20)
# #      points(y=fitted.values(rq((res_adjfin)~(x_sample),tau=.25)),x=x_sample,col="green",pch=20)
# #      points(y=fitted.values(rq((res_adjfin)~(x_sample),tau=.75)),x=x_sample,col="green",pch=20)
# #      points(y=fitted.values(rq((res_adjfin)~(x_sample),tau=.025)),x=x_sample,col="green",pch=20)
# #      points(y=fitted.values(rq((res_adjfin)~(x_sample),tau=.975)),x=x_sample,col="green",pch=20)
# 
# 
# samp <- res_adjfin  
# slrqfit <- rq(samp~1,tau=.5)
# slrqfithigh <- rq(samp~1,tau=.56)
# slrqfitlow <- rq(samp~1,tau=.44)
# s2x1 <- sum(x_sample^2)-sum(x_sample)^2/n_obs
# r12 <- 0  
# se1adj <- sqrt((s2x1)*(1-r12)/n_obs)
# sl_h <- slrqfithigh$coefficients[1]/se1adj/t_adj*median(x_sample)
# sl_l <- slrqfitlow$coefficients[1]/se1adj/t_adj*median(x_sample)
# abline(h=sl_h,col="green",lty=2)
# abline(h=sl_l,col="green",lty=2)
# abline(h=sqrt(y_high_int_rq^2+sl_h^2),col="green",lty=2)
# abline(h=-sqrt(y_low_int_rq^2+sl_l^2),col="green",lty=2)




}

# setting weight values of sample
wgt_filereg <- rep(1,length(x_sample))#y_sample#
wgt_fileregexp <- rep(1,400)

# END MANUAL SETTING 1


bcov <- 0
bincov <- 0
e_symmcov <- 0

bcov2 <- 0
bincov2 <- 0
e_symmcov2 <- 0

by0_cov <- 0
e_symmy0_cov <- 0
e_symmy0_cov_u <- 0
e_symmy0_cov_w <- 0
biny0_cov <- 0
biny0_cov_u <- 0
biny0_cov_w <- 0



  
    
for (j in 1:length(taus_set))  {


  
  taus <- taus_set[j]
  
  boot_cov <- 0
  evd_symm_cov <- 0
  bin_cov <- 0
  
  boot_cov2 <- 0
  evd_symm_cov2 <- 0
  bin_cov2 <- 0
  
  boot_y0_cov <- 0
  boot_y0_cov_u <- 0
  boot_y0_cov_w <- 0
  evd_symm_y0_cov <- 0
  evd_symm_y0_cov_u <- 0
  evd_symm_y0_cov_w <- 0
  bin_y0_cov <- 0
  bin_y0_cov_u <- 0
  bin_y0_cov_w <- 0
  
  
  boot_low <- 0
  boot_high <- 0
  evd_low <- 0
  evd_high <- 0
  evd_symm_low <- 0
  evd_symm_high <- 0
  bin_low <- 0
  bin_high <- 0
  
  boot_low2 <- 0
  boot_high2 <- 0
  evd_low2 <- 0
  evd_high2 <- 0
  evd_symm_low2 <- 0
  evd_symm_high2 <- 0
  bin_low2 <- 0
  bin_high2 <- 0
  
  boot_y0_low <- 0
  boot_y0_high <- 0
  evd_y0_low <- 0
  evd_y0_high <- 0
  evd_symm_y0_low <- 0
  evd_symm_y0_high <- 0
  evd_symm_y0_low_u <- 0
  evd_symm_y0_high_u <- 0
  evd_symm_y0_low_w <- 0
  evd_symm_y0_high_w <- 0
  bin_y0_low <- 0
  bin_y0_high <- 0
  bin_y0_low_u <- 0
  bin_y0_high_u <- 0
  bin_y0_low_w <- 0
  bin_y0_high_w <- 0
  
  samp_var <- 0
  y_var <- 0
  samp_var_theta <- 0  
  samp_r1 <- 0
  samp_low <- 0
  samp_high <- 0
  pt_est <- 0
  pt_est_b2 <- 0
  pt_est_b1 <- 0
  pt_est_b0 <- 0
  
  irep <- 1

  x_sample <- x_samp[1,]
  x_sample2 <- x_samp2[1,]
  y_sample <- y_samp[1,]
  model_fit_rq <- rq((y_sample)~(x_sample+x_sample2),tau=taus,weights=wgt_filereg)
  model_fit_rq
  rq_beta0 <- model_fit_rq$coefficients[1]
  rq_beta1 <- model_fit_rq$coefficients[2]
  rq_beta2 <- model_fit_rq$coefficients[3]
  summary(model_fit_rq,se="boot",R=600)
  samp <- y_sample-(rq_beta0+rq_beta1*x_sample+rq_beta2*x_sample2)
  samp_int <- samp
  
  n_samp <- length(samp)
  
  # quantreg estimates of X distribution values for given $\theta$ values
  
  
if(method_df == "perc_scale") {
  q_var <- sqrt(taus*(1-taus)/(n_samp-3))} else
    {q_var <- sqrt(taus*(1-taus)/(n_samp))}
  
  p_low <- qnorm(0.025,taus,q_var)#*sqrt(2))
  p_high <- qnorm(0.975,taus,q_var)#*sqrt(2))

  if(method_df == "perc_scale")  {b_low <- qbinom(.025,(sampsize-3)/1,taus)/((sampsize-3)/1)
  b_high <- qbinom(.975,(sampsize-3)/1,taus)/((sampsize-3)/1)}  else  {b_low <- qbinom(.025,sampsize,taus)/(sampsize)
  b_high <- qbinom(.975,sampsize,taus)/(sampsize)}
  
    p_low <- max(p_low,0)
    p_high <- min(p_high,1)
  
  
  print("binomial CI | evd CI  , in percentile scale")
  print("lower bounds")
  print(paste(b_low,p_low))
  print("upper bounds")
  print(paste(b_high,p_high))
  
  
  for (reps in 1:replicates) {
    
    x_sample <- x_samp[reps,]
    x_sample2 <- x_samp2[reps,]
    y_sample <- y_samp[reps,]
    model_fit_rq <- rq((y_sample)~(x_sample+x_sample2),tau=taus,weights=wgt_filereg)
    model_fit_rq
    rqr_beta0 <- model_fit_rq$coefficients[1]
    rqr_beta1 <- model_fit_rq$coefficients[2]
    rqr_beta2 <- model_fit_rq$coefficients[3]
    std_errors <- summary(model_fit_rq,se="boot",R=600)$coefficients[c(4:6)]
    beta2_std <- std_errors[3]
    beta1_std <- std_errors[2]
    beta0_std <- std_errors[1]
    
if(adj_het == "adjHet")   {
  
 # trying to adjust for heterogeniety
    yres <- resid(model_fit_rq)
    y2 <- abs(yres)
    x2 <- x_sample^2
    res_model <- rq(y2~poly(x_sample,6),tau=0.5)
   # res_model <- rq(y2~x_sample+x_sample^3+x_sample^4+x_sample2+x_sample2^2+x_sample2^3,tau=0.5)
    var_est <- fitted.values(res_model)
   q_maxes <- rq(yres~1,tau=c(.025,.975))$coefficients[1:2]
   
   trimdatax <- subset(x_sample,!(yres > q_maxes[2] | yres < q_maxes[1] | var_est == 0 )) 
   trimdatay <- subset(y2,!(yres > q_maxes[2] | yres < q_maxes[1] | var_est == 0 )) 
    trimdatavar_est <- subset(abs(var_est),!(yres > q_maxes[2] | yres < q_maxes[1] | var_est == 0  )) 
    var_adj <- var_est/(suppressMessages(sintegral(trimdatax,trimdatay)$value)/
                          suppressMessages(sintegral(trimdatax,trimdatay/trimdatavar_est)$value))
    pts_adj <- ifelse(!(var_est ==0 | abs(yres/var_adj) > max(abs(yres))  ), 1, 0)
  #   print(paste(reps,sum(pts_adj)))

          res_adj <- ifelse(pts_adj,yres/var_adj,yres)
      wgt_adj <- ifelse(pts_adj,1/var_adj,1)
      wgt_adj <- ifelse(wgt_adj > 0,wgt_adj,0)
      
     res_adjfin <- res_adj#*max(y2)/max(abs(res_adj))
     
#     plot(y=y2,x=x_sample)
#      points(y=var_est,x=x_sample,col="red")
# 
#    plot(y=yres,x=x_sample,ylim=c(-max(y2,abs(res_adj)),max(y2,abs(res_adj))))
#     points(y=res_adj,x=x_sample,col="red")
#     points(y=res_adjfin,x=x_sample,col="green")

    #      plot(y=res_adjfin,x=x_sample)
#      points(y=fitted.values(rq((res_adjfin)~(x_sample),tau=.5)),x=x_sample,col="red",pch=20)
#      points(y=fitted.values(rq((res_adjfin)~(x_sample),tau=.1)),x=x_sample,col="green",pch=20)
#      points(y=fitted.values(rq((res_adjfin)~(x_sample),tau=.9)),x=x_sample,col="green",pch=20)
#      points(y=fitted.values(rq((res_adjfin)~(x_sample),tau=.25)),x=x_sample,col="green",pch=20)
#      points(y=fitted.values(rq((res_adjfin)~(x_sample),tau=.75)),x=x_sample,col="green",pch=20)
#      points(y=fitted.values(rq((res_adjfin)~(x_sample),tau=.025)),x=x_sample,col="green",pch=20)
#      points(y=fitted.values(rq((res_adjfin)~(x_sample),tau=.975)),x=x_sample,col="green",pch=20)
     med_est_x <- rq((x_samp[reps,])~(1),tau=.5,weights=wgt_adj^1)
     med_pt_est_x_wgt[reps] <- med_est_x$coefficients[1]
     med_est_x2 <- rq((x_samp2[reps,])~(1),tau=.5,weights=wgt_adj^1)
     med_pt_est_x2_wgt[reps] <- med_est_x2$coefficients[1]
     
    
    # var_est <- ifelse(var_est > 0, var_est, 0.0001)
#    model_fit_rq2 <- rq(y_sample~x_sample,tau=taus,weights=1/var_est)
#    samp <- resid(model_fit_rq2)
    samp <- res_adjfin  
}
    else {  samp <- y_sample-(rqr_beta0+rqr_beta1*x_sample+rqr_beta2*x_sample2) }
    
    samp_var[reps] <- var(samp)
    y_var[reps] <- var(y_sample)
    
#     #redo modelling as weighted regression
#     model_fit_rq <- rq((y_sample)~(x_sample),tau=taus,weights=1/var_est)
#     model_fit_rq
#     rqr_beta0 <- model_fit_rq$coefficients[1]
#     rqr_beta1 <- model_fit_rq$coefficients[2]
#     # rqr_beta2 <- model_fit_rq$coefficients[3]
#     std_errors <- summary(model_fit_rq,se="boot",R=600)$coefficients[c(3:4)]
#     # beta2_std <- std_errors[3]
#     beta1_std <- std_errors[2]
#     beta0_std <- std_errors[1]
    
    
    rho <- function(u,taur=.5) u*(taur - (u < 0))
    samp_r1[reps] <- 1 - model_fit_rq$rho/rq(y_sample~1,tau=taus,weights=wgt_filereg)$rho
    
    
    n_row <- length(taus)
    
    qr_dat <- matrix(c(1:1*irep),nrow=irep,ncol=1,
                     dimnames=list(taus,c("beta0")))
    qp_lstd <- matrix(c(1:1*irep),nrow=irep,ncol=1,
                      dimnames=list(taus,c("lstd")))
    qp_ustd <- matrix(c(1:1*irep),nrow=irep,ncol=1,
                      dimnames=list(taus,c("ustd")))
    backtrans_std <- matrix(c(1:3*irep),nrow=irep,ncol=3,
                            dimnames=list(taus,c("lower_std_dist","upper_std_dist","average_std_value")))
    qrboot_dat <- matrix(c(1:4*irep),nrow=irep,ncol=4,
                         dimnames=list(taus,c("beta0","std_error","t_value","Pr")))
    
    
    
    for (i in 1:length(taus)) {
      

      model_fit <- rq((samp)~(1),tau=taus[i],method="fn")
      #      print(summary(model_fit,se="boot",R=600))
      qr_dat[i,] <- c(summary(model_fit)$coefficients[1])
      qrboot_dat[i,] <- c(summary(model_fit,se="boot",R=600)$coefficients[1:4])
      
      qp_lstd[i] <- taus[i]-q_var[i]
      qp_ustd[i] <- taus[i]+q_var[i]
      taustd <- c(qp_lstd[i],qp_ustd[i])
      
      #       plot(y=samp,x=y_sample)
      #       mod2 <- rq((samp)~(y_sample),tau=taus[i],method="fn",weights=wgt_filereg)
      #       
      model_fit_lstd <- summary(rq((samp)~(1),tau=max(qp_lstd[i],0.00005),method="fn"))$coefficients[1]
      model_fit_ustd <- summary(rq((samp)~(1),tau=min(qp_ustd[i],.99995),method="fn"))$coefficients[1]
      
      backtrans_std[i,] <- c(model_fit_lstd-qr_dat[i],model_fit_ustd-qr_dat[i],(model_fit_ustd-model_fit_lstd)/2)
      
      
      
    }
    

    back_pts1 <- rq((samp)~(1),tau=p_low)$coefficients[1]
    
    back_pts2 <- rq((samp)~(1),tau=p_high)$coefficients[1]

    back_pts1a <- sum(rq((y_sample-(rqr_beta0+rqr_beta1*x_sample+rqr_beta2*x_sample2))~(x_sample+x_sample2),
                         tau=p_low)$coefficients[1:3]*(c(1,1*med_pt_est_x_wgt[reps],1*med_pt_est_x2_wgt[reps])))
    
    back_pts3a <- sum(rq((y_sample-(rqr_beta0+rqr_beta1*x_sample+rqr_beta2*x_sample2))~(x_sample+x_sample2),
                         tau=taus[i])$coefficients[1:3]*(c(1,1*med_pt_est_x_wgt[reps],1*med_pt_est_x2_wgt[reps])))

    back_pts2a <- sum(rq((y_sample-(rqr_beta0+rqr_beta1*x_sample+rqr_beta2*x_sample2))~(x_sample+x_sample2),
                         tau=p_high)$coefficients[1:3]*(c(1,1*med_pt_est_x_wgt[reps],1*med_pt_est_x2_wgt[reps])))
    
    back_pts3 <- summary(rq((samp)~(1),tau=b_low))$coefficients[1]
    
    back_pts4 <- summary(rq((samp)~(1),tau=b_high))$coefficients[1]

    back_pts5a <- sum(rq((y_sample-(rqr_beta0+rqr_beta1*x_sample+rqr_beta2*x_sample2))~(x_sample+x_sample2),
                         tau=b_low)$coefficients[1:3]*(c(1,1*med_pt_est_x_wgt[reps],1*med_pt_est_x2_wgt[reps])))
    
    back_pts6a <- sum(rq((y_sample-(rqr_beta0+rqr_beta1*x_sample+rqr_beta2*x_sample2))~(x_sample+x_sample2),
                         tau=b_high)$coefficients[1:3]*(c(1,1*med_pt_est_x_wgt[reps],1*med_pt_est_x2_wgt[reps])))
    
    
    
    
    bk_pt <- c(back_pts1,back_pts2,back_pts3,
               back_pts4,back_pts1a,back_pts2a,back_pts3a,
               back_pts5a,back_pts6a)
    

   s2x1x2 <- sum(x_sample*x_sample2)-sum(x_sample)*sum(x_sample2)/n_samp
    s2x1 <- sum(x_sample^2)-sum(x_sample)^2/n_samp
    s2x2 <- sum(x_sample2^2)-sum(x_sample2)^2/n_samp
       
      r12 <- s2x1x2/sqrt(s2x1*s2x2)
    # r12 <- 0  
     
    #     popvar1 <- sqrt((s2x1)/(n_samp))
    #     popvar2 <- sqrt((s2x2)/(n_samp))
    #     dfadj <- (n_samp/(n_samp-3))
    
    if(method_df == "perc_scale")   
    {se1adj <- sqrt((s2x1)*(1-r12)/n_samp)
    se2adj <- sqrt((s2x2)*(1-r12)/n_samp)}
      else
      {se1adj <- sqrt((s2x1)*(1-r12)/(1/(n_samp-3)))/n_samp
      se2adj <- sqrt((s2x2)*(1-r12)/(1/(n_samp-3)))/n_samp}
    
        # var adj for beta1/2 plus df adj for res var
    ### df correction in original frame
    #    se1adj <- sqrt((s2x1)*(1-r12)/(1/(n_samp-2)))/n_samp
    
    ### df correction in percentile frame
    #se1adj <- sqrt((s2x1)*(1-r12)/n_samp)
      
    
    boot_low[reps] <- rqr_beta1-t_adj*beta1_std# 3 for 2 parameter model, 2 for one pa
    boot_high[reps] <- rqr_beta1+t_adj*beta1_std# 3 for 2 parameter model, 2 for one pa
    evd_low[reps] <- rqr_beta1+bk_pt[1]/se1adj
    evd_high[reps] <- rqr_beta1+bk_pt[2]/se1adj
        #evd_ci <- abs((bk_pt[1]-bk_pt[2])/2)
    evd_ci <- max(abs(bk_pt[1]),abs(bk_pt[2]))
    # print(paste("beta1",(t_adj*beta1_std)*se1adj,evd_ci,sqrt(2)*evd_ci))
    # print(paste(evd_ci,se1adj,r12))
    evd_mean_ci <- mean(abs(bk_pt[1]),abs(bk_pt[2]))
    evd_symm_low[reps] <- rqr_beta1-evd_ci/se1adj
    evd_symm_high[reps] <- rqr_beta1+evd_ci/se1adj
        #bin_ci <- abs((bk_pt[3]-bk_pt[4])/2)
    bin_ci <- max(abs(bk_pt[3]),abs(bk_pt[4]))
    bin_mean_ci <- mean(abs(bk_pt[3]),abs(bk_pt[4]))
    bin_low[reps] <- rqr_beta1-bin_ci/se1adj
    bin_high[reps] <- rqr_beta1+bin_ci/se1adj
    evd_hom_ci <- max(abs(bk_pt[7]-bk_pt[5]),abs(bk_pt[7]-bk_pt[6]))
    bin_hom_ci <- max(abs(bk_pt[7]-bk_pt[8]),abs(bk_pt[7]-bk_pt[9]))
    #     bin_ci <- max(abs(bk_pt[5]),abs(bk_pt[6]))
#     bin_low[reps] <- rqr_beta1-bin_ci/se1adj
#     bin_high[reps] <- rqr_beta1+bin_ci/se1adj
    
    ifelse(l_mk[j] >= boot_low[reps] & l_mk[j] <= boot_high[reps],
           boot_cov[reps] <-1, boot_cov[reps] <- 0)
    
    ifelse(l_mk[j] >= evd_symm_low[reps] & l_mk[j] <= evd_symm_high[reps],
           evd_symm_cov[reps] <- 1,evd_symm_cov[reps] <- 0)
    
    ifelse(l_mk[j] >= bin_low[reps] & l_mk[j] <= bin_high[reps],
           bin_cov[reps] <- 1,bin_cov[reps] <- 0)
 
    # print(paste("beta2",(t_adj*beta2_std)*se2adj,evd_ci))
    
    boot_low2[reps] <- rqr_beta2-t_adj*beta2_std# 3 for 2 parameter model, 2 for one pa
    boot_high2[reps] <- rqr_beta2+t_adj*beta2_std# 3 for 2 parameter model, 2 for one pa
    #    evd_ci <- abs((bk_pt[1]-bk_pt[2])/2)
    evd_low2[reps] <- rqr_beta2+bk_pt[1]/se2adj
    evd_high2[reps] <- rqr_beta2+bk_pt[2]/se2adj
    #evd_ci <- max(abs(bk_pt[1]),abs(bk_pt[2]))
    evd_symm_low2[reps] <- rqr_beta2-evd_ci/se2adj
    evd_symm_high2[reps] <- rqr_beta2+evd_ci/se2adj
    #    bin_ci <- abs((bk_pt[3]-bk_pt[4])/2)
    #bin_ci <- max(abs(bk_pt[3]),abs(bk_pt[4]))
    bin_low2[reps] <- rqr_beta2-bin_ci/se2adj
    bin_high2[reps] <- rqr_beta2+bin_ci/se2adj
    
    ifelse(l2_mk[j] >= boot_low2[reps] & l2_mk[j] <= boot_high2[reps],
           boot_cov2[reps] <-1, boot_cov2[reps] <- 0)
    
    ifelse(l2_mk[j] >= evd_symm_low2[reps] & l2_mk[j] <= evd_symm_high2[reps],
           evd_symm_cov2[reps] <- 1,evd_symm_cov2[reps] <- 0)
    
    ifelse(l2_mk[j] >= bin_low2[reps] & l2_mk[j] <= bin_high2[reps],
           bin_cov2[reps] <- 1,bin_cov2[reps] <- 0)
    
    samp_low[reps] <- min(samp)
    samp_high[reps] <- max(samp)
    
    pt_est[reps] <- qr_dat
    pt_est_b1[reps] <- rqr_beta1 
    pt_est_b0[reps] <- rqr_beta0 
    pt_est_b2[reps] <- rqr_beta2
    
    #    denxvar2 <- sqrt(n_samp/(n_samp-2)+(pt_est_x[reps])^2/(se1adj)^2)
    #denxvar2 <- sqrt(1+(pt_est_x[reps])^2/(se1adj)^2)
    
    
    
    if(method_df == "perc_scale")   
    {denxevd0 <- sqrt(1)
    denxbin0 <- sqrt(1)
#     denxevdu <- sqrt(1+(med_pt_est_x[reps]*evd_ci/evd_hom_ci)^2/(se1adj)^2)
#     denxbinu <- sqrt(1+(med_pt_est_x[reps]*bin_ci/bin_hom_ci)^2/(se1adj)^2)
#     denxevdw <- sqrt(0+(med_pt_est_x[reps]*evd_ci/evd_hom_ci)^2/(se1adj)^2)
#     denxbinw <- sqrt(0+(med_pt_est_x_wgt[reps]*bin_ci/bin_hom_ci)^2/(se1adj)^2)
    denxevdu <- sqrt(evd_hom_ci^2+(med_pt_est_x[reps]*evd_ci/se1adj)^2+(med_pt_est_x2[reps]*evd_ci/se2adj)^2)
    denxbinu <- sqrt(bin_hom_ci^2+(med_pt_est_x[reps]*bin_ci/se1adj)^2+(med_pt_est_x2[reps]*bin_ci/se2adj)^2)
    denxevdw <- sqrt(evd_hom_ci^2+(med_pt_est_x_wgt[reps]*evd_mean_ci/se1adj)^2+
                       (med_pt_est_x2_wgt[reps]*evd_mean_ci/se2adj)^2)
    denxbinw <- sqrt(bin_hom_ci^2+(med_pt_est_x_wgt[reps]*bin_mean_ci/se1adj)^2+
                       (med_pt_est_x2_wgt[reps]*bin_mean_ci/se2adj)^2)
    }
#     +
#         (med_pt_est_x2[reps])^2/(se2adj)^2)}
      else
      {denxvar2 <- sqrt(n_samp/(n_samp-2)+(med_pt_est_x[reps])^2/(se1adj)^2)
      denxbin2 <- sqrt(n_samp/(n_samp-2)+(med_pt_est_x[reps])^2/(se1adj)^2)}
#     +
#         (med_pt_est_x2[reps])^2/(se2adj)^2)}
    

        ### df correction in original frame
    #        denxvar2 <- sqrt(n_samp/(n_samp-2))
    
    ### df correction in percentile frame
    #denxvar2 <- 1
    
    
    #    denxboot2 <- 1
    #    denxboot2 <- sqrt(1+(pt_est_x[reps])^2/((uncond_evd_high_x2std[reps]-
    #                                               uncond_evd_low_x2std[reps])/2/t_adj)^2/(n_samp-1)/(n_samp-1))
#     denxevd2 <- denxvar2
#     
#     denxbin2 <- denxvar3
    
    boot_y0_low[reps] <- rqr_beta0-t_adj*beta0_std
    boot_y0_high[reps] <- rqr_beta0+t_adj*beta0_std
    evd_y0_low[reps] <- rqr_beta0+bk_pt[1]*denxevdu
    evd_y0_high[reps] <- rqr_beta0+bk_pt[2]*denxevdu
    evd_symm_y0_low[reps] <- rqr_beta0-evd_hom_ci*denxevd0
    evd_symm_y0_high[reps] <- rqr_beta0+evd_hom_ci*denxevd0
    bin_y0_low[reps] <- rqr_beta0-bin_hom_ci*denxbin0
    bin_y0_high[reps] <- rqr_beta0+bin_hom_ci*denxbin0
    evd_symm_y0_low_u[reps] <- rqr_beta0-denxevdu
    evd_symm_y0_high_u[reps] <- rqr_beta0+denxevdu
#     evd_symm_y0_low_w[reps] <- rqr_beta0-evd_hom_ci*denxevdw
#     evd_symm_y0_high_w[reps] <- rqr_beta0+evd_hom_ci*denxevdw
        evd_symm_y0_low_w[reps] <- rqr_beta0-denxevdw
        evd_symm_y0_high_w[reps] <- rqr_beta0+denxevdw
    bin_y0_low_u[reps] <- rqr_beta0-denxbinu
    bin_y0_high_u[reps] <- rqr_beta0+denxbinu
    bin_y0_low_w[reps] <- rqr_beta0-denxbinw
    bin_y0_high_w[reps] <- rqr_beta0+denxbinw
    
    if(method_df == "perc_scale") {  samp_var_theta[reps] <- (bin_ci/t_adj)^2}
     else 
      {samp_var_theta[reps] <- (bin_ci/t_adj)^2*n_samp/(n_samp-3)}


        ### df correction in original frame
#    samp_var_theta[reps] <- (bin_ci/t_adj)^2*n_samp/(n_samp-2)
    
    ### df correction in percentile frame
#    samp_var_theta[reps] <- (bin_ci/t_adj)^2
    
    
    ifelse(m_mk[j] >= boot_y0_low[reps] & m_mk[j] <= boot_y0_high[reps],
           boot_y0_cov[reps] <-1, boot_y0_cov[reps] <- 0)
    
    ifelse(m_mk[j] >= evd_symm_y0_low[reps] & m_mk[j] <= evd_symm_y0_high[reps],
           evd_symm_y0_cov[reps] <- 1,evd_symm_y0_cov[reps] <- 0)
    
    ifelse(m_mk[j] >= bin_y0_low[reps] & m_mk[j] <= bin_y0_high[reps],
           bin_y0_cov[reps] <- 1,bin_y0_cov[reps] <- 0)

    ifelse(m_mk[j] >= evd_symm_y0_low_u[reps] & m_mk[j] <= evd_symm_y0_high_u[reps],
           evd_symm_y0_cov_u[reps] <- 1,evd_symm_y0_cov_u[reps] <- 0)
    
    ifelse(m_mk[j] >= evd_symm_y0_low_w[reps] & m_mk[j] <= evd_symm_y0_high_w[reps],
           evd_symm_y0_cov_w[reps] <- 1,evd_symm_y0_cov_w[reps] <- 0)
    
    ifelse(m_mk[j] >= bin_y0_low_u[reps] & m_mk[j] <= bin_y0_high_u[reps],
           bin_y0_cov_u[reps] <- 1,bin_y0_cov_u[reps] <- 0)
    
    ifelse(m_mk[j] >= bin_y0_low_w[reps] & m_mk[j] <= bin_y0_high_w[reps],
           bin_y0_cov_w[reps] <- 1,bin_y0_cov_w[reps] <- 0)
    
    
    
  }
  
  bcov[j] <- mean(boot_cov)
  n_reps <- length(boot_cov)
  bincov[j] <- mean(bin_cov)
  e_symmcov[j] <- mean(evd_symm_cov)
  
  print("")
  print("coverage performance beta1 and mean quantile regression estimate")
  print("q pt | bootstrap | binomial | evdsymm | quantreg est | pop value")
  print(paste(taus,bcov[j],bincov[j],e_symmcov[j],mean(pt_est_b1),l_mk[j]))
  
  bcov2[j] <- mean(boot_cov2)
  #n_reps <- length(boot_cov)
  bincov2[j] <- mean(bin_cov2)
  e_symmcov2[j] <- mean(evd_symm_cov2)
  
  print("")
  print("coverage performance beta2 and mean quantile regression estimate")
  print("q pt | bootstrap | binomial | evdsymm | quantreg est | pop value")
  print(paste(taus,bcov2[j],bincov2[j],e_symmcov2[j],mean(pt_est_b2),l2_mk[j]))
  
  by0_cov[j] <- mean(boot_y0_cov)
  n_reps <- length(boot_y0_cov)
  biny0_cov[j] <- mean(bin_y0_cov)
  biny0_cov_u[j] <- mean(bin_y0_cov_u)
  biny0_cov_w[j] <- mean(bin_y0_cov_w)
  e_symmy0_cov[j] <- mean(evd_symm_y0_cov)
  e_symmy0_cov_u[j] <- mean(evd_symm_y0_cov_u)
  e_symmy0_cov_w[j] <- mean(evd_symm_y0_cov_w)
  
  print("")
  print("coverage performance beta0 and mean quantile regression estimate")
  print("q pt | bootstrap | binomial | bin_u | bin_w | evdsymm | evdsymm_u |  evdsymm_w | quantreg est | pop value")
  print(paste(taus,by0_cov[j],biny0_cov[j],biny0_cov_u[j],biny0_cov_w[j],e_symmy0_cov[j],e_symmy0_cov_u[j],e_symmy0_cov_w[j],mean(pt_est_b0),m_mk[j]))
  
  
  boot_inc <- 0
  evd_symm_inc <- 0
  bin_inc <- 0
 
  boot_inc2 <- 0
  evd_symm_inc2 <- 0
  bin_inc2 <- 0
  
  
  boot_y0_inc <- 0
  evd_symm_y0_inc <- 0
  bin_y0_inc <- 0
  
  for (i in 1:n_reps) {
    boot_inc[i] <- mean(boot_cov[1:i])
    evd_symm_inc[i] <- mean(evd_symm_cov[1:i])
    bin_inc[i] <- mean(bin_cov[1:i])

        boot_inc2[i] <- mean(boot_cov2[1:i])
    evd_symm_inc2[i] <- mean(evd_symm_cov2[1:i])
    bin_inc2[i] <- mean(bin_cov2[1:i])
    
    
    boot_y0_inc[i] <- mean(boot_y0_cov[1:i])
    evd_symm_y0_inc[i] <- mean(evd_symm_y0_cov[1:i])
    bin_y0_inc[i] <- mean(bin_y0_cov[1:i])
    
  }
  
  
# }

boot_low_r2std <- boot_low
boot_high_r2std <- boot_high
bin_low_r2std <- bin_low
bin_high_r2std <- bin_high
evd_symm_low_r2std <- evd_symm_low
evd_symm_high_r2std <- evd_symm_high
evd_low_r2std <- evd_low
evd_high_r2std <- evd_high

mean(boot_cov)
mean(evd_symm_cov)
mean(bin_cov)

            


# END MANUAL SETTING 6

# hist(rqr_beta0,breaks=50)
# hist(rqr_beta1,breaks=50)


plot(boot_inc,ylim=c(.8,1))
points(bin_inc,col="green")
points(evd_symm_inc,col="brown")

mean(pt_est_b1-boot_low)/t_adj
mean(pt_est_b1-evd_symm_low)/t_adj
mean(pt_est_b1-bin_low)/t_adj
mean(pt_est_b1-evd_low)/t_adj
mean(pt_est_b1-evd_high)/t_adj

sum(ifelse(boot_high < l_mk[j],1,0))/length(boot_high)
sum(ifelse(boot_low > l_mk[j],1,0))/length(boot_high)
sum(ifelse(evd_symm_high < l_mk[j],1,0))/length(boot_high)
sum(ifelse(evd_symm_low > l_mk[j],1,0))/length(boot_high)
sum(ifelse(bin_high < l_mk[j],1,0))/length(boot_high)
sum(ifelse(bin_low > l_mk[j],1,0))/length(boot_high)

hist(pt_est_b1,breaks=100)
# hist(boot_high,breaks=100)
# hist(boot_low,breaks=100)
# hist(evd_symm_high,breaks=100)
# hist(evd_symm_low,breaks=100)
# hist(bin_high,breaks=100)
# hist(bin_low,breaks=100)



mean(boot_y0_cov)
mean(evd_symm_y0_cov)
mean(bin_y0_cov)

           

# hist(rqr_beta0,breaks=50)
# hist(rqr_beta1,breaks=50)


plot(boot_y0_inc,ylim=c(.8,1))
points(evd_symm_y0_inc,col="brown")
points(bin_y0_inc,col="green")

mean(pt_est_b0-boot_y0_low)/t_adj
mean(pt_est_b0-evd_symm_y0_low)/t_adj
mean(pt_est_b0-bin_y0_low)/t_adj
mean(pt_est_b0-evd_y0_low)/t_adj
mean(pt_est_b0-evd_y0_high)/t_adj
            
sum(ifelse(boot_y0_high < m_mk[j],1,0))/length(boot_high)
sum(ifelse(boot_y0_low > m_mk[j],1,0))/length(boot_high)
sum(ifelse(evd_symm_y0_high < m_mk[j],1,0))/length(boot_high)
sum(ifelse(evd_symm_y0_low > m_mk[j],1,0))/length(boot_high)
sum(ifelse(bin_y0_high < m_mk[j],1,0))/length(boot_high)
sum(ifelse(bin_y0_low > m_mk[j],1,0))/length(boot_high)

hist(pt_est_b0,breaks=100)
# hist(boot_y0_high,breaks=100)
# hist(boot_y0_low,breaks=100)
# hist(evd_symm_y0_high,breaks=100)
# hist(evd_symm_y0_low,breaks=100)
# hist(bin_y0_high,breaks=100)
# hist(bin_y0_low,breaks=100)

mea_lm <- 0
for (i in 1:replicates) mea_lm[i] <- summary(lm(y_samp[i,]~x_samp[i,]+x_samp2[i,]))$r.squared
mean(mea_lm)

rsq <- 1- samp_var/y_var;print(mean(rsq))
#rsq_theta <- 1- samp_var_theta/y_var_theta;print(mean(rsq_theta))
print(mean(samp_r1))
# hist(rsq,breaks=30);#hist(rsq_theta,breaks=30);
# hist(samp_r1,breaks=30)
# hist(pt_est,breaks=30)

print(j)
yyers  <- c(taus,l_mk[j],mean(pt_est_b1),
            mean(boot_cov),mean(evd_symm_cov),mean(bin_cov),
            mean(pt_est_b1-boot_low)/t_adj,
            mean(pt_est_b1-evd_symm_low)/t_adj,
            mean(pt_est_b1-bin_low)/t_adj,
            mean(pt_est_b1-evd_low)/t_adj,
            mean(pt_est_b1-evd_high)/t_adj,
            m_mk[j],mean(pt_est_b0),
            mean(boot_y0_cov),
            mean(evd_symm_y0_cov_w),mean(bin_y0_cov_w),
            mean(evd_symm_y0_cov),mean(bin_y0_cov),
            mean(evd_symm_y0_cov_u),mean(bin_y0_cov_u),
            mean(pt_est_b0-boot_y0_low)/t_adj,
            mean(pt_est_b0-evd_symm_y0_low_w)/t_adj,
            mean(pt_est_b0-evd_symm_y0_low)/t_adj,
            mean(pt_est_b0-bin_y0_low_w)/t_adj,
            mean(med_pt_est_x),
            mean(med_pt_est_x_wgt),mean(rsq),mean(samp_r1),
            mean(med_pt_est_x2),
            mean(med_pt_est_x2_wgt),
            l2_mk[j],mean(pt_est_b2),
            mean(boot_cov2),mean(evd_symm_cov2),mean(bin_cov2),
            mean(pt_est_b2-boot_low2)/t_adj,
            mean(pt_est_b2-evd_symm_low2)/t_adj,
            mean(pt_est_b2-bin_low2)/t_adj,
            mean(pt_est_b2-evd_low2)/t_adj,
            mean(pt_est_b2-evd_high2)/t_adj)

output[j,] <- yyers
#ifelse(j > 1,app_flag <- TRUE,app_flag <- FALSE)


}

print(output)
write.csv(output,file=paste("/home/bear/Projects/quantile_regression/",
                            sampsize,"_",filename,".csv"))


output[2,]

# h_list <- 0
# for (i in 1:1000) h_list[i] <- (sqrt(var(exp(rnorm(sampsize,7.5,.5)^1))))
# hist(h_list,breaks=100)
# mean(h_list)
# sqrt(var(exp(rnorm(1000000,7.5,.5)^1)))
