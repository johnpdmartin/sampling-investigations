library(Runuran)
library(quantreg)
sampsize <- 1000
replicates <- 1000
x_samp <- matrix(rep(1,sampsize*replicates),nrow=replicates,ncol=sampsize)
x_samp2 <- matrix(rep(1,sampsize*replicates),nrow=replicates,ncol=sampsize)
y_samp <- matrix(rep(1,sampsize*replicates),nrow=replicates,ncol=sampsize)
med_pt_est_x <- matrix(rep(1,sampsize*replicates),nrow=replicates,ncol=sampsize)
med_pt_est_x2 <- matrix(rep(1,sampsize*replicates),nrow=replicates,ncol=sampsize)
output <- matrix(rep(1,34*15),nrow=15,ncol=34)
filename <- "twoterms_hetguassian_1_dfadj_in_perc_scale" #  _perc_scale" # _org_scale" #
set.seed(592)
# set quantile point of interest
exp_theta <- 0.5
method_df <- "perc_scale" # "perc_scale" # "org_scale" #  
multi <- 1
multi2 <- 0.5
lambda <- 1;shift <- -0;pcon <- 0
x_sample <- runif(sampsize,-180,180)^1+(shift)#exp(rnorm(sampsize,0,2))#runif(sampsize,-10,10)#seq(0.025,.975,length.out=20)
x_sample2 <- (x_sample-shift)^2
y_sample <- multi*(x_sample-shift)+multi2*x_sample2+(10+x_sample-(shift))^pcon*rnorm(sampsize,0,lambda)#1*urlaplace(sampsize,0,lambda) #
n_obs <- length(x_sample)
t_adj <- abs(qt(.025,n_obs))
#quantile(urlaplace(1000000,0,lambda),taus)


#lambda <- 600

model_fit_rq <- rq((y_sample)~(x_sample+x_sample2),tau=exp_theta)
model_fit_rq
rq_beta0 <- model_fit_rq$coefficients[1]
rq_beta1 <- model_fit_rq$coefficients[2]
rq_beta2 <- model_fit_rq$coefficients[3]
summary(model_fit_rq,se="boot",R=600)

plot(y=y_sample,x=x_sample)
points(y=(rq_beta0+rq_beta1*x_sample+rq_beta2*x_sample2),x=x_sample,col="red")

# full mode 95% CI

# START MANUAL SETTING 1

#lambda <- 1
sampsize <- length(y_sample)
# taus_set <- c(0.025,.05,.1,.2,.25,.3,.4,.5,.6,.7,.75,.8,.9,.95,.975)
taus_set <- c(exp_theta)

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
biny0_cov <- 0

# START MANUAL SETTING 2

for (reps_setup in 1:replicates) {
  
  x_samp[reps_setup,] <- runif(sampsize,-180,180)^1+(shift)#exp(rnorm(sampsize,0,2))#runif(sampsize,-10,10)#seq(0.025,.975,length.out=20)
  x_samp2[reps_setup,] <- (x_samp[reps_setup,]-shift)^2
  y_samp[reps_setup,] <- multi*(x_samp[reps_setup,]-shift)+multi2*x_samp2[reps_setup,]+(10+x_samp[reps_setup,]-(shift))^pcon*rnorm(sampsize,0,lambda) #

  med_est_x <- rq((x_samp[reps_setup,])~(1),tau=.5,weights=wgt_filereg)
  med_pt_est_x[reps_setup,] <- med_est_x$coefficients[1]

    med_est_x2 <- rq((x_samp2[reps_setup,])~(1),tau=.5,weights=wgt_filereg)
  med_pt_est_x2[reps_setup,] <- med_est_x2$coefficients[1]
  
}  
  
    
for (j in 1:15)  {
#for (j in 1:1)  {
  
  # END MANUAL SETTING 2
  
  taus <- taus_set[j]
  
  boot_cov <- 0
  evd_symm_cov <- 0
  bin_cov <- 0
  
  boot_cov2 <- 0
  evd_symm_cov2 <- 0
  bin_cov2 <- 0
  
  boot_y0_cov <- 0
  evd_symm_y0_cov <- 0
  bin_y0_cov <- 0
  
  
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
  bin_y0_low <- 0
  bin_y0_high <- 0
  
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
  
  
  # START MANUAL SETTING 3
  
  set.seed(592)
  x_sample <- runif(sampsize,-180,180)^1+(shift)#exp(rnorm(sampsize,0,2))#runif(sampsize,-10,10)#seq(0.025,.975,length.out=20)
  x_sample2 <- (x_sample-shift)^2
  y_sample <- multi*(x_sample-shift)+multi2*x_sample2+(10+x_sample-(shift))^pcon*rnorm(sampsize,0,lambda) #
  model_fit_rq <- rq((y_sample)~(x_sample+x_sample2),tau=taus,weights=wgt_filereg)
  model_fit_rq
  rq_beta0 <- model_fit_rq$coefficients[1]
  rq_beta1 <- model_fit_rq$coefficients[2]
  rq_beta2 <- model_fit_rq$coefficients[3]
  summary(model_fit_rq,se="boot",R=600)
  samp <- y_sample-(rq_beta0+rq_beta1*x_sample+rq_beta2*x_sample2)
  samp_int <- samp
  
  # STOP MANUAL SETTING 3
  
  
  n_samp <- length(samp)
  
  # quantreg estimates of X distribution values for given $\theta$ values
  
  
  # START MANUAL SETTING 4
  
  
  # l_mk <- qpois(taus,lambda)
  
  #l_mk <- qnorm(taus,0,lambda)
  l_mk <- multi
  l2_mk <- multi2

  x_large <- c(x_samp)
  x2_large <- c(x_samp2)
  y_large <- c(y_samp)
  
  
  bmk_model <- rq(y_large~x_large+x2_large,tau=taus,method="fn")
  print(bmk_model) 
  m_mk <- bmk_model$coefficient[1]
  l_mk <- bmk_model$coefficient[2]
  l2_mk <- bmk_model$coefficient[3]
  print("estimated popn regression slopes and intercept")
  print(taus_set)
  print(round(m_mk,3))
  print(round(l_mk,3))
  print(round(l2_mk,3))
  
  
  # END MANUAL SETTING 4
  
if(method_df == "perc_scale") {
  q_var <- sqrt(taus*(1-taus)/(n_samp-3))} else
    {q_var <- sqrt(taus*(1-taus)/(n_samp))}
  
  
  ### df correction in original frame
  #  q_var <- sqrt(taus*(1-taus)/(n_samp))
  
  ### df correction in percentile frame
#  q_var <- sqrt(taus*(1-taus)/(n_samp-2))
  
  p_low <- qnorm(0.025,taus,q_var)
  p_high <- qnorm(0.975,taus,q_var)

  if(method_df == "perc_scale")  {b_low <- qbinom(.025,sampsize-3,taus)/(sampsize-3)
  b_high <- qbinom(.975,sampsize-3,taus)/(sampsize-3)}  else  {b_low <- qbinom(.025,sampsize,taus)/(sampsize)
  b_high <- qbinom(.975,sampsize,taus)/(sampsize)}
  
  ### df correction in original frame
  #   b_low <- qbinom(.025,sampsize,taus)/(sampsize)
  #   b_high <- qbinom(.975,sampsize,taus)/(sampsize)
  
  ### df correction in percentile frame
#   b_low <- qbinom(.025,sampsize-2,taus)/(sampsize-2)
#   b_high <- qbinom(.975,sampsize-2,taus)/(sampsize-2)
  
    p_low <- max(p_low,0)
    p_high <- min(p_high,1)
  
  
  print("binomial CI | evd CI  , in percentile scale")
  print("lower bounds")
  print(paste(b_low,p_low))
  print("upper bounds")
  print(paste(b_high,p_high))
  
  # START MANUAL SETTING 5
  
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
    samp <- y_sample-(rqr_beta0+rqr_beta1*x_sample+rqr_beta2*x_sample2)
    samp_var[reps] <- var(samp)
    y_var[reps] <- var(y_sample)
    
    rho <- function(u,taur=.5) u*(taur - (u < 0))
    samp_r1[reps] <- 1 - model_fit_rq$rho/rq(y_sample~1,tau=taus,weights=wgt_filereg)$rho
    
    # END MANUAL SETTING 5   
    
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
      
      model_fit <- rq((samp)~(1),tau=taus[i],method="fn",weights=wgt_filereg)
      #      print(summary(model_fit,se="boot",R=600))
      qr_dat[i,] <- c(summary(model_fit)$coefficients[1])
      qrboot_dat[i,] <- c(summary(model_fit,se="boot",R=600)$coefficients[1:4])
      
      qp_lstd[i] <- taus[i]-q_var[i]
      qp_ustd[i] <- taus[i]+q_var[i]
      taustd <- c(qp_lstd[i],qp_ustd[i])
      
      #       plot(y=samp,x=y_sample)
      #       mod2 <- rq((samp)~(y_sample),tau=taus[i],method="fn",weights=wgt_filereg)
      #       
      model_fit_lstd <- summary(rq((samp)~(1),tau=max(qp_lstd[i],0.00005),method="fn",weights=wgt_filereg))$coefficients[1]
      model_fit_ustd <- summary(rq((samp)~(1),tau=min(qp_ustd[i],.99995),method="fn",weights=wgt_filereg))$coefficients[1]
      
      backtrans_std[i,] <- c(model_fit_lstd-qr_dat[i],model_fit_ustd-qr_dat[i],(model_fit_ustd-model_fit_lstd)/2)
      
      
      
    }
    

    back_pts1 <- rq((samp)~(1),tau=p_low,weights=wgt_filereg)$coefficients[1]
    
    back_pts2 <- rq((samp)~(1),tau=p_high,weights=wgt_filereg)$coefficients[1]
    
    back_pts3 <- summary(rq((samp)~(1),tau=b_low,weights=wgt_filereg))$coefficients[1]
    
    back_pts4 <- summary(rq((samp)~(1),tau=b_high,weights=wgt_filereg))$coefficients[1]
    
    
    
    
    bk_pt <- c(back_pts1,back_pts2,back_pts3,
               back_pts4)
    

           s2x1x2 <- sum(x_sample*x_sample2)-sum(x_sample)*sum(x_sample2)/n_samp
    s2x1 <- sum(x_sample^2)-sum(x_sample)^2/n_samp
           s2x2 <- sum(x_sample2^2)-sum(x_sample2)^2/n_samp
       
           r12 <- s2x1x2/sqrt(s2x1*s2x2)
    # r12 <- 0.935 
    
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
    #    evd_ci <- abs((bk_pt[1]-bk_pt[2])/2)
    evd_low[reps] <- rqr_beta1+bk_pt[1]/se1adj
    evd_high[reps] <- rqr_beta1+bk_pt[2]/se1adj
    evd_ci <- mean(abs(bk_pt[1]),abs(bk_pt[2]))
    evd_symm_low[reps] <- rqr_beta1-evd_ci/se1adj
    evd_symm_high[reps] <- rqr_beta1+evd_ci/se1adj
    #    bin_ci <- abs((bk_pt[3]-bk_pt[4])/2)
    bin_ci <- mean(abs(bk_pt[3]),abs(bk_pt[4]))
    bin_low[reps] <- rqr_beta1-bin_ci/se1adj
    bin_high[reps] <- rqr_beta1+bin_ci/se1adj
    
    ifelse(l_mk >= boot_low[reps] & l_mk <= boot_high[reps],
           boot_cov[reps] <-1, boot_cov[reps] <- 0)
    
    ifelse(l_mk >= evd_symm_low[reps] & l_mk <= evd_symm_high[reps],
           evd_symm_cov[reps] <- 1,evd_symm_cov[reps] <- 0)
    
    ifelse(l_mk >= bin_low[reps] & l_mk <= bin_high[reps],
           bin_cov[reps] <- 1,bin_cov[reps] <- 0)
    
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
    
    ifelse(l2_mk >= boot_low2[reps] & l2_mk <= boot_high2[reps],
           boot_cov2[reps] <-1, boot_cov2[reps] <- 0)
    
    ifelse(l2_mk >= evd_symm_low2[reps] & l2_mk <= evd_symm_high2[reps],
           evd_symm_cov2[reps] <- 1,evd_symm_cov2[reps] <- 0)
    
    ifelse(l2_mk >= bin_low2[reps] & l2_mk <= bin_high2[reps],
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
    {denxvar2 <- sqrt(1+(med_pt_est_x[reps])^2/(se1adj)^2+
        (med_pt_est_x2[reps])^2/(se2adj)^2)}
      else
      {denxvar2 <- sqrt(n_samp/(n_samp-3)+(shift-med_pt_est_x[reps])^2/(se1adj)^2+
        (shift^2-med_pt_est_x2[reps])^2/(se2adj)^2)}
    

        ### df correction in original frame
    #        denxvar2 <- sqrt(n_samp/(n_samp-2))
    
    ### df correction in percentile frame
    #denxvar2 <- 1
    
    
    #    denxboot2 <- 1
    #    denxboot2 <- sqrt(1+(pt_est_x[reps])^2/((uncond_evd_high_x2std[reps]-
    #                                               uncond_evd_low_x2std[reps])/2/t_adj)^2/(n_samp-1)/(n_samp-1))
    denxevd2 <- denxvar2
    
    denxbin2 <- denxvar2
    
    boot_y0_low[reps] <- rqr_beta0-t_adj*beta0_std
    boot_y0_high[reps] <- rqr_beta0+t_adj*beta0_std
    evd_y0_low[reps] <- rqr_beta0+bk_pt[1]*denxevd2
    evd_y0_high[reps] <- rqr_beta0+bk_pt[2]*denxevd2
    evd_symm_y0_low[reps] <- rqr_beta0-evd_ci*denxevd2
    evd_symm_y0_high[reps] <- rqr_beta0+evd_ci*denxevd2
    bin_y0_low[reps] <- rqr_beta0-bin_ci*denxbin2
    bin_y0_high[reps] <- rqr_beta0+bin_ci*denxbin2

    if(method_df == "perc_scale") {  samp_var_theta[reps] <- (bin_ci/t_adj)^2}
     else 
      {samp_var_theta[reps] <- (bin_ci/t_adj)^2*n_samp/(n_samp-3)}


        ### df correction in original frame
#    samp_var_theta[reps] <- (bin_ci/t_adj)^2*n_samp/(n_samp-2)
    
    ### df correction in percentile frame
#    samp_var_theta[reps] <- (bin_ci/t_adj)^2
    
    
    ifelse(m_mk >= boot_y0_low[reps] & m_mk <= boot_y0_high[reps],
           boot_y0_cov[reps] <-1, boot_y0_cov[reps] <- 0)
    
    ifelse(m_mk >= evd_symm_y0_low[reps] & m_mk <= evd_symm_y0_high[reps],
           evd_symm_y0_cov[reps] <- 1,evd_symm_y0_cov[reps] <- 0)
    
    ifelse(m_mk >= bin_y0_low[reps] & m_mk <= bin_y0_high[reps],
           bin_y0_cov[reps] <- 1,bin_y0_cov[reps] <- 0)
    
    
  }
  
  bcov[j] <- mean(boot_cov)
  n_reps <- length(boot_cov)
  bincov[j] <- mean(bin_cov)
  e_symmcov[j] <- mean(evd_symm_cov)
  
  print("")
  print("coverage performance beta1 and mean quantile regression estimate")
  print("q pt | bootstrap | binomial | evdsymm | quantreg est | pop value")
  print(paste(taus,bcov[j],bincov[j],e_symmcov[j],mean(pt_est_b1),l_mk))
  
  bcov2[j] <- mean(boot_cov2)
  #n_reps <- length(boot_cov)
  bincov2[j] <- mean(bin_cov2)
  e_symmcov2[j] <- mean(evd_symm_cov2)
  
  print("")
  print("coverage performance beta1 and mean quantile regression estimate")
  print("q pt | bootstrap | binomial | evdsymm | quantreg est | pop value")
  print(paste(taus,bcov2[j],bincov2[j],e_symmcov2[j],mean(pt_est_b2),l2_mk))
  
  by0_cov[j] <- mean(boot_y0_cov)
  n_reps <- length(boot_y0_cov)
  biny0_cov[j] <- mean(bin_y0_cov)
  e_symmy0_cov[j] <- mean(evd_symm_y0_cov)
  
  print("")
  print("coverage performance beta0 and mean quantile regression estimate")
  print("q pt | bootstrap | binomial | evdsymm | quantreg est | pop value")
  print(paste(taus,by0_cov[j],biny0_cov[j],e_symmy0_cov[j],mean(pt_est_b0),m_mk))
  
  
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

sum(ifelse(boot_high < multi,1,0))/length(boot_high)
sum(ifelse(boot_low > multi,1,0))/length(boot_high)
sum(ifelse(evd_symm_high < multi,1,0))/length(boot_high)
sum(ifelse(evd_symm_low > multi,1,0))/length(boot_high)
sum(ifelse(bin_high < multi,1,0))/length(boot_high)
sum(ifelse(bin_low > multi,1,0))/length(boot_high)

hist(pt_est_b1,breaks=100)
hist(boot_high,breaks=100)
hist(boot_low,breaks=100)
hist(evd_symm_high,breaks=100)
hist(evd_symm_low,breaks=100)
hist(bin_high,breaks=100)
hist(bin_low,breaks=100)



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
            
sum(ifelse(boot_y0_high < m_mk,1,0))/length(boot_high)
sum(ifelse(boot_y0_low > m_mk,1,0))/length(boot_high)
sum(ifelse(evd_symm_y0_high < m_mk,1,0))/length(boot_high)
sum(ifelse(evd_symm_y0_low > m_mk,1,0))/length(boot_high)
sum(ifelse(bin_y0_high < m_mk,1,0))/length(boot_high)
sum(ifelse(bin_y0_low > m_mk,1,0))/length(boot_high)

hist(pt_est_b0,breaks=100)
hist(boot_y0_high,breaks=100)
hist(boot_y0_low,breaks=100)
hist(evd_symm_y0_high,breaks=100)
hist(evd_symm_y0_low,breaks=100)
hist(bin_y0_high,breaks=100)
hist(bin_y0_low,breaks=100)

mea_lm <- 0
for (i in 1:replicates) mea_lm[i] <- summary(lm(y_samp[i,]~x_samp[i,]+x_samp2[i,]))$r.squared
mean(mea_lm)

rsq <- 1- samp_var/y_var;print(mean(rsq))
#rsq_theta <- 1- samp_var_theta/y_var_theta;print(mean(rsq_theta))
print(mean(samp_r1))
hist(rsq,breaks=30);#hist(rsq_theta,breaks=30);
hist(samp_r1,breaks=30)
hist(pt_est,breaks=30)

print(j)
output[j,] <- c(taus,l_mk,mean(pt_est_b1),
            mean(boot_cov),mean(evd_symm_cov),mean(bin_cov),
            mean(pt_est_b1-boot_low)/t_adj,
            mean(pt_est_b1-evd_symm_low)/t_adj,
            mean(pt_est_b1-bin_low)/t_adj,
            mean(pt_est_b1-evd_low)/t_adj,
            mean(pt_est_b1-evd_high)/t_adj,
            m_mk,mean(pt_est_b0),
            mean(boot_y0_cov),mean(evd_symm_y0_cov),mean(bin_y0_cov),
            mean(pt_est_b0-boot_y0_low)/t_adj,
            mean(pt_est_b0-evd_symm_y0_low)/t_adj,
            mean(pt_est_b0-bin_y0_low)/t_adj,
            mean(pt_est_b0-evd_y0_low)/t_adj,
            mean(pt_est_b0-evd_y0_high)/t_adj,
            mean(mea_lm),mean(rsq),mean(samp_r1),
            l2_mk,mean(pt_est_b2),
            mean(boot_cov2),mean(evd_symm_cov2),mean(bin_cov2),
            mean(pt_est_b2-boot_low2)/t_adj,
            mean(pt_est_b2-evd_symm_low2)/t_adj,
            mean(pt_est_b2-bin_low2)/t_adj,
            mean(pt_est_b2-evd_low2)/t_adj,
            mean(pt_est_b2-evd_high2)/t_adj)


#ifelse(j > 1,app_flag <- TRUE,app_flag <- FALSE)


}

print(output)
write.csv(output,file=paste("/home/bear/Projects/quantile_regression/",
                            sampsize,"_",filename,".csv"))


output[2,]
mean(med_pt_est_x)
mean(med_pt_est_x2)