get_limits = function(lambda_g0_0 = .001,effect0=2,lambda_g1_0 = .0005,effect_delta=5,
                      p_A_g0 = .8,p_A_g1 = .5,p_g = .5){
  lambda_g0_1 = lambda_g0_0*effect0
  lambda_g1_1 = lambda_g1_0*effect0*effect_delta
  cc_limit = (lambda_g1_1*p_A_g1*(1 - p_A_g1)*p_g + lambda_g0_1*p_A_g0*(1 - p_A_g0)*(1-p_g))/
    (lambda_g1_0*p_A_g1*(1 - p_A_g1)*p_g + lambda_g0_0*p_A_g0*(1 - p_A_g0)*(1-p_g))
  
  rct_limit = (lambda_g1_1*p_g + lambda_g0_1*(1-p_g))/
    (lambda_g1_0*p_g + lambda_g0_0*(1-p_g))
  
  cohort_est_limit =  ((lambda_g1_1*p_A_g1*p_g + lambda_g0_1*p_A_g0*(1-p_g))/(p_A_g1*p_g + p_A_g0*(1-p_g)))/
    ((lambda_g1_0*(1-p_A_g1)*p_g + lambda_g0_0*(1-p_A_g0)*(1-p_g))/((1-p_A_g1)*p_g + (1-p_A_g0)*(1-p_g)))
  
  return(c(cc_limit,rct_limit,cohort_est_limit))
}

#only effect 0, effect 1, p_A_g_0, and p_A_g1 impact heterogeneity bias from case crossover
#but comparison to bias of a naive cohort estimator may also be of interest, and that may also depend on lambda_gi_j and p_g

#independently vary: lambda_g0_0, lambda_g1_0, effect0, effect1/effect0,p_A_g_0,p_A_g1,p_g
#lambda_g0_0, lambda_g1_0: .0005,.001
#effect1/effect0: 1-10
#effect0: 1-5
#p_A_g0, p_A_g1: .05 - 1
lambda_g0_0 = c(.0005,.001)
effect0 = 1:5
lambda_g1_0 = c(.0005,.001)
effect_delta = 1:10
p_A_g0 = (1:19)*.05
p_A_g1 = (1:19)*.05

sensitivity_grid = expand.grid(list(lambda_g0_0=lambda_g0_0,effect0=effect0,
                                    lambda_g1_0=lambda_g1_0,effect_delta=effect_delta,
                                    p_A_g0=p_A_g0,p_A_g1=p_A_g1))

results = lapply(1:nrow(sensitivity_grid),function(i) get_limits(lambda_g0_0 = sensitivity_grid$lambda_g0_0[i],effect0=sensitivity_grid$effect0[i],lambda_g1_0 = sensitivity_grid$lambda_g1_0[i],
                                                                 effect_delta=sensitivity_grid$effect_delta[i],
                                                                 p_A_g0 = sensitivity_grid$p_A_g0[i],p_A_g1 = sensitivity_grid$p_A_g1[i],p_g = .5))

results = t(data.frame(results))
results = data.frame(results)
names(results) = c('CC','RCT','Cohort')
rownames(results) = NULL
sensitivity_grid = cbind(sensitivity_grid,results)

par(mfrow=c(1,3))
library(scales)
plot(sensitivity_grid$CC/sensitivity_grid$RCT, sensitivity_grid$Cohort/sensitivity_grid$RCT,col=alpha("red",.05),
     xlab='Multiplicative Case Crossover Bias',ylab = "Multiplicative Cohort Bias",ylim=c(0,2.6),xlim=c(0,2.5))
library(MASS)
truehist(sensitivity_grid$CC[sensitivity_grid$effect_delta==10 & abs(abs(sensitivity_grid$p_A_g1-.5)-abs(sensitivity_grid$p_A_g0-.5))>=.1]/sensitivity_grid$RCT[sensitivity_grid$effect_delta==10 & abs(abs(sensitivity_grid$p_A_g1-.5)-abs(sensitivity_grid$p_A_g0-.5))>=.1],xlab="Multiplicative Case-Crossover Bias",main="")
truehist((sensitivity_grid$CC/sensitivity_grid$RCT)/(sensitivity_grid$Cohort/sensitivity_grid$RCT),xlab="Multiplicative Case-Crossover Bias / Multiplicative Cohort Bias",main="")


rownames(sensitivity_grid) <- c()
sensitivity_grid$cc_bias = sensitivity_grid$CC/sensitivity_grid$RCT
library(dplyr)
library(plsgenomics)
for(i in 1:length(effect0)){
  for(j in 1:length(effect_delta)){
    cc_only = distinct(sensitivity_grid[sensitivity_grid$effect0==effect0[i] & sensitivity_grid$effect_delta==effect_delta[j] & sensitivity_grid$lambda_g0_0==.0005 & sensitivity_grid$lambda_g1_0==.0005,c("effect0","effect_delta","p_A_g0","p_A_g1","cc_bias")])
    x = matrix(cc_only$cc_bias,nrow=20)
    heatmap.2( x, Rowv=FALSE, Colv=FALSE, dendrogram='none')
  }
}
