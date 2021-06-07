N=100000
sim_pat = function(p_U=.001,beta_U=.45,beta_last_U=.45,p_A=.5,int=0,hr=2,hours=24){
  U= rbinom(hours,1,p_U)
  if(sum(U)==0){
    return(c(NA,NA))
  }
  A= rbinom(hours,1,p_A)
  last_U = c(0,U[1:(hours-1)])
  Y_0 = rbinom(hours,1,pmin(1/hr,beta_U*U + beta_last_U*last_U + int))
  Y_1 = rbinom(hours,1,hr*pmin(1/hr,beta_U*U + beta_last_U*last_U + int))
  Y = A*Y_1 + (1-A)*Y_0
  if(sum(Y)==0){
    return(c(NA,NA))
  }else{
    time = min(which(Y==1))
    if(time==1){
      return(c(NA,NA))
    }else{
      return(c(A[time-1],A[time]))
    }
  }
}
ests = rep(NA,100)
for(i in 1:100){
  data = replicate(N,sim_pat())
  D = t(data.frame(data))
  D = D[!is.na(D[,1]),]
  mh = sum(D[,2]==1&D[,1]==0)/sum(D[,2]==0&D[,1]==1)
  ests[i] = mh
}
