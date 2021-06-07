N=100000
sim_pat = function(p_U=.001,beta_U=1-.55^(1/3600),beta_last_U=1-.55^(1/3600),p_A=.5,int=0,hr=3.85,hours=24,Time=3600*hours){
  U_hours = rbinom(hours,1,p_U)
  U = rep(U_hours,each=3600)
  if(sum(U)==0){
    return(c(NA,NA))
  }
  A_hours = rbinom(hours,1,p_A)
  A = rep(A_hours,each=3600)
  last_U_hours = c(0,U_hours[1:(hours-1)])
  last_U = rep(last_U_hours,each=3600)
  Y_0 = rbinom(Time,1,beta_U*U + beta_last_U*last_U + int)
  Y_1 = rbinom(Time,1,hr*(beta_U*U + beta_last_U*last_U + int))
  Y = A*Y_1 + (1-A)*Y_0
  if(sum(Y)==0){
    return(c(NA,NA))
  }else{
    time = min(which(Y==1))
    if(time<=3600){
      return(c(NA,NA))
    }else{
      return(c(A[time-3600],A[time]))
    }
  }
}
ests = rep(NA,100)
for(i in 1:100){
  data = replicate(N,sim_pat())
  mh = sum(data[2,]==1&data[1,]==0,na.rm=T)/sum(data[2,]==0&data[1,]==1,na.rm=T)
  ests[i] = mh
}

U_hours = rbinom(hours,1,p_U)
U = rep(U_hours,each=3600)
sum(U)
