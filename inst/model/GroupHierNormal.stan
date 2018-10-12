data {
  int<lower=1> N_obs; //total number of obs per season-region
  int<lower=1> N_subj; //total number of season-region combination
  int<lower=1> dim_space; // number of eigen functions
  int<lower=1> N_group; // number of groups
  real W[N_subj, N_obs]; // observed value of flu
  real E[N_obs, dim_space]; // eigenfns of smoothed data
  int<lower=1> group[N_subj]; // group indicator variable
}

parameters {
  real<lower=0> eps; // SD for data model
  real<lower=0> alphaL; // SD for t prior
  vector<lower=0>[N_group] tau; // SD for beta or Global Shrinkage Param
  vector<lower=0>[dim_space] lamb[N_group]; // SD for beta or Local Shrinkage Param
  vector[dim_space] beta[N_group];
}

transformed parameters {
  vector[N_obs] m[N_subj]; // Underlying smooth function
  for(i in 1:N_subj){
    for(t in 1:N_obs){
      m[i,t] = beta[group[i],1]*E[t,1]+
	             beta[group[i],2]*E[t,2]+
	             beta[group[i],3]*E[t,3]+
	             beta[group[i],4]*E[t,4]+
	             beta[group[i],5]*E[t,5]+
	             beta[group[i],6]*E[t,6]+
	             beta[group[i],7]*E[t,7]+
	             beta[group[i],8]*E[t,8]+
	             beta[group[i],9]*E[t,9]+
	             beta[group[i],10]*E[t,10]+
	             beta[group[i],11]*E[t,11]+
	             beta[group[i],12]*E[t,12]+
	             beta[group[i],13]*E[t,13]+
	             beta[group[i],14]*E[t,14]+
	             beta[group[i],15]*E[t,15]+
	             beta[group[i],16]*E[t,16]+
	             beta[group[i],17]*E[t,17]+
	             beta[group[i],18]*E[t,18]+
	             beta[group[i],19]*E[t,19]+
	             beta[group[i],20]*E[t,20]+
	             beta[group[i],21]*E[t,21]+
	             beta[group[i],22]*E[t,22]+
	             beta[group[i],23]*E[t,23]+
	             beta[group[i],24]*E[t,24]+
	             beta[group[i],25]*E[t,25]+
	             beta[group[i],26]*E[t,26]+
	             beta[group[i],27]*E[t,27]+
	             beta[group[i],28]*E[t,28]+
	             beta[group[i],29]*E[t,29]+
	             beta[group[i],30]*E[t,20];
    }
  }
}

model {
  
  eps ~ student_t(4, 0, 1); //Prior on model SD
  alphaL ~ inv_gamma(1, 1); //Prior on std of t dist
  
  
  for(g in 1:N_group){
    for(j in 1:dim_space){
      lamb[g,j] ~ student_t(4, 0, alphaL); //Prior on  SD
    }
  }
  
  
  for(i in 1:N_subj){
    for(t in 1:N_obs){
      W[i,t]~normal(m[i,t], eps);
    }
    
    for (j in 1:dim_space){
      beta[group[i],j] ~ normal(0, lamb[group[i],j]);
	  }

  }
}

