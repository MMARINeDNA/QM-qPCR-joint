data {
  
  int Nplates;
  int Nobs;
  int NSamples; //number of unique biol samples, overall
  int NstdSamples; //number of unique biol samples with known concentrations (standards)
  int plate_idx[Nobs];
  int std_idx[NstdSamples]; //index relative to NSamples; which ones are the standards?
  int unkn_idx[NSamples-NstdSamples];
  int plateSample_idx[Nobs]; //index of unique combinations of plate and biological sample
  
  
  vector[Nobs] y; //Ct observations
  int z[Nobs]; //indicator; z = 1 if a Ct was observed, z = 0 otherwise
  vector[NstdSamples] known_concentration;
  
  real stdCurvePrior_intercept[2];
  real stdCurvePrior_slope[2];

  //vector[NSamples-NstdSamples] dilutionFactor;
  
}

transformed data {
  vector[NstdSamples] known_conc;
    known_conc = log10(known_concentration);
}

parameters {
  vector<lower=0>[Nplates] beta_std_curve_0;
  vector<upper=0>[Nplates] beta_std_curve_1;
  real<upper=0> gamma_0; //intercept to scale variance w the mean
  // real<lower=0> gamma_1; //slope to scale variance w the mean
  // vector<upper=0>[Nplates]  gamma_0; //intercept to scale variance w the mean
  vector<upper=0>[Nplates]  gamma_1; //slope to scale variance w the mean
  // vector[Nplates] phi_0;
  // real phi_0;
  // vector[Nplates] phi_1;
  real phi_1;
  // real phi_0;
  // real phi_1;
  // real beta_std_curve_0_hyperMean;
  // real beta_std_curve_1_hyperMean;
  // real<lower=0> hyperSigma_0;
  // real<lower=0> hyperSigma_1;
  
  // real<upper=0> gamma_0_hyperMean;
  // real gamma_1_hyperMean;
  // real<lower=0> gamma_hyperSigma_0;
  // real gamma_hyperSigma_1;

  // real phi_1_hyperMean;
  // real phi_hyperSigma_1;

  
  vector[NSamples-NstdSamples] envir_concentration;
  
}

transformed parameters{
  vector[NSamples] Concentration;
  vector[NSamples] mu;
  vector[NSamples] sigma;
  vector[NSamples] theta; //bernoulli param

  //slot knowns and unknowns into a common vector, where unknowns are treated as params and knowns are treated as data 
  Concentration[std_idx] = known_conc; //log10 scale
  Concentration[unkn_idx] = envir_concentration;
  
  for(i in 1:Nobs){
    mu[plateSample_idx[i]] = beta_std_curve_0[plate_idx[i]] + 
                              beta_std_curve_1[plate_idx[i]] * Concentration[plateSample_idx[i]];

//may need to change the scale of gamma_0; model can't disting low from very low

    sigma[plateSample_idx[i]] = exp(gamma_0  
                                       + gamma_1[plate_idx[i]]*Concentration[plateSample_idx[i]]);
                                  
    theta[plateSample_idx[i]] = inv_logit(phi_1*Concentration[plateSample_idx[i]]); //phi_0 +
                                  
  }
}


model {
   for(i in 1:Nobs){
     z[i]   ~ bernoulli(theta[plateSample_idx[i]]);
      // z[i]   ~ bernoulli( inv_logit(theta[plateSample_idx[i]]) ) ;
    }
    
    for(i in 1:Nobs){
      if (z[i]==1){ //if Ct observed, then compute likelihood
        y[i] ~ normal(mu[plateSample_idx[i]], sigma[plateSample_idx[i]]);   
      }
    }

  //beta params hierarchical
  // beta_std_curve_0_hyperMean ~ normal(stdCurvePrior_intercept[1], stdCurvePrior_intercept[2]);
  // beta_std_curve_1_hyperMean ~ normal(stdCurvePrior_slope[1], stdCurvePrior_slope[2]);
  // hyperSigma_0 ~ gamma(.5,.5);
  // hyperSigma_1 ~ gamma(1,1);
  beta_std_curve_0 ~ normal(stdCurvePrior_intercept[1], stdCurvePrior_intercept[2]);
  beta_std_curve_1 ~ normal(stdCurvePrior_slope[1], stdCurvePrior_slope[2]);
  
  //gamma params hierarchical
  // gamma_0_hyperMean ~ normal(-5,1);
  // gamma_1_hyperMean ~ normal(0,3);
  // gamma_hyperSigma_0 ~ gamma(1,1);
  // gamma_hyperSigma_1 ~ gamma(1,1);
  // gamma_0 ~ normal(gamma_0_hyperMean,gamma_hyperSigma_0);
  gamma_1 ~ normal(0,5);

  gamma_0 ~ normal(-2,1);
  // gamma_1 ~ normal(0,1);
  

  
  envir_concentration ~ normal(0, 2); //log10 scale
  
  //phi can be shared; doesn't need hierarchy
  // phi_1_hyperMean ~ normal(5, 2);
  // phi_hyperSigma_1 ~ gamma(1, 1);
  // phi_1 ~ normal(phi_1_hyperMean, phi_hyperSigma_1);
  
  phi_1 ~ normal(5, 2);  
  // phi_0 ~ normal(0, 2);
  
}

generated quantities {

  vector[Nplates] efficiency;

  for (i in 1:Nplates){
    efficiency[i] = 10^(-1 / beta_std_curve_1[i]) - 1;  
  }
  
}


