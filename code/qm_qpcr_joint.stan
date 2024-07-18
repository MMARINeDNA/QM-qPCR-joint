data { 
  // DATA FOR qPCR PART OF THE MODEL
  
  int Nplates; // number of PCR plates
  int Nobs_qpcr; // number of field observations for qPCR
  int NSamples_qpcr; //number of unique biol samples, overall
  int NstdSamples; //number of unique biol samples with known concentrations (standards)
  int plate_idx[Nobs_qpcr]; //index denoting which PCR plate each sample is on
  int std_idx[NstdSamples]; //index relative to NSamples; which ones are the standards?
  int unkn_idx[NSamples_qpcr-NstdSamples]; //index relative to total samples; which ones are the unknown/field samples?
  int plateSample_idx[Nobs_qpcr]; //index of unique combinations of plate and biological sample
 
  vector[Nobs_qpcr] y; //Ct observations
  int z[Nobs_qpcr]; //indicator; z = 1 if a Ct was observed, z = 0 otherwise
  vector[NstdSamples] known_concentration; //known concentration (copies/vol) in standards
  
  real stdCurvePrior_intercept[2]; // prior on the intercept of the std curves
  real stdCurvePrior_slope[2]; // prior on the slope of the std curves

  //vector[NSamples-NstdSamples] dilutionFactor;
  
  //END DATA FOR qPCR
  
  // DATA FOR METABARCODING PART OF THE MODEL
  int N_species; // Number of species in data
  int N_obs_mb_samp;  // Number of observed samples, also the number of groups for qPCR samps to link to MB samps
  int N_obs_mb_samp_small;  // Number of observed samples for individual sites.
  int N_obs_mock; // Number of observed mock samples

  // Observed data of community matrices
  int sample_data[N_obs_mb_samp,N_species];
  // Observed data of mock community matrices
  int mock_data[N_obs_mock,N_species];  
  // True proportions for mock community
  matrix[N_obs_mock,N_species] alr_mock_true_prop ;
 // matrix[N_obs_mock_small,N_species] alr_mock_true_prop_small ;
    
  // Design matrices: field samples
  int N_b_samp_col; // Number of species (except the reference?)
  matrix[N_obs_mb_samp,N_b_samp_col]  model_matrix_b_samp; // all samples design matrix
  matrix[N_obs_mb_samp_small,N_b_samp_col]  model_matrix_b_samp_small; // all samples with replicates collapsed

  vector[N_obs_mb_samp]  model_vector_a_samp; // Npcr cycles for each sample, replicates included
  vector[N_obs_mb_samp_small]  model_vector_a_samp_small; // Npcr cycles without replicates

  // Design matrices: mock community samples
  vector[N_obs_mock]  model_vector_a_mock;

  // Priors
  real alpha_prior[2]; // Parameters of normal distribution for prior on alphas
  real tau_prior[2]; // Parameters of gamma distribution for prior on tau (observation precision)
  
  // END DATA FOR METABARCODING
  
  // DATA FOR LINKING QM AND QPCR
  int N_mb_link; //How many qpCR samples have a match in a MB sample
  int mb_link_idx[N_mb_link]; // index linking qpcr samples to mb samples
  int mb_link_sp_idx; // the index for the species linking QM to qPCR (usually hake)

}


transformed data {
  
  vector[NstdSamples] known_conc; // log known concentration of qPCR standards
  matrix[N_obs_mb_samp,N_b_samp_col+1] model_matrix_samp; // QM design matrix (samples by species), where the last column denotes Npcr cycles
  
  known_conc = log(known_concentration);
    
  model_matrix_samp = append_col(model_matrix_b_samp, model_vector_a_samp); //extend design matrix to include Npcr cycles as the last column
  
  // The model_matrix will be multiplied by the relative abundances (log(conc) relative to qpcr reference) and efficiencies (alphas) to obtain copy numbers per species
}

parameters {
  // for QM part
  //real<lower=0> tau_base; // single overdispersion sd for multinomial.
  vector<lower=0>[N_species-1] tau; // single overdispersion sd for multinomial.
  vector[N_species-1] alpha_raw;
  vector[N_obs_mb_samp] eta_samp_raw[N_species-1]; //overdispersion
  vector[N_obs_mock] eta_mock_raw[N_species-1]; //overdispersion
  
  // for qPCR part
  vector<lower=0>[Nplates] beta_std_curve_0; // intercept of standard curve
  vector<upper=0>[Nplates] beta_std_curve_1; // slope of standard curve
  real<upper=0> gamma_0; //intercept to scale variance of standard curves w the mean
  vector<upper=0>[Nplates]  gamma_1; //slopes to scale variance of standards curves w the mean
  real phi_0;
  real<lower=0> phi_1;
  vector[NSamples_qpcr-NstdSamples] envir_concentration; // DNA concentration in unknown samples
  
  //for linking 
  matrix[N_obs_mb_samp,N_species-1] log_B_raw; // estimated true copy numbers by sample, minus the qPCR link species (hake)

}

transformed parameters {
  // for qPCR part
  vector[NSamples_qpcr] Concentration; //this is the concentration of standards and field samples (unknowns) combined
  vector[NSamples_qpcr] Ct; // estimated Ct for all qPCR samples
  vector[NSamples_qpcr] sigma; // SD of Ct values
  vector[NSamples_qpcr] theta; // Bernoulli param, probability of amplification
  
  // for QM part
  vector[N_species] alpha; // vector of coefficients (log-efficiencies relative to reference taxon)
  //vector[N_obs_mb_samp] eta_samp[N_species]; // overdispersion coefficients
  vector[N_obs_mock] eta_mock[N_species]; // overdispersion coefficients
  matrix[N_obs_mb_samp,N_species] mu_samp; // estimates of read counts, in log space
  matrix[N_obs_mock,N_species] mu_mock; // estimates of read counts, in log space
  
  // for linking
  matrix[N_obs_mb_samp,N_species] log_B; // estimated true copy numbers by sample, including the link species
  
 // local variables declaration
  matrix[N_obs_mb_samp,N_species] logit_val_samp;
  matrix[N_obs_mock,N_species] logit_val_mock;
  matrix[N_species,N_obs_mb_samp] prob_samp_t;
  matrix[N_species,N_obs_mock] prob_mock_t;
  
  // qPCR standard curves
  //slot knowns and unknowns into a common vector, where unknowns are treated as params and knowns are treated as data 
  
  // concentration_stds
  // concentration_unks
  Concentration[std_idx] = known_conc; //log scale
  Concentration[unkn_idx] = envir_concentration;
  
  for(i in 1:Nobs_qpcr){
    Ct[plateSample_idx[i]] = beta_std_curve_0[plate_idx[i]] + 
                              beta_std_curve_1[plate_idx[i]] * Concentration[plateSample_idx[i]];

    sigma[plateSample_idx[i]] = exp(gamma_0 + gamma_1[plate_idx[i]]*Concentration[plateSample_idx[i]]);
                                  
    theta[plateSample_idx[i]] = inv_logit(phi_0 + phi_1*exp(Concentration[plateSample_idx[i]]));
                                  
  }
  
  // Link to QM
  for(i in 1:N_species){
    for(j in 1:N_obs_mb_samp){
      if(i==mb_link_sp_idx){ // if index is equal to link species, fill in qpcr estimate
        log_B[j,i] = log(Concentration[plateSample_idx[mb_link_idx[j]]]); 
      }else if(i < mb_link_sp_idx){ //if index is less than link sp. index, fill from log_B_raw
        log_B[j,i] = log_B_raw[j,i];
      }else{ //finally, if index is greater than link sp. index, fill from log_B_raw but shifted by one because of missing column for link species
        log_B[j,i] = log_B_raw[j,i-1];
      }
    }
  };
  
  // QM MODEL PIECES
  
  // Fixed effects components
  alpha[1:(N_species-1)] = alpha_prior[1] + alpha_raw * alpha_prior[2]; 
        // non-centered param beta ~ normal(alpha_prior[1], alpha_prior[2])
  alpha[N_species] = 0; // final species is zero (reference species)

  //tau = rep_vector(tau_base,N_species-1);
  eta_mock[N_species] = rep_vector(0.0,N_obs_mock); // final species is zero (reference species)
  //eta_samp[N_species] = rep_vector(0.0,N_obs_mb_samp); // final species is zero (reference species)
  // random effects vector of vectors
  for (l in 1:(N_species-1)) {
    eta_mock[l] = eta_mock_raw[l] * tau[l] ; // non-centered param eta_mock ~ normal(0,tau)
    //eta_samp[l] = eta_samp_raw[l] * tau[l] ; // non-centered param eta_samp ~ normal(0,tau)
  }

// from qPCR estimates, alphas, and etas we can calculate sample-specific mu
  for (n in 1:N_species) {
    logit_val_samp[,n] = model_matrix_samp * append_row((log_B[,n] - log_B[,N_species]),alpha[n]); 
                            //+eta_samp[n];
    logit_val_mock[,n] = alr_mock_true_prop[,n] + 
                              model_vector_a_mock * alpha[n] + 
                              eta_mock[n];
  }
  for(m in 1:N_obs_mb_samp){
    prob_samp_t[,m] = softmax(transpose(logit_val_samp[m,])); // proportion of each taxon in field samples
  }
  for(m in 1:N_obs_mock){
    prob_mock_t[,m] = softmax(transpose(logit_val_mock[m,])); // proportion of each taxon in mocks
  }
  
  for (n in 1:N_species) {
    mu_samp[,n] = transpose(prob_samp_t)[,n] ; 
    mu_mock[,n] = transpose(prob_mock_t)[,n] ;
  // if(n==1){print("log_prob 1 ",log_prob);}
  }

    // print("MU_SAMP_row",mu_samp[1,]);
    // print("SUM_MU_SAMP",sum(mu_samp[1,]));
    // 
    // print("MU_SAMP_col",mu_samp[,1]);
    // print("SUM_MU_SAMP_col",sum(mu_samp[,1]));
}

model {
  
// qPCR part
  for(i in 1:Nobs_qpcr){
     z[i]   ~ bernoulli(theta[plateSample_idx[i]]);
      // z[i]   ~ bernoulli( inv_logit(theta[plateSample_idx[i]]) ) ;
    }
    
    for(i in 1:Nobs_qpcr){
      if (z[i]==1){ //if Ct observed, then compute likelihood
        y[i] ~ normal(Ct[plateSample_idx[i]], sigma[plateSample_idx[i]]);   
      }
    }

  //beta standard curve params
  beta_std_curve_0 ~ normal(stdCurvePrior_intercept[1], stdCurvePrior_intercept[2]);
  beta_std_curve_1 ~ normal(stdCurvePrior_slope[1], stdCurvePrior_slope[2]);
  
  //gamma params for scaling variance on the standards
  gamma_1 ~ normal(0,5);
  gamma_0 ~ normal(-2,1);
  
  for(i in 1:(N_species-1)){
   log_B_raw[,i] ~ normal(0,2);
  }
  
  envir_concentration ~ normal(0, 2); //log10 scale

  phi_0 ~ std_normal(); //need to fix this
  phi_1 ~ normal(5, 2);
  
  // QM part
  for(i in 1:N_obs_mb_samp){ 
    sample_data[i,] ~  multinomial(transpose(mu_samp[i,])); // Multinomial sampling of mu (proportions in field samples)
  }
  for(i in 1:N_obs_mock){
    mock_data[i,]   ~  multinomial(transpose(mu_mock[i,])); // Multinomial sampling of mu (proportions in mocks)
  }
  // Priors
  for(i in 1:(N_species-1)){
    eta_samp_raw[i] ~ std_normal(); // N(0,tau)
    eta_mock_raw[i] ~ std_normal(); // N(0,tau)
  }
  alpha_raw ~ std_normal(); // prior of normal(alpha_prior[1],alpha_prior[2]);
  tau ~ gamma(tau_prior[1],tau_prior[2]); //
}

