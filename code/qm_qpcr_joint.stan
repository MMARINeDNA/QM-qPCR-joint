data { 
  // DATA FOR qPCR PART OF THE MODEL
  
  int Nplates; // number of PCR plates
  int Nobs_qpcr; // number of field observations for qPCR
  int NSamples_qpcr; //number of unique biological samples, overall
  int NstdSamples; //number of unique biological samples with known concentrations (standards)
  int plate_idx[Nobs_qpcr]; //index denoting which PCR plate each field sample is on
  int std_plate_idx[NstdSamples];//index denoting which PCR plate each standard sample is on
  int qpcr_sample_idx[Nobs_qpcr]; //index linking observations to unique biological samples
 
  vector[Nobs_qpcr] y_unk; //Ct for field observations
  int z_unk[Nobs_qpcr]; //indicator for field obs; z = 1 if a Ct was observed, z = 0 otherwise
  vector[NstdSamples] y_std; //Ct for standards
  int z_std[NstdSamples]; //indicator for standards; z = 1 if a Ct was observed, z = 0 otherwise
  vector[NstdSamples] known_concentration; //known concentration (copies/vol) in standards
  
  real stdCurvePrior_intercept[2]; // prior on the intercept of the std curves
  real stdCurvePrior_slope[2]; // prior on the slope of the std curves

  //Covariates and offsets
  vector[NSamples_qpcr] log_dil; //log dilution
  vector[NSamples_qpcr] wash_idx;//design matrix for wash effect
  real wash_prior[2]; //priors for wash offset ~ N(wash_offset_prior[1],wash_offset_prior[2]) 
  
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
  int N_b_samp_col; // Number of samples
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
  int mb_link_idx[N_mb_link]; // index: which qpcr samples (plateSample_idx) does each MB sample correspond to?
  int mb_link_sp_idx; // the index for the species linking QM to qPCR (usually hake)

}


transformed data {
  vector[NstdSamples] known_conc; // log known concentration of qPCR standards
  matrix[N_obs_mb_samp,N_b_samp_col+1] model_matrix_samp; // QM design matrix (samples by species), where the last column denotes Npcr cycles
  
  known_conc = log(known_concentration);
    
  model_matrix_samp = append_col(model_matrix_b_samp, model_vector_a_samp); //extend MB design matrix to include Npcr cycles as the last column
  
  // The model_matrix will be multiplied by the relative abundances (log(conc) relative to qpcr reference) and efficiencies (alphas) to obtain copy numbers per species
}

parameters {
  // for QM part
  //real<lower=0> tau_base; // single overdispersion sd for multinomial.
  // vector<lower=0>[N_species-1] tau; // single overdispersion sd for multinomial.
  // vector[N_species-1] alpha_raw;
  // //vector[N_obs_mb_samp] eta_samp_raw[N_species-1]; //overdispersion
  // vector[N_obs_mock] eta_mock_raw[N_species-1]; //overdispersion
  
  // for qPCR part
  vector<lower=0>[Nplates] beta_std_curve_0; // intercept of standard curve
  vector<upper=0>[Nplates] beta_std_curve_1; // slope of standard curve
  real gamma_0; //intercept to scale variance of standard curves w the mean
  real<upper=0> gamma_1; //slopes to scale variance of standards curves w the mean
  real phi_0;
  real<lower=0> phi_1;
  real wash_effect; //estimate of the EtOH wash effect
  vector[NSamples_qpcr] unk_conc_raw; // log DNA concentration in field samples
  //for linking 
  // matrix[N_obs_mb_samp,N_species] log_D_raw; // estimated true copy numbers by sample
}

transformed parameters {
  // for qPCR part
  vector[Nobs_qpcr] Ct; // estimated Ct for all unknown qPCR samples
  vector[NSamples_qpcr] unk_conc; //log DNA concentration in field samples after adjusting for covariates and offsets
  vector[NstdSamples] Ct_std; //estimated Ct for standards
  vector[NstdSamples] sigma_std; // SD of Ct values, standards
  vector[NstdSamples] theta_std; // Bernoulli param, probability of amplification, standards
  vector[Nobs_qpcr] sigma_samp; //SD of Ct values, field samples
  vector[Nobs_qpcr] theta_samp; //Probability of amplification, field samples
  
  // for QM part
  // vector[N_species] alpha; // vector of coefficients (log-efficiencies relative to reference taxon)
  // //vector[N_obs_mb_samp] eta_samp[N_species]; // overdispersion coefficients
  // vector[N_obs_mock] eta_mock[N_species]; // overdispersion coefficients
  // matrix[N_obs_mb_samp,N_species] mu_samp; // estimates of read counts, in log space
  // matrix[N_obs_mock,N_species] mu_mock; // estimates of read counts, in log space
  // 
  // for linking
  // matrix[N_obs_mb_samp,N_species] log_D; // estimated true copy numbers by sample, including the link species
  
 
 { // local variables declaration
  matrix[N_obs_mb_samp,N_species] logit_val_samp;
  matrix[N_obs_mock,N_species] logit_val_mock;
  matrix[N_species,N_obs_mb_samp] prob_samp_t;
  matrix[N_species,N_obs_mock] prob_mock_t;
  
  // qPCR standard curves
  for(i in 1:NstdSamples){
    Ct_std[i] = beta_std_curve_0[std_plate_idx[i]] + 
                              beta_std_curve_1[std_plate_idx[i]] * known_conc[i];

    //sigma_std[i] = exp(gamma_0 + gamma_1[std_plate_idx[i]]*known_conc[i]);
    sigma_std[i] = exp(gamma_0 + gamma_1*known_conc[i]);
    
    theta_std[i] = inv_logit(phi_0 + phi_1*exp(known_conc[i]));
  }
  // qPCR unknowns
  
  unk_conc = unk_conc_raw + log_dil + wash_effect * wash_idx;
  
  // for(i in 1:NSamples_qpcr){
  //   unk_conc[i] = unk_conc_raw[i]+log_dil[i]+wash_effect*wash_idx[i];
  // }
  for(i in 1:Nobs_qpcr){
    Ct[i] = beta_std_curve_0[plate_idx[i]]+beta_std_curve_1[plate_idx[i]]*unk_conc[qpcr_sample_idx[i]];
    
    sigma_samp[i] = exp(gamma_0 + gamma_1*(unk_conc[qpcr_sample_idx[i]]));
    
    theta_samp[i] = inv_logit(phi_0+phi_1*exp(unk_conc[qpcr_sample_idx[i]]));
  }
 } // end local variables
  // Link to QM
//   for(i in 1:N_species){
//     for(j in 1:N_obs_mb_samp){
//       if(i==mb_link_sp_idx){ // if index is equal to link species, fill in qpcr estimate
//         log_D[j,i] = unk_conc[qpcr_sample_idx[mb_link_idx[j]]]; 
//       }else{ // otherwise, fill from log_D_raw
//         log_D[j,i] = log_D_raw[j,i];
//       }
//     }
//   };
// //  print("rows",log_D[7,]);
//  print("columns",log_D[,2]);
  
  // QM MODEL PIECES
  
  // Fixed effects components
//   alpha[1:(N_species-1)] = alpha_prior[1] + alpha_raw * alpha_prior[2]; 
//         // non-centered param beta ~ normal(alpha_prior[1], alpha_prior[2])
//   alpha[N_species] = 0; // final species is zero (reference species)
// 
//   //tau = rep_vector(tau_base,N_species-1);
//   eta_mock[N_species] = rep_vector(0.0,N_obs_mock); // final species is zero (reference species)
//   //eta_samp[N_species] = rep_vector(0.0,N_obs_mb_samp); // final species is zero (reference species)
//   // random effects vector of vectors
//   for (l in 1:(N_species-1)) {
//     eta_mock[l] = eta_mock_raw[l] * tau[l] ; // non-centered param eta_mock ~ normal(0,tau)
//     //eta_samp[l] = eta_samp_raw[l] * tau[l] ; // non-centered param eta_samp ~ normal(0,tau)
//   }
// 
// // from qPCR estimates, alphas, and etas we can calculate sample-specific mu
//   for (n in 1:N_species) {
//     logit_val_samp[,n] = model_matrix_samp * append_row((log_D[,n] - log_D[,N_species]),alpha[n]); 
//                             //+eta_samp[n];
//     logit_val_mock[,n] = alr_mock_true_prop[,n] + 
//                               model_vector_a_mock * alpha[n] + 
//                               eta_mock[n];
//   }
//   for(m in 1:N_obs_mb_samp){
//     prob_samp_t[,m] = softmax(transpose(logit_val_samp[m,])); // proportion of each taxon in field samples
//   }
//   for(m in 1:N_obs_mock){
//     prob_mock_t[,m] = softmax(transpose(logit_val_mock[m,])); // proportion of each taxon in mocks
//   }
//   
//   for (n in 1:N_species) {
//     mu_samp[,n] = transpose(prob_samp_t)[,n] ; 
//     mu_mock[,n] = transpose(prob_mock_t)[,n] ;
//   // if(n==1){print("log_prob 1 ",log_prob);}
//   }
// 
//     // print("MU_SAMP_row",mu_samp[1,]);
    // print("SUM_MU_SAMP",sum(mu_samp[1,]));
    // 
    // print("MU_SAMP_col",mu_samp[,1]);
    // print("SUM_MU_SAMP_col",sum(mu_samp[,1]));
}

model {
  
  
  // qPCR part
  for(i in 1:NstdSamples){
    z_std[i] ~ bernoulli(theta_std[i]);
    if(z_std[i]==1){
      y_std[i] ~ normal(Ct_std[i],sigma_std[i]);
    }
  }

  for(i in 1:Nobs_qpcr){
     z_unk[i]   ~ bernoulli(theta_samp[i]);
     if (z_unk[i]==1){ //if Ct observed, then compute likelihood
        y_unk[i] ~ normal(Ct[i], sigma_samp[i]);   
      }
    }

// print("Obs= ",y_unk[1:5]," ct =",Ct[1:5],";sigma = ",sigma_samp[1:5]);

  //beta standard curve params
  beta_std_curve_0 ~ normal(stdCurvePrior_intercept[1], stdCurvePrior_intercept[2]);
  beta_std_curve_1 ~ normal(stdCurvePrior_slope[1], stdCurvePrior_slope[2]);
  
  //gamma params for scaling variance on the standards
  gamma_1 ~ normal(0,1);
  gamma_0 ~ normal(-2,1);
  
  unk_conc ~ normal(0,10); //log scale
  
  // for(i in 1:(N_species)){
  //   // ONLY set a prior for the species that ARE NOT the qPCR link species (hake)
  //   // The values for the link species will come from the qPCR part of the joint model
  //   // (which will use prior information from envir_concentration, above)
  //   if(i!=mb_link_sp_idx){
  //     log_D_raw[,i] ~ normal(0,10); 
  //   }
  // }
  
  phi_0 ~ normal(0.58,0.2); //assuming Poisson from bottles to replicates (pipetting)
  phi_1 ~ normal(5, 2);
  
  wash_effect ~ normal(wash_prior[1],wash_prior[2]) ;
  
  // QM part
  // for(i in 1:N_obs_mb_samp){
  //   sample_data[i,] ~  multinomial(transpose(mu_samp[i,])); // Multinomial sampling of mu (proportions in field samples)
  // }
  // for(i in 1:N_obs_mock){
  //   mock_data[i,]   ~  multinomial(transpose(mu_mock[i,])); // Multinomial sampling of mu (proportions in mocks)
  // }
  // // Priors
  // for(i in 1:(N_species-1)){
  //   // eta_samp_raw[i] ~ std_normal(); // N(0,tau)
  //   eta_mock_raw[i] ~ std_normal(); // N(0,tau)
  // // }
  // alpha_raw ~ std_normal(); // prior of normal(alpha_prior[1],alpha_prior[2]);
  // tau ~ gamma(tau_prior[1],tau_prior[2]); //
  // }
}

