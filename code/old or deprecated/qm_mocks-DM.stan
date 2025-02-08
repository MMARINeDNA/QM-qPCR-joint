
data { 
  
  // DATA FOR qPCR PART OF THE MODEL

  // DATA FOR METABARCODING PART OF THE MODEL
  
  int N_species; // Number of species in data
  int N_obs_mock; // Number of observed mock samples

  // Observed data of mock community matrices
  array[N_obs_mock,N_species] int mock_data;  
  // True proportions for mock community in log ratios
  matrix[N_obs_mock,N_species] alr_mock_true_prop ;
    
  // Design matrices: mock community samples
  vector[N_obs_mock]  model_vector_a_mock;

  // Priors
  array[2] real alpha_prior;// Parameters of normal distribution for prior on alphas
  // real dm_alpha0_mock; // if you want a fixed Dirichlet alpha0 value
  // real tau_prior[2]; // Parameters of gamma distribution for prior on tau (observation precision)
}

transformed data {
}

parameters {
  
  // for QM part
  vector[N_species-1] alpha_raw; // log-efficiencies of PCR in MB
  real log_dm_alpha0_mock; //log-scale alpha param for the Dirichlet multinomial, mocks
}

transformed parameters {

  // for QM part
  vector[N_species] alpha; // vector of efficiency coefficients (log-efficiencies relative to reference taxon)
  // vector[N_obs_mb_samp] eta_samp[N_species]; // overdispersion coefficients
  // vector[N_obs_mock] eta_mock[N_species]; // overdispersion coefficients
  real dm_alpha0_mock; // alpha param for the Dirichlet multinomial, mocks
  matrix[N_obs_mock,N_species] logit_val_mock; //species proportions in metabarcoding, logit
  matrix[N_obs_mock,N_species] prop_mock; // proportion of each taxon in field samples

  // QM MODEL PIECES
  // Fixed effects components
  alpha[1:(N_species-1)] = alpha_prior[1] + alpha_raw * alpha_prior[2];
        // non-centered param beta ~ normal(alpha_prior[1], alpha_prior[2])
  alpha[N_species] = 0; // final species is zero (reference species)

  for (n in 1:N_species) {
    logit_val_mock[,n] = alr_mock_true_prop[,n] + model_vector_a_mock .* alpha[n]; //+eta_mock[n]
  }
  
  dm_alpha0_mock = exp(log_dm_alpha0_mock);
  
  for(m in 1:N_obs_mock){
    prop_mock[m,] = to_row_vector(softmax(to_vector(logit_val_mock[m,]))); // proportion of each taxon in mocks
  }  
 
}

model{
  for(i in 1:(N_species-1)){
    alpha_raw[i] ~ std_normal();
  }
  log_dm_alpha0_mock ~ normal(8,4); // prior on log of Dirichlet multinomial alpha0 for mock communities
  
  // QM Likelihoods
  for(i in 1:N_obs_mock){
    mock_data[i,]   ~  dirichlet_multinomial(to_vector(prop_mock[i,])*dm_alpha0_mock); // Multinomial sampling of mu (proportions in mocks)
  }
}

