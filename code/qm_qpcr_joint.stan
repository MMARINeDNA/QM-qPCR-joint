data { // 
  
  // DATA FOR METABARCODING PART OF THE MODEL
  int N_species; // Number of species in data
  int N_obs_samp;  // Number of observed samples 
  int N_obs_samp_small;  // Number of observed samples for individual sites.
  int N_obs_mock; // Number of observed mock samples

  // Observed data of community matrices
  int sample_data[N_obs_samp,N_species];
  // Observed data of mock community matrices
  int mock_data[N_obs_mock,N_species];  
  // True proportions for mock community
  matrix[N_obs_mock,N_species] alr_mock_true_prop ;
 // matrix[N_obs_mock_small,N_species] alr_mock_true_prop_small ;
    
  // Design matrices: field samples
  int N_b_samp_col; //Number of
  matrix[N_obs_samp,N_b_samp_col]  model_matrix_b_samp;
  matrix[N_obs_samp_small,N_b_samp_col]  model_matrix_b_samp_small;

  vector[N_obs_samp]  model_vector_a_samp;
  vector[N_obs_samp_small]  model_vector_a_samp_small;

  // Design matrices: mock community samples
  //int N_b_mock_col;
  //matrix[N_obs_mock,N_b_mock_col]  model_matrix_b_mock;
  vector[N_obs_mock]  model_vector_a_mock;
 //vector[N_obs_mock_small]  model_vector_a_mock_small;

  // Priors
  real alpha_prior[2]; // Parameters of normal distribution for prior on alphas
  real beta_prior[2]; // Parameters of normal distribution for prior on betas
  real tau_prior[2]; // Parameters of gamma distribution for prior on tau (observation precision)
  
  // END DATA FOR METABARCODING
  
  // DATA FOR qPCR PART OF THE MODEL
    
  int Nplates; // number of PCR plates
  int Nobs; // number of field observations for qPCR
  int NSamples_qpcr; //number of unique biol samples, overall
  int NstdSamples; //number of unique biol samples with known concentrations (standards)
  int plate_idx[Nobs]; //index denoting which PCR plate each sample is on
  int std_idx[NstdSamples]; //index relative to NSamples; which ones are the standards?
  int unkn_idx[NSamples_qpcr-NstdSamples]; //index relative to total samples; which ones are the unknown/field samples?
  int plateSample_idx[Nobs]; //index of unique combinations of plate and biological sample
  
  
  vector[Nobs] y; //Ct observations
  int z[Nobs]; //indicator; z = 1 if a Ct was observed, z = 0 otherwise
  vector[NstdSamples] known_concentration;
  
  real stdCurvePrior_intercept[2]; // prior on the intercept of the std curves
  real stdCurvePrior_slope[2]; // prior on the slope of the std curves

  //vector[NSamples-NstdSamples] dilutionFactor;
}


transformed data {
  vector[NstdSamples] known_conc; // log known concentration (qPCR standards)
  matrix[N_obs_samp,N_b_samp_col+1] model_matrix_samp; // QM design matrix (observations by species), where the last column denotes 
  
  known_conc = log10(known_concentration);
    
  model_matrix_samp = append_col(model_matrix_b_samp, model_vector_a_samp);
}