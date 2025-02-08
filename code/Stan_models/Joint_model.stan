
data { 
  
  // DATA FOR qPCR PART OF THE MODEL
  
  int Nplates; // number of PCR plates
  int Nobs_qpcr; // number of field observations for qPCR
  int NSamples_qpcr; //number of unique biological samples, overall
  int NstdSamples; //number of unique biological samples with known concentrations (standards)
  array[Nobs_qpcr] int plate_idx; //index denoting which PCR plate each field sample is on
  array[NstdSamples] int std_plate_idx;//index denoting which PCR plate each standard sample is on
  real beta_std_curve_0_offset;
 
 
  vector[Nobs_qpcr] y_unk; //Ct for field observations
  array[Nobs_qpcr] int z_unk; //indicator for field obs; z = 1 if a Ct was observed, z = 0 otherwise
  vector[NstdSamples] y_std; //Ct for standards
  array[NstdSamples] int z_std; //indicator for standards; z = 1 if a Ct was observed, z = 0 otherwise
  vector[NstdSamples] known_concentration; //known concentration (copies/vol) in standards
  
  array[2] real stdCurvePrior_intercept; // prior on the intercept of the std curves
  array[2] real stdCurvePrior_slope; // prior on the slope of the std curves

  //Covariates and offsets
  vector[Nobs_qpcr] X_offset_tot; //log dilution and log volume offsets together in one vector.

  int N_station_depth;
  matrix[NSamples_qpcr,N_station_depth] X_station_depth_tube;// covariate design matrices
  matrix[Nobs_qpcr,N_station_depth] X_station_depth_obs; //covariate design matrices
  
  int N_bio_rep_RE;
  int N_bio_rep_param;
  int N_bio_rep_idx;
  array[N_bio_rep_idx] int bio_rep_idx; //index of biological replicates
  matrix[NSamples_qpcr,N_bio_rep_RE] X_bio_rep_tube;// covariate design matrices for unique samples
  matrix[Nobs_qpcr,N_bio_rep_RE] X_bio_rep_obs;// covariate design matrices for observations

  vector[Nobs_qpcr] wash_idx;//design matrix for wash effect

  // END DATA FOR qPCR
  
  // DATA FOR METABARCODING PART OF THE MODEL
  
  int N_species; // Number of species in data
  int N_obs_mb_samp;  // Number of observed samples, also the number of groups for qPCR samps to link to MB samps
  int N_obs_mock; // Number of observed mock samples



  // Observed data of community matrices
  array[N_obs_mb_samp, N_species] int sample_data;
  // Observed data of mock community matrices

  // True proportions for mock community in log ratios

  real wash_effect; //estimate of the EtOH wash effect  
  // Read in fixed alphas from mock community.
  vector[N_species] alpha_fix;
    
  // Design matrices: field samples

  vector[N_obs_mb_samp]  model_vector_a_samp; // Npcr cycles for each sample, replicates included
  
  // Design matrices: mock community samples
  

  // Identify a reference species for each observation (most abundant species in each sample)
  array[N_obs_mb_samp] int ref_sp_idx;

  // END DATA FOR METABARCODING
  
  // DATA FOR LINKING QM AND QPCR
  int qpcr_mb_link_sp_idx; // the index for the species linking QM to qPCR (usually hake)
  array[N_obs_mb_samp] int tube_link_idx; //index linking observations to unique biological samples
  real log_D_mu; //prior on mean for log_D_raw, where log_D = log_D_mu + log_D_raw*log_D_scale
  real log_D_scale; //prior on variance param for log_D_raw, where log_D = log_D_mu + log_D_raw*log_D_scale
}

transformed data {
  
  vector[NstdSamples] log_known_conc; // log known concentration of qPCR standards
  vector[N_species] alpha = alpha_fix;
  log_known_conc = log(known_concentration);
}

parameters {
  
  // for qPCR part
  real mean_hake; //global mean hake concentration
  
  vector[Nplates] beta_std_curve_0; // intercept of standard curve
  vector[Nplates] beta_std_curve_1; // slope of standard curve
  real gamma_0; //intercept to scale variance of standard curves w the mean
  real<upper=0> gamma_1; //slopes to scale variance of standards curves w the mean
  real logit_phi_0; // logit scale reduction to the pres-abs model. (phi_0 must be less than 1,greater than 0) 
  
  vector[N_station_depth] log_D_station_depth; // log DNA concentration in field samples in each tube
  real<lower=0> log_D_sigma; //variance on log_D_station_depth
  vector[N_bio_rep_param] bio_rep_param; // log DNA concentration in field samples
    
  real<lower=0> tau_bio_rep; //random effect between biological replicates.
  
  //for linking 
  matrix[N_obs_mb_samp,(N_species-1)] log_D_raw; // estimated true DNA concentration by sample (centered)

}

transformed parameters {
  // for qPCR part
  vector[Nobs_qpcr] Ct; // estimated Ct for all unknown qPCR samples
  vector[Nobs_qpcr] unk_conc_qpcr; //log DNA concentration in field samples observed in qPCR after adjusting for covariates and offsets
  vector[NSamples_qpcr] log_D_station_depth_tube; //log DNA concentration in field samples (tubes)
  vector[NstdSamples] Ct_std; //estimated Ct for standards
  vector[NstdSamples] sigma_std; // SD of Ct values, standards
  vector[NstdSamples] logit_theta_std; // Bernoulli param, probability of amplification, standards
  vector[Nobs_qpcr] sigma_samp; //SD of Ct values, field samples
  vector[Nobs_qpcr] logit_theta_samp; //Probability of amplification, field samples
  vector[N_bio_rep_RE] bio_rep_RE; // log DNA concentration in field samples

  real phi_0 ; // Bound 

  // for QM part
  matrix[N_obs_mb_samp,N_species] logit_val_samp; //species proportions in metabarcoding, logit
  matrix[N_obs_mb_samp,N_species] prop_samp; // proportion of each taxon in field samples= softmax(transpose(logit_val_samp[m,]));
  
  // for linking
  matrix[N_obs_mb_samp,N_species] log_D; // estimated true copy numbers by sample, including the link species
 
  // qPCR standard curves (vectorized)
  Ct_std = beta_std_curve_0_offset+beta_std_curve_0[std_plate_idx] +
                              beta_std_curve_1[std_plate_idx] .* log_known_conc;
  sigma_std = exp(gamma_0 + gamma_1 .* log_known_conc);
  
  
  phi_0 = inv_logit(logit_phi_0) ;
  {// local variable
  vector[NstdSamples] p_tmp; // log known concentration of qPCR standards
    p_tmp = - 2 * known_concentration * phi_0;
    logit_theta_std = log1m_exp(p_tmp) - p_tmp ;
  }
    
  // qPCR unknowns 
    {// locals for making sum-to-0 random effects.
      int count_tot;
      int count_par;
      real bio_rep_sum;
    // random effect of biological replicate 
    // This does depend on the stations and tubes being in order from small to large.
    count_tot = 0;
    count_par = 0;
    for(j in 1:N_bio_rep_idx){
      bio_rep_sum = 0 ;   
      for(k in 1:bio_rep_idx[j]){
        count_tot = count_tot + 1;
        if(k < bio_rep_idx[j]){
          count_par = count_par + 1;
          bio_rep_RE[count_tot] = bio_rep_param[count_par] * tau_bio_rep ;
          bio_rep_sum = bio_rep_sum + bio_rep_RE[count_tot];
        }else if(bio_rep_idx[j]==1){
          bio_rep_RE[count_tot] = 0 ;
        }else{
          bio_rep_RE[count_tot] = -bio_rep_sum;
        }
          } // end k loop
        } // end j loop
      } // end local variables.

  /// THIS IS THE LATENT STATE THAT WILL BE NEEDED TO CONNECT TO THE MB DATA
  log_D_station_depth_tube = mean_hake + X_station_depth_tube * log_D_station_depth +
                            X_bio_rep_tube * bio_rep_RE ;

  /// THIS IS THE LATENT STATE CONNECTS TO THE QPCR OBSERVATIONS
  unk_conc_qpcr = mean_hake + X_station_depth_obs * log_D_station_depth + 
                      X_bio_rep_obs * bio_rep_RE +
                      wash_idx * wash_effect +
                      X_offset_tot ;

  // Vectorized predictions
  Ct = (beta_std_curve_0_offset + beta_std_curve_0[plate_idx]) + beta_std_curve_1[plate_idx].*unk_conc_qpcr;
  sigma_samp = exp(gamma_0 + gamma_1 .* unk_conc_qpcr );
  
  { // local variable}
    vector[Nobs_qpcr] p_tmp; // log known concentration of qPCR standards
  
    p_tmp = - 2 * exp(unk_conc_qpcr) * phi_0;
    logit_theta_samp = log1m_exp(p_tmp) - p_tmp ;
  }
  
  //Link to QM
  for(i in 1:N_species){
    for(j in 1:N_obs_mb_samp){
      if(i==qpcr_mb_link_sp_idx){ // if index is equal to link species (hake), fill in qpcr estimate
        log_D[j,i] = log_D_station_depth_tube[tube_link_idx[j]];
      }else{ // otherwise, fill from log_D_raw
        if(i<qpcr_mb_link_sp_idx){
          log_D[j,i] = log_D_mu+log_D_raw[j,i]*log_D_scale;
        }else{
          log_D[j,i] = log_D_mu+log_D_raw[j,(i-1)]*log_D_scale;
        }
      }
    }
  }

  // QM MODEL PIECES
  
  // Make a vector for the reference species D and for alpha
 {// local variables for making reference species vectors
      vector[N_obs_mb_samp] log_D_ref;
      vector[N_obs_mb_samp] alpha_ref;
      
  for(i in 1:N_obs_mb_samp){
     log_D_ref[i] = log_D[i,ref_sp_idx[i]];
     alpha_ref[i] = alpha[ref_sp_idx[i]];
  }
  
  for (n in 1:N_species) {
    logit_val_samp[,n] = (log_D[,n] - log_D_ref) + model_vector_a_samp.*(alpha[n] - alpha_ref);
  //  logit_val_mock[,n] = alr_mock_true_prop[,n] + model_vector_a_mock .* alpha[n]; //+eta_mock[n]
  }
 }
 
  for(m in 1:N_obs_mb_samp){
    prop_samp[m,] = to_row_vector(softmax(to_vector(logit_val_samp[m,]))); // proportion of each taxon in field samples
  }
}

model{
  
  // qPCR part
  z_std ~ bernoulli_logit(logit_theta_std); 
  
  for(i in 1:NstdSamples){
    if(z_std[i]==1){ //if Ct observed, then compute likelihood
      y_std[i] ~ normal(Ct_std[i],sigma_std[i]);
    }
  }
    
  z_unk   ~ bernoulli_logit(logit_theta_samp);
  
  for(i in 1:Nobs_qpcr){
     if (z_unk[i]==1){ //if Ct observed, then compute likelihood
        y_unk[i] ~ normal(Ct[i], sigma_samp[i]);   
      }
    }

  //beta standard curve params
  beta_std_curve_0 ~ normal(stdCurvePrior_intercept[1]-beta_std_curve_0_offset, stdCurvePrior_intercept[2]);
  beta_std_curve_1 ~ normal(stdCurvePrior_slope[1], stdCurvePrior_slope[2]);
  
  //gamma params for scaling variance on the standards
  gamma_1 ~ std_normal();
  gamma_0 ~ std_normal();
  
  bio_rep_param ~ std_normal(); 
  tau_bio_rep ~ normal(0,0.2);
  
  for(i in 1:(N_species-1)){
    // ONLY set a prior for the species that ARE NOT the qPCR link species (hake)
    // The values for the link species will come from the qPCR part of the joint model
    log_D_raw[,i] ~ std_normal();
  }
  
  logit_phi_0 ~ normal(3,1); //assuming Poisson from bottles to replicates this should be large >3
  
  mean_hake ~normal(2,2); //global mean hake concentration
  log_D_sigma ~ normal(0,3); //variance on log_D_station_depth
  log_D_station_depth ~ normal(0,log_D_sigma); //log scale
  
  // if you're only using the Dirichlet for the mocks...
  for(i in 1:N_obs_mb_samp){
    sample_data[i,] ~  multinomial_logit(to_vector(logit_val_samp[i,])); // Multinomial sampling of mu (proportions in field samples)
    //sample_data[i,] ~  dirichlet_multinomial(to_vector(prop_samp[i,])*dm_alpha0_mock) ;
  }
}

