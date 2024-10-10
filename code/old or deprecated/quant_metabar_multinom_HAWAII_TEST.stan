data { // 
  // LATENT STATE DATA
    int N_species_main; // number of unique species included in any dataset
    int N_site_main;    // number of unique sites included in any dataset
    
  // Size Data
    // This is a padded out data frame for predicted mean size of each species
    // in each site.  It is read in as data and not modified... Can't have NA or 0s
    matrix[N_site_main,N_species_main] log_weight_coef; //
  
  // Covariate Data
   int N_strata ;
   matrix[N_site_main,N_strata] model_matrix_strata ; 
  
   vector[N_site_main] Xs_1 ; 
   int N_Zs_1 ;
   matrix[N_site_main,N_Zs_1] Zs_1 ;     
  
  // COUNT DATA and HELPER FILES
    int N_species_count; // number of unique species included in count dataset
    int N_site_count;   // number of unique sites included in count dataset
    int N_omega_offset; //
  
    array[N_site_count,N_species_count] int Y_count; // matrix for count data
  
    array[N_species_count] int omega_count_idx ; // vector of 0s and index for omega_offsets 
  
    // Mapping large matrix of sites and species to those observed in count data
    array[N_site_count] int site_count_idx ; // This is an index including the rows 
                                        // from the main that are observed in mifish

    matrix[N_species_main,N_species_count] M_main_to_count; // 0 or 1, 1 means map species include the data in the count data
      // post multiple X_matrix by M_main_to_count to get Species in Count observations.
      // produces matrix with rows I_count by columns J_count 
   
  // METABARCODING DATA and HELPER FILES
    // WHICH SITES YOU OBSERVE Mifish AT -
    int N_species_mifish; // Number of species in data
    int N_site_mifish; // Number of sites in data
    int N_obs_mifish ;  // Number of observed samples 
    int N_obs_mifish_mock ; // Number of observed mock samples

      int N_bio_rep_col ; // Number of random effects columns in design matrix for biological replicates
      int N_bio_rep_param ; // number of unique parameters to estimate for biological replicates
      int max_bio_rep ; // maximum number of biological replicates at any site
      array[N_site_mifish] int max_bio_rep_idx ; // largest number of bio rep at each site.
      array[N_bio_rep_param,2] int bio_rep_idx ; // Index for biological replicate parameters

      matrix[N_obs_mifish,N_bio_rep_col ] model_matrix_bio_rep ; // Matrix used for random effects 

    array[N_obs_mifish, N_species_mifish] int Z_mifish ;
    array[N_species_mifish, N_obs_mifish] int Z_mifish_t ;
    array[N_obs_mifish_mock, N_species_mifish] int Z_mifish_mock ; 
    array[N_species_mifish,N_obs_mifish_mock ] int Z_mifish_mock_t ;
   
    // True proportions for mock community
    matrix[N_obs_mifish_mock,N_species_mifish] alr_mock_true_prop ;
   
   // Design matrices: field samples
    matrix[N_obs_mifish,N_site_mifish] model_matrix_mifish_site; //
   
    vector[N_obs_mifish]  model_vector_a_mifish;
    array[N_site_mifish] int site_mifish_idx ; // This is an index including the rows 
                                        // from the main that are observed in mifish

   // Design matrices: mock community samples
    vector[N_obs_mifish_mock] model_vector_a_mifish_mock;
    
    // Mapping matrix, all species to mifish.
    matrix[N_species_main,N_species_mifish] M_main_to_mifish;

  // PRIORS
  // ADD THESE IN AS DATA.  CURRENTLY HARD CODED IN MODEL BLOCK
  
  // Count Priors
    /// ADD.
  
  // Metabarcoding Priors
  array[2] real alpha_prior; // Parameters of normal distribution for prior on alphas
  array[2] real tau_prior; // Parameters of gamma distribution for prior on tau (observation precision)
}

transformed data {
  //vector[N_species_count] omega_count_base;
  matrix[N_obs_mifish,N_site_mifish +1] model_matrix_mifish;
  
  // FOR FIXING ALPHA AT 0
  //vector[N_species_mifish]  alpha ;
  //alpha = rep_vector(0,N_species_mifish);
  
  //int Z_mifish_t[N_species_mifish,N_obs_mifish] ;

  //omega_count_base = rep_vector(0,N_species_count);
  model_matrix_mifish = append_col(model_matrix_mifish_site, model_vector_a_mifish);
  //Z_mifish_t = transpose(Z_mifish) ;
}

parameters {
  //Negative bionomial overdispersion
  real<lower=0> phi;
  // real<lower=0> phi_alpha;
  // real<lower=0> phi_beta;
  // Coefficients for each strata
  array[N_species_main] vector[N_strata] beta_strata;
  // Coefficients for smoothes.
  vector[N_species_main] beta_Xs_1;
  array[N_species_main] vector[N_Zs_1] beta_Zs_1;
  
  vector<lower=0>[N_species_main] sds_1;  // standard deviations of spline coefficients

 // Random effect for each species in each location.
  // vector<lower=0>[N_species_main] sigma_log_X;
  // matrix[N_site_main,N_species_main] epsilon_X_raw ;
  
  // VISUAL COUNT PARAMETERS
  // real<lower=0> sigma_count[N_species_count] ; // Overdispersion parameter for each species.
  vector[N_omega_offset]  omega_logit; // catchability offsets for each species in logit-space
  // matrix[N_site_count,N_species_count] epsilon_count_raw ;

  // MIFISH PARAMETERS
  array[N_species_mifish] real<lower=0> tau_B_mifish; // random effect accounting for the translating counts to biomass.
  matrix[N_site_mifish,N_species_mifish] gamma_B_mifish_raw; //realizations of tau_B_mifish
  
  // real tau_0; // intercept for overdispersion as a function of biomass for multinomial.
  // real tau_1; // slope for overdispersion as a function of biomass for multinomial.
  // real<lower=0>tau_mock; // overdispersion among mock community samples.
   
  vector[N_species_mifish-1] alpha_raw;
  array[N_species_mifish-1] real<lower=0> tau_bio_rep ;
  //real<lower=0> tau_bio_rep ;
  
  array[N_species_mifish-1] vector[N_bio_rep_param] kappa_mifish ; // random effect for bio replicates in mifish
  //array[N_species_mifish-1] vector[N_bio_rep_param] kappa_mifish_raw ; // random effect for bio replicates in mifish
  // matrix[N_obs_mifish,N_species_mifish-1] eta_mifish_raw;
  // matrix[N_obs_mifish_mock,N_species_mifish-1] eta_mifish_mock_raw;
}

transformed parameters {
  
  // log abundance (fish density) at each site.
  matrix[N_site_main,N_species_main] log_X; 
  matrix[N_site_main,N_species_main] log_B; 

  // Actual spline coefficients.(vector for each species of vectors)
  array[N_species_main] vector[N_Zs_1] s_1;

  // COUNT DATA TRANSFORMED PARAMETERS
  matrix[N_site_count,N_species_count] log_X_count ;
  //matrix[N_site_count,N_species_count] lambda  ;
  //matrix[N_site_count,N_species_count] epsilon_count ;
  vector[N_species_count]  omega_count_long;

  // METABARCODING TRANSFORMED PARAMETERS
  matrix[N_site_mifish,N_species_mifish] log_B_mifish; // log abundance (fish density) at each site
  vector[N_species_mifish] alpha; // coefficients
  
  matrix[N_site_mifish,N_species_mifish] gamma_B_mifish; //realizations of tau_B_mifish
  
  //array[N_species_mifish] vector[N_bio_rep_param] kappa_mifish ; // random effect for bio replicates in mifish
  
  array[N_species_mifish] matrix[N_site_mifish,  max_bio_rep] kappa_mifish_mat ; 
  
  
  // matrix[N_obs_mifish,N_species_mifish] eta_mifish; // overdispersion coefficients
  // matrix[N_obs_mifish_mock,N_species_mifish] eta_mifish_mock; // overdispersion coefficients
  matrix[N_species_mifish,N_obs_mifish] mu_samp;   // estimates, in log space
  matrix[N_species_mifish,N_obs_mifish_mock] mu_mock;   // estimates, in log space
  
  // Compute actual spline coefficients
  for(i in 1: N_species_main){
    s_1[i] = sds_1[i] * beta_Zs_1[i];
  }

  // Vectorized true log-abundance as a function of strata and other covariates
  for(i in 1:N_species_main){
    log_X[,i] = model_matrix_strata * beta_strata[i] + // Strata effects
                      Xs_1 * beta_Xs_1[i] + // first smooth
                      Zs_1 * s_1[i] ; //+ epsilon_X_raw[,i]*sigma_log_X[i] ; 
  }
  
  // Biomass latent variable.
   log_B = log_X + log_weight_coef;
  
  // COUNT MODEL SECTION
  /// make sure the 
    // { // and trims to include only sites with mifish and species groups observable by count.
    // matrix[N_site_count,N_species_main] log_X_temp ;
    //   for(i in 1:N_site_count){
    //     log_X_temp[i,] =  ;
    //   }
    //    
    // }
    log_X_count = log_X[site_count_idx,] * M_main_to_count; 
  // Catchability offsets
  for(i in 1:N_species_count){
    if(omega_count_idx[i]>0){
      omega_count_long[i] =  -log(1+exp(- omega_logit[omega_count_idx[i]]));
    }else{
      omega_count_long[i] = 0 ;
    }
  }
  
    // print("log_X ",log_X[,1]);
    
    // Make predictions for counts including each species-specific catchability
    // for( i in 1:N_species_count){
    //  // epsilon_count[,i] = epsilon_count_raw[,i] * sigma_count[i] ;
    //   lambda[,i] = log_X_count[,i] +  omega_count_idx[i] * omega_count[i] ;//+ epsilon_count[,i] ;
    //   // 
    //   // print("omega ",omega_count_idx[i] * omega_count[i]) ;
    //   // print("log_X ",log_X_count[,i]) ;
    //   // print("lambda ",lambda[,i]) ;
    // }

  //MIFISH MODEL SECTION
  // Initiate local variable
    { // trim to include only sites with mifish and species groups observable by count.
    matrix[N_site_mifish,N_species_main] log_B_temp ;
      for(i in 1:N_site_mifish){
        log_B_temp[i,] = log_B[site_mifish_idx[i],]  ;
      }
      log_B_mifish = log_B_temp * M_main_to_mifish ;
    } // end local variable

    for(i in 1:(N_species_mifish)){
      gamma_B_mifish[,i] = gamma_B_mifish_raw[,i] * tau_B_mifish[i]; //
      log_B_mifish [,i] = log_B_mifish [,i] + gamma_B_mifish[,i]; //
    }

  { // local variables delcaration
  matrix[N_obs_mifish,N_species_mifish] v_samp;
  matrix[N_obs_mifish_mock,N_species_mifish] v_mock;
  // matrix[N_species_mifish,N_obs_mifish] prob_samp_t;
  // matrix[N_species_mifish,N_obs_mifish_mock] prob_mock_t;
  
    // Fixed effects components.
    alpha[1:(N_species_mifish-1)] = alpha_prior[1] + alpha_raw * alpha_prior[2]; 
          // non-centered param alpha ~ normal(alpha_prior[1], alpha_prior[2])
    alpha[N_species_mifish] = 0; // final species is zero (reference species)
    
    // // random effect of biological replicate 
    // for(i in 1:N_species_mifish-1){
    //   kappa_mifish[i] = kappa_mifish_raw[i] * tau_bio_rep[i];
    // }          // non-centered param beta ~ normal(alpha_prior[1], alpha_prior[2])
    // kappa_mifish[N_species_mifish] = rep_vector(0.0,N_bio_rep_param); // final species is zero (reference species)

    // Create Kappa Mifish Matrix to ensure sum to zero constraint on random effects.
    // Make last value needed in each row equal to the negative of the sum of the other components.
    for(i in 1:N_species_mifish-1){
      kappa_mifish_mat[i] = rep_matrix(0,N_site_mifish, max_bio_rep);
      for(j in 1:N_bio_rep_param){
        kappa_mifish_mat[i,bio_rep_idx[j,1],bio_rep_idx[j,2]] = kappa_mifish[i,j];
      }
      for(k in 1:N_site_mifish){ // This calculates the last random effect.
        kappa_mifish_mat[i,k,max_bio_rep_idx[k]] = -sum(kappa_mifish_mat[i,k,:(max_bio_rep_idx[k]-1)]) ;
      }
    }          // non-centered param beta ~ normal(alpha_prior[1], alpha_prior[2])
    kappa_mifish_mat[N_species_mifish] = rep_matrix(0,N_site_mifish, max_bio_rep);
    
    // add over dispersion have a normal distribution.
    // eta_mifish_mock[N_species_mifish] = rep_vector(0.0,N_obs_mifish_mock); // final species is zero (reference species)
    // eta_mifish[,N_species_mifish] = rep_vector(0.0,N_obs_mifish); // final species is zero (reference species)
    //     
    // for (n in 1:(N_species_mifish-1)) {
    //     // eta_mifish_mock[n] = eta_mifish_mock_raw[n] * tau_mock ;
    //     eta_mifish[,n] = eta_mifish_raw[,n] .* exp(tau_0 + tau_1 * model_matrix_mifish_site * (log_B_mifish[,n] - log_B_mifish[,N_species_mifish])) ;
    // }
  
  // from betas, alphas, and etas we can calculate sample-specific mu
    for (n in 1:N_species_mifish) {
      v_samp[,n] =  model_matrix_mifish * append_row(log_B_mifish[,n] - log_B_mifish[,N_species_mifish],alpha[n]) +
                              model_matrix_bio_rep * to_vector(kappa_mifish_mat[n]) ; 
                              //eta_mifish[,n];
      v_mock[,n] = alr_mock_true_prop[,n] + 
                                 model_vector_a_mifish_mock * alpha[n] ;  
      //                           eta_mock[n];
    }

    mu_samp = transpose(v_samp);
    mu_mock = transpose(v_mock);
    for(m in 1:N_obs_mifish){
      mu_samp[,m] = softmax(mu_samp[,m]);
      // Enforce small amount in each species in each sample to avoid 0
      // errors
      
    }
    // print("mu_samp[,1] ", mu_samp[,1]);
    // print(".1  ");
    // print("mu_samp[,3] ", mu_samp[,3]);
    // print(".2  ");

    for(m in 1:N_obs_mifish_mock){
       mu_mock[,m] = softmax(mu_mock[,m]);
    }
    
    // print("alpha ",alpha);
    // print("alr_mock_true_prop",alr_mock_true_prop);
    // print("v mock",v_mock);
    // print("mu mock",mu_mock);
    // print("mu row ",sum(mu_samp[1,]));
    // print("mu col ",mu_samp[,1]);

    // print(mu_samp[1,]);

    // print("MU_SAMP_row",mu_samp[1,]);
    // print("SUM_MU_SAMP",sum(mu_samp[1,]));
    // 
    // print("MU_SAMP_col",mu_samp[,1]);
    // print("SUM_MU_SAMP_col",sum(mu_samp[,1]));

  } // end local variables declaration
  // print("mu 1 ",mu[1]);
  // print(N_pcr_samp*log(1+alpha[1]));
}

model {
  // Likelihoods for Count data:
  for(i in 1:N_species_count){
    Y_count[,i] ~ neg_binomial_2( exp(log_X_count[,i] +  omega_count_long[i]),phi ) ;
    // Y_count[,i] ~ poisson_log( log_X_count[,i] +  omega_count_idx[i] * omega_count[i] ) ;
    // print("count ",i) ;
  }

  // Likelihoods and priors for Metabarcoding:
  for(i in 1:N_obs_mifish){
    Z_mifish_t[,i] ~  multinomial(mu_samp[,i]);
    
  }
  // print("mu_samp 1",mu_samp[,1]) ;
  // print("mu_samp 2 ",mu_samp[,2]) ;
  
  for(i in 1:N_obs_mifish_mock){
     Z_mifish_mock_t[,i]   ~  multinomial(mu_mock[,i]);  
  }
  
  // Priors ////////////////////////////////////
  // Latent variable level priors
  // Strata and smooth Priors
  for(i in 1:N_species_main){
    beta_strata[i] ~ normal(0,3) ;
    beta_Xs_1[i] ~  normal(0,3) ;
    beta_Zs_1[i] ~ std_normal();
    target += student_t_lpdf(sds_1 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  }
  
  // sigma_log_X ~ gamma(10,10) ;

  // VISUAL COUNT Priors
  // sigma_count ~ gamma(1,1) ; // Overdispersion parameter for each species.
  omega_logit ~ normal(2, 4) ;// catchability offsets for each species.
  // for(i in 1:N_species_count){
  //   epsilon_count_raw[,i] ~ std_normal() ; // non-centered parameterization
  // }

  // MIFISH parameter Priors
  for(i in 1:(N_species_mifish-1)){
    kappa_mifish[i] ~ normal(0,tau_bio_rep[i]) ;
  }
   
  // for(i in 1:(N_species_mifish-1)){
  // eta_mifish_raw[,i] ~ std_normal(); // N(0,tau)
  //   eta_mifish_mock_raw[,i] ~ std_normal(); // N(0,tau)
  // }
  
  //alpha_raw ~ std_normal(); // prior of normal(alpha_prior[1],alpha_prior[2]);
  tau_B_mifish ~ gamma(2,4) ; // random effect accounting for the translating counts to biomass.
  for(i in 1:(N_species_mifish)){
    gamma_B_mifish_raw[,i] ~ std_normal(); //realizations of tau_B_mifish... accounts for random variation in biomass at a given strata. among strat
  } 

  // tau_mock ~ gamma(tau_prior[1],tau_prior[2]); //
  // tau_0 ~ normal(0,0.1); //
  // tau_1 ~ std_normal(); //
  tau_bio_rep ~ gamma(5,10) ; // standard deviation for each species for biological replicates
  phi ~ gamma(2,2) ;
}

generated quantities {
  // MAKE MU PREDICTIONS.
  // MAKE PROPER omega_count
  // MAKE LAMBDA PREDICTIONS
  matrix[N_site_count,N_species_count] lambda  ;
  matrix[N_site_main,N_species_main] pred_X_smooth  ;
  matrix[N_site_main,N_species_main] pred_X_fixed  ;
  matrix[N_site_main,N_species_main] pred_X_all  ;

  for( i in 1:N_species_count){
    lambda[,i] = log_X_count[,i] +  omega_count_long[i] ;
  }

  for( i in 1:N_species_main){
    pred_X_smooth[,i] = Xs_1 * beta_Xs_1[i] + // first smooth
                      Zs_1 * s_1[i] ; //+ epsilon_X_raw[,i]*sigma_log_X[i] ; 
    pred_X_fixed[,i] = model_matrix_strata * beta_strata[i] ; // Strata effects
    pred_X_all[,i] = pred_X_smooth[,i] + pred_X_fixed[,i];
  }

    
  //matrix[N_obs_samp_small,N_species] int_samp_small; 
  //     // estimates, in proportion space one for each field sample.
  // 
  //   { // local variables delcaration
  //   matrix[N_obs_samp_small,N_species] logit_val_samp_ppd;
  //   matrix[N_species,N_obs_samp_small] prob_samp_t_ppd;
  // 
  //     for(n in 1:N_species) {
  //       logit_val_samp_ppd[,n] = model_matrix_b_samp_small * beta[n] ; //+
  //     }
  //     for(m in 1:N_obs_samp_small){
  //       prob_samp_t_ppd[,m] = softmax(transpose(logit_val_samp_ppd[m,]));
  //     }
  //     for(n in 1:N_species) {
  //       int_samp_small[,n] = transpose(prob_samp_t_ppd)[,n] ;
  //     }
  //   } // End local variable declaration.

}
