// GLED's MODEL for linking trawl catches

functions {
  real[] average_by_idx(int N, int K, real[] input, int[] index) {
    real output[K];
    int count[K];
    count = rep_array(0, K);
    output = rep_array(0.0, K);
    for (i in 1:N) {
      count[index[i]] = count[index[i]] + 1;
      output[index[i]] = output[index[i]] + input[i];
    }
    for (i in 1:K) {
      // if(count[i]!=0){
        output[i] = output[i] / count[i];
        // }
    }
    return output;
  }
}
data {
  //////////////////////////////////////////// Intigers
  int N_yest_yj_M; // Number of years and station for the entire joint Model
  int N_ye_y_M; // Number of years for the entire joint Model
  //
    int N_sp_i_MC1; // Number of species Model Compartment 1
  int N_sp_i_MC2; // Number of species Model Compartment 2
  int N_spye_iy_MC1; // Number of species * years in Model Compartment 1
  int N_spyest_iyj_MC1; // Number of species * years * stations in Model Compartment 1
  //
    int N_obs_Z_M1; // Number of observation in Model 1
  int N_obs_Z_M2; // Number of observation in Model 2
  int N_obs_W_M3; // Number of observation in Model 3
  int N_obs_W_M4; // Number of observation in Model 4
  int N_obs_Y_M5; // Number of observation in Model 5
  int N_obs_Y_M6; // Number of observation in Model 6
  //////////////////////////////////////////// Data
  //
    int Z_M1[N_obs_Z_M1]; // number of fish caught in Model 1
  real E_M1[N_obs_Z_M1]; // Trawl effort (km2) in Model 1
  //
    int W_M3[N_obs_W_M3]; // Number of positive droplets in Model 3
  int U_M3[N_obs_W_M3]; // Number of total droplets in Model 3
  real S_M3[N_obs_W_M3]; // Standards concentration in Model 3
  int W_M4[N_obs_W_M4]; // Number of positive droplets in Model 4
  int U_M4[N_obs_W_M4]; // Number of total droplets in Model 4
  //
    int Z_M2[N_sp_i_MC2-1,N_obs_Z_M2]; // umber of fish caught in Model 2
  real E_M2[N_obs_Z_M2]; // Trawl effort (km2) in Model 2
  //
    vector[N_sp_i_MC2] alr_M5; // Additive log ratio of mock initial concentration in Model 5
  int Y_M5[N_sp_i_MC2,N_obs_Y_M5]; // Metabarcoding reads in Model 5
  int Y_M6[N_sp_i_MC2,N_obs_Y_M6]; // Metabarcoding reads in Model 6
  int NPCR; //Number of PCR cycles in Model 5 and 6
  ////////////////////////////////////////////  Indexes
  //
    int idx_spyest_to_spye_MC1[N_spyest_iyj_MC1]; // Jumping index from iyj to iy in Model Compartment 1
  int idx_yest_yj_to_ye_y_MC2[N_yest_yj_M]; // Jumping index from yj to y in Model Compartment 2
  //
    int idx_spye_iy_M1[N_obs_Z_M1]; // Species * Year index in Model 1
  int idx_sp_i_M3[N_obs_W_M3]; // Species index in Model 3
  int idx_sp_i_M4[N_obs_W_M4]; // Species index in Model 4
  int idx_spyest_iyj_M4[N_obs_W_M4]; // Species * Year * Station index in in Model 4
  int idx_stye_yj_M6[N_obs_Y_M6]; // Year * Station index in in Model 6
  int idx_refsp_i_MC1;
  ////////////////////////////////////////////  Parameters
  //
    real X2_sd;
  real tau_eta_mock_mu;
  real tau_eta_mock_sd;
  real tau_eta_samp_mu;
  real tau_eta_samp_sd;
  // ////////////////////////////////////////////  Data for Generate Quantities
  // int gq_samp_R[N_obs_Y_M6]; //Total number of reads that multinomial distribution draws from in Model 6
  // int gq_mock_R[N_obs_Y_M5]; //Total number of reads that multinomial distribution draws from in Model 6
}
parameters {
  //////////////////////////////////////////// Model Compartment 1
  vector[N_spyest_iyj_MC1] X1;
  vector[N_sp_i_MC1] theta_1;
  vector<lower=0>[N_spye_iy_MC1] phi_1;
  vector[N_sp_i_MC1] beta_0;
  vector[N_sp_i_MC1] beta_1;
  //////////////////////////////////////////// Model Compartment 2
  vector<lower=0>[N_sp_i_MC2-1] tau_eta_M5;
  vector<lower=0>[N_sp_i_MC2-1] tau_eta_M6;
  vector[N_sp_i_MC2-1] alpha_raw;
  matrix[N_sp_i_MC2-1,N_obs_Y_M6] X2_raw;
  matrix[N_sp_i_MC2-1,N_obs_Y_M5] eta_raw_M5;
  matrix[N_sp_i_MC2-1,N_obs_Y_M6] eta_raw_M6;
  matrix<lower=0>[N_sp_i_MC2-1,N_ye_y_M] phi_2;
  vector[N_sp_i_MC2-1] theta_2;
  //
}
transformed parameters{
  //////////////////////////////////////////// Declaration of transformed parameters
  vector[N_yest_yj_M] C_i_ref = X1[1:N_yest_yj_M]+theta_1[idx_refsp_i_MC1]; //This can be automated better
  matrix[N_sp_i_MC2,N_obs_Y_M6] gamma_M6;
  gamma_M6[N_sp_i_MC2,] = to_row_vector(rep_vector(0.0,N_obs_Y_M6));
  matrix[N_sp_i_MC2,N_obs_Y_M5] gamma_M5;
  matrix[N_sp_i_MC2,N_obs_Y_M6] psi_M6;
  matrix[N_sp_i_MC2,N_obs_Y_M5] psi_M5;
  matrix[N_sp_i_MC2,N_obs_Y_M6] eta_M6 = rep_matrix(0.0,N_sp_i_MC2,N_obs_Y_M6);
  matrix[N_sp_i_MC2,N_obs_Y_M5] eta_M5 = rep_matrix(0.0,N_sp_i_MC2,N_obs_Y_M5);
  vector[N_sp_i_MC2] alpha;
  matrix[N_sp_i_MC2-1,N_yest_yj_M] X2;
  matrix[N_sp_i_MC2-1,N_ye_y_M] V2;
  vector[N_spye_iy_MC1] V1;
  vector[N_obs_W_M3] omega_M3;
  vector[N_obs_W_M4] omega_M4;
  //////////////////////////////////////////// Adjusting priors
  // Adjusting Alpha prior
  alpha[1:(N_sp_i_MC2-1)] = alpha_raw * 0.01;
  alpha[N_sp_i_MC2] = 0;
  // Adjusting X2 prior !!!!
    for (i in 1:N_obs_Y_M6){
      X2[,idx_stye_yj_M6[i]] = X2_raw[,i]*X2_sd;
    }
  // Adjusting Eta priors for samples
  for (i in 1:(N_sp_i_MC2-1)) {
    eta_M6[i,] = eta_raw_M6[i,] * tau_eta_M6[i];
  }
  // Adjusting Eta priors for mock
  for (i in 1:(N_sp_i_MC2-1)) {
    eta_M5[i,] = eta_raw_M5[i,] * tau_eta_M5[i];
  }
  // Averaging X1 across stations -> V1
  V1 = to_vector(average_by_idx(N_spyest_iyj_MC1, N_spye_iy_MC1, to_array_1d(X1), idx_spyest_to_spye_MC1));
  // Averaging X2 across stations -> V2
  for(i in 1:N_sp_i_MC2-1){
    for(j in 1:N_ye_y_M){
      V2[i,] = to_row_vector(average_by_idx(N_yest_yj_M, N_ye_y_M, to_array_1d(X2[i,]), idx_yest_yj_to_ye_y_MC2));
    }
  }
  //
    //////////////////////////////////////////// Model Compartment 1
  // Model 3
  for (i in 1:N_obs_W_M3){
    omega_M3[i] = beta_0[idx_sp_i_M3[i]] + (beta_1[idx_sp_i_M3[i]] * S_M3[i]);
  }
  // Model 4
  for (i in 1:N_obs_W_M4){
    omega_M4[i] = beta_0[idx_sp_i_M4[i]] + (beta_1[idx_sp_i_M4[i]]*
                                              (X1[idx_spyest_iyj_M4[i]]+theta_1[idx_sp_i_M4[i]])); //
  }
  //////////////////////////////////////////// Model Compartment 2
  // Model 5 (Gamma)
  for (i in 1:N_obs_Y_M5){
    for(j in 1:N_sp_i_MC2){
      gamma_M5[j,i] = alr_M5[j]+(NPCR*(alpha[j]))+eta_M5[j,i];
    }
  }
  // Model 6 (Gamma)
  for (i in 1:N_obs_Y_M6){
    for(j in 1:N_sp_i_MC2-1){
      gamma_M6[j,i] = ((X2[j,idx_stye_yj_M6[i]]+theta_2[j])-C_i_ref[idx_stye_yj_M6[i]])+(NPCR*(alpha[j]))+eta_M6[j,i];
    }
  }
  // Model 5 (Psi)
  for (i in 1:N_obs_Y_M5){
    psi_M5[,i] = softmax(gamma_M5[,i]);
  }
  // Model 6 (Psi)
  for (i in 1:N_obs_Y_M6){
    psi_M6[,i] = softmax(gamma_M6[,i]);
  }
}
model {
  //////////////////////////////////////////// Model compartment 1
  // Model 1
  for (i in 1:N_obs_Z_M1) {
    Z_M1[i] ~ neg_binomial_2(exp(V1[idx_spye_iy_M1[i]])*E_M1[i],
                             phi_1[idx_spye_iy_M1[i]]);
  }
  // Model 2
  for (i in 1:N_sp_i_MC2-1) {
    for (j in 1:N_obs_Z_M2) {
      Z_M2[i,j] ~ neg_binomial_2(exp(V2[i,idx_yest_yj_to_ye_y_MC2[j]])*E_M2[j],
                                 phi_2[i,idx_yest_yj_to_ye_y_MC2[j]]);
    }
  }
  // Model 3
  W_M3 ~ binomial(U_M3, inv_cloglog(omega_M3));
  // Model 4
  W_M4 ~ binomial(U_M4, inv_cloglog(omega_M4));
  //////////////////////////////////////////// Model compartment 2
  // Model 5
  for (i in 1:N_obs_Y_M5){
    Y_M5[,i] ~ multinomial(psi_M5[,i]);
  }
  // Model 6
  for (i in 1:N_obs_Y_M6){
    Y_M6[,i] ~ multinomial(psi_M6[,i]);
  }
  //////////////////////////////////////////// Priors
  beta_0 ~ normal(0, 10);
  beta_1 ~ normal(0, 10);
  theta_1 ~ normal(0,10); //
    phi_1 ~ gamma(50,1); //
    theta_2 ~ normal(0,10); //
    for(i in 1:N_sp_i_MC2-1){
      phi_2[i,] ~ gamma(50,1); //
    }
  // Alpha Prior
  alpha_raw ~ std_normal();
  // Eta Prior
  for(i in 1:(N_sp_i_MC2-1)){
    eta_raw_M5[i,] ~ std_normal();
    eta_raw_M6[i,] ~ std_normal();
  }
  tau_eta_M6 ~ gamma(tau_eta_samp_mu,tau_eta_samp_sd);
  tau_eta_M5 ~ gamma(tau_eta_mock_mu,tau_eta_mock_sd);
  for(i in 1:(N_sp_i_MC2-1)){
    X2_raw[i,] ~ std_normal();
  }
}
generated quantities{
  matrix[N_sp_i_MC2-1,N_yest_yj_M] C2;
  vector[N_spyest_iyj_MC1] C1;
  // New_C1
  for (i in 1:N_obs_W_M4){
    C1[idx_spyest_iyj_M4[i]]=X1[idx_spyest_iyj_M4[i]]+theta_1[idx_sp_i_M4[i]]; //
  }
  // New_C2
  for (i in 1:N_obs_Y_M6){
    for(j in 1:N_sp_i_MC2-1){
      C2[j,idx_stye_yj_M6[i]] = X2[j,idx_stye_yj_M6[i]]+theta_2[j];
    }
  }
  // // Declaration of variables
  // int Z_M1_new[N_obs_Z_M1_new];
  // int Z_M2_new[N_sp_i_MC2_new-1,N_obs_Z_M2_new];
  //
    // // Generated quantities for Trawl Model 1 and Model 2
  // for (i in 1:N_obs_Z_M1_new) {
    //   Z_M1_new[i] = neg_binomial_2_rng(exp(V1_new[idx_spye_iy_M1_new[i]])*E_M1_new[i],
                                          //   phi_1_new[idx_spye_iy_M1_new[i]]);
    // }
  // for (i in 1:N_sp_i_MC2_new-1) {
    //   for (j in 1:N_obs_Z_M2_new) {
      //     Z_M2_new[i,j] = neg_binomial_2_rng(exp(V2_new[i,idx_yest_yj_to_ye_y_MC2_new[j]])*E_M2_new[j],
                                                //     phi_2_new[i,idx_yest_yj_to_ye_y_MC2_new[j]]);
      //   }
  } 