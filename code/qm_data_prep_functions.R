### THIS SCRIPT INCLUDES HELPER FUNCTIONS FOR THE QM/qPCR JOINT MODEL TO HELP TRANSLATE/FEED INTO STAN
### INCLUDING FUNCTIONS TO HELP WITH FORMATTING, DATA TRANSFORMS AND INITIAL VALUES

suppressPackageStartupMessages(library(mgcv))

# code for parsing smoothers
source(here("code","smoothers.R"),local=TRUE)

# qPCR to Stan ------------------------------------------------------------------------------
# Slightly more complex than MB because we have to do a bit more formatting; optionally create smoothers here
# Inputs: qPCR unknown samples, and the qPCR standards

format_qPCR_data <- function(qPCR_unknowns, qPCR_standards,tube_dat,
                             unk_covariates=NULL,cov_type=NULL,unk_smoothes=NULL,
                             #unk_random = NULL,
                             unk_offsets=NULL,offset_type = rep("log",length(unk_offsets))){
  
  # Make matrices for qPCR smooths and covariates
  # This is just for marginal effects, not interactions.
  FORM <- list()
  X_cov <- list()
  SM_FORM <- list()
  SM <-list()
  X_offset <- NULL
  
  if(is.null(unk_covariates) ==FALSE){
    for(i in 1:length(unk_covariates)){
      if(cov_type[i] =="continuous"){
        FORM[[unk_covariates[i]]] <- paste("hake_Ct ~ 0+", unk_covariates[i])
        model_frame   <- model.frame(FORM[[unk_covariates[i]]], qPCR_unknowns)  
        X_cov[[unk_covariates[i]]] <- model.matrix(as.formula(FORM[[i]]), model_frame)
      }else if(cov_type[i] =="factor"){
        FORM[[unk_covariates[i]]] <- paste("hake_Ct ~ 0 + factor(", unk_covariates[i],")")
        model_frame   <- model.frame(FORM[[unk_covariates[i]]], qPCR_unknowns)  
        X_cov[[unk_covariates[i]]] <- model.matrix(as.formula(FORM[[unk_covariates[i]]]), model_frame)
      }
    }
  }
  if(is.null(unk_smoothes)==FALSE){
    for(i in 1:length(unk_smoothes)){
      SM_FORM[[unk_smoothes[i]]] <- paste0("hake_Ct ~ s(", unk_smoothes[i],")")
      # Model form for smoothes
      #FORM.smoothes <- "copies_ul ~ s(bottom.depth.consensus,by=year,k=4)"
      SM[[unk_smoothes[i]]] <- parse_smoothers(eval(SM_FORM[[unk_smoothes[i]]]) ,qPCR_unknowns)
      # Objects you care about in SM are:
      # Zs basis function matrices
      # Xs; // smoother linear effect matrix
      # SM$basis_out is the basis function
      
      # This is for making predictions to new data using the old basis 
      # new_smooth_pred <- parse_smoothers(eval(SM_FORM[[i]]),data=qPCR_unknowns,
      #                              #newdata= NEWDATA,
      #                              basis_prev = SM$basis_out)
      
      # other things that are somewhat convenient (mostly ported from TMB, not all relevant.)
      # n_bs     <- ncol(SM$Xs)
      # b_smooth_start <- SM$b_smooth_start
      # n_smooth <- length(b_smooth_start)
      # b_smooth <- if (SM$has_smooths) rep(0,sum(SM$sm_dims)) else array(0) 
      # has_smooths <- SM$has_smooths
    }
  }
  
  if(is.null(unk_offsets) ==FALSE){
    for(i in 1:length(unk_offsets)){
      if(offset_type[i] =="log"){
        if(i == 1){
          X_offset <-qPCR_unknowns[,unk_offsets[i]] %>% log() 
        } else{
          X_offset <- cbind(X_offset, qPCR_unknowns[,unk_offsets[i]] %>% log())
        }
        colnames(X_offset)[i] = paste0("log_",unk_offsets[i])
      }else{print("Offset type not supported at present")}
      
    }
    # Collapse multiple offsets into single vector
    X_offset_tot = rowSums(X_offset)
  }
  
  # pull covars (not sure how to generalize this yet)
  # X_cov <- map_df(X_cov,~bind_cols(as.numeric))
  # X_offset <- bind_cols(X_offset)
  
  # format unknowns
  qPCR_unk <- qPCR_unknowns %>% 
    # pick columns we care about and rename
    # select(qPCR, well,tubeID,type,task,IPC_Ct,inhibit.val,inhibit.bin,wash_idx_obs=wash_idx,
    #        dilution,depth_cat,Ct=hake_Ct,copies_ul=hake_copies_ul) %>% 
    rename(wash_idx_obs=wash_idx,
           Ct=hake_Ct,copies_ul=hake_copies_ul) %>% 
    mutate(z=ifelse(Ct=="Undetermined",0,1)) %>%
    mutate(Ct=str_replace_all(Ct,"Undetermined",'')) %>% 
    mutate(Ct=as.numeric(Ct) %>% round(2)) %>% 
    mutate(Ct = ifelse(is.na(Ct), 99, Ct)) %>%  # Stan doesn't like NAs
    filter(task=="UNKNOWN",type=="unknowns") %>%  #this shouldn't really filter anything out (i.e., the data should have been cleaned before this step)
    # INDEX OF UNIQUE PLATES
    mutate(plate_idx=match(qPCR,unique(qPCR))) %>% 
    # INDEX OF UNIQUE SAMPLES
    mutate(tube_idx=match(tubeID,unique(tubeID))) 
  # WASH AND DILUTION EFFECTS
  #mutate(log_dilution=log(dilution))
  
  # as of 08.21.24, 31 unique plates, 1818 unique samples
  
  #standards
  qPCR_std <- qPCR_standards %>% 
    # pick columns we care about and rename
    select(qPCR, well,tubeID,type,task=hake_task,IPC_Ct,Ct=hake_Ct,copies_ul=hake_copies_ul) %>% 
    mutate(z=ifelse(Ct=="Undetermined",0,1)) %>%
    mutate(Ct=str_replace_all(Ct,"Undetermined",'')) %>% 
    mutate(Ct=as.numeric(Ct) %>% round(2))%>% 
    mutate(Ct = ifelse(is.na(Ct), 99, Ct)) %>%  # Stan doesn't like NAs
    # there are 7 samples where the task is "UNKNOWN" instead of "STANDARD", but that is a lab error. let's fix these
    mutate(copies_ul=case_when(
      task=="STANDARD" ~ copies_ul,
      task=="UNKNOWN"&tubeID=="E00"~1,
      task=="UNKNOWN"&tubeID=="E01"~10,
      task=="UNKNOWN"&tubeID=="5"~5
    )) %>% 
    # then we can call them all standards
    mutate(task="STANDARD") %>% 
    #had an issue with non-unique sample names between stds and unks because of '5' being used as sample id for standards with conc. of 5 copies
    mutate(tubeID=ifelse(tubeID=="5","5C",tubeID)) %>% 
    # add the plate index from the unknowns
    left_join(distinct(qPCR_unk,qPCR,plate_idx),by=join_by(qPCR)) %>% 
    # finally, remove plates that are in standards but not qpcr (H8,H17,H24, which we checked are plates with errors and full sets of controls with no unknowns)
    filter(!is.na(plate_idx))
  
  qPCRdata <- list(qPCR_unk = qPCR_unk,
                   qPCR_std = qPCR_std,
                   tube_dat = tube_dat,
                   FORM = FORM,
                   X_cov = X_cov,
                   X_offset = X_offset,
                   X_offset_tot = X_offset_tot,
                   SM_FORM = SM_FORM,
                   SM = SM)
  return(qPCRdata)
}


# make Stan data inputs for qPCR, including making design matrices

prepare_stan_data_qPCR <- function(qPCRdata){
  
  unk <- qPCRdata$qPCR_unk
  std <- qPCRdata$qPCR_std
  tube_dat <- qPCRdata$tube_dat
  
  ## Design matrix for qPCR data
  
  form <- "depth_cat ~ 0 + factor(station_depth_idx)"
  model_frame   <- model.frame(form, tube_dat)  
  X_station_depth_tube <- model.matrix(as.formula(form), model_frame)
  
  ### RANDOM EFFECT DESIGN MATRICES
  # Make a matrix for random effect associated with each station-depth combination at observation level 
  form <- "year ~ 0 + factor(tubeID): factor(n_tube_station_depth)"
  model_frame   <- model.frame(form, unk)
  A <- model.matrix(as.formula(form), model_frame)
  # get rid of factor levels (columns) that == 0
  X_bio_rep_obs <- A[,which(colSums(A)>0)]
  
  ### Make random effect that sums to zero for the tubes.
  bio_rep_dat <- tube_dat %>% group_by(station_depth_idx) %>% mutate(bio_rep_idx = rep(1:n())) %>%  ungroup()
  bio_rep_dat2 <- bio_rep_dat %>% filter(bio_rep_idx == n_tube_station_depth)
  form <- "depth_cat ~ 0 + factor(tubeID): factor(n_tube_station_depth)"
  model_frame   <- model.frame(form, bio_rep_dat)  
  X_bio_rep_tube <- model.matrix(as.formula(form), model_frame)
  X_bio_rep_tube <- X_bio_rep_tube[,which(colSums(X_bio_rep_tube)>0)]
  
  # figure out how many parameters you actually need to estimate for random effects: RE > param
  N_bio_rep_RE    <- nrow(bio_rep_dat)
  N_bio_rep_param <- bio_rep_dat %>% filter(bio_rep_idx != n_tube_station_depth) %>% nrow()
  bio_rep_idx <- bio_rep_dat2 %>%  pull(bio_rep_idx)
  N_bio_rep_idx <- length(bio_rep_idx)
  
  stan_qPCR_data <- list(
    qpcr_unk = unk, # unknown field samples
    qpcr_std = std, # standards
    Nplates = length(unique(unk$qPCR)), # number of plates
    Nobs_qpcr = nrow(unk), # number of field observations
    NSamples_qpcr = max(unk$tube_idx), # number of unique tubes
    NstdSamples = nrow(std), # number of standard samples
    plate_idx = unk$plate_idx, # index linking field samples to qPCR plates
    std_plate_idx = std$plate_idx, # index linking standards to qPCR plates
    tube_idx = tube_dat$tube_idx, # index identifying unique tubes
    y_unk = unk$Ct, # cycles, unknowns
    z_unk = unk$z, # did it amplify? unknowns
    y_std = std$Ct, # cycles, standards
    z_std = std$z, # did it amplify? standards
    known_concentration = std$copies_ul,# known copy number from standards
    beta_std_curve_0_offset =36,
    stdCurvePrior_intercept = c(39, 3), #normal distr, mean and sd ; hyperpriors
    stdCurvePrior_slope = c(-1.442695, 0.05), #normal distr, mean and sd ; hyperpriors 
    # hard coded covariates and offsets- COULD GENERALIZE THIS LATER
    wash_prior = c(-1,1),
    X_offset_tot =  qPCRdata$X_offset_tot, # dilution and other offsets combined
    # RE Matrices
    X_bio_rep_tube =X_bio_rep_tube,
    X_bio_rep_obs =X_bio_rep_obs,
    X_station_depth_tube=X_station_depth_tube,
    N_bio_rep_RE = N_bio_rep_RE,
    N_bio_rep_param = N_bio_rep_param,
    bio_rep_idx =bio_rep_idx,
    N_bio_rep_idx = N_bio_rep_idx
    # COULD ADD SMOOTHS HERE EVENTUALLY
  )
  
  
  if("wash_idx"%in% names(qPCRdata$X_cov)){
    wash_idx = c(qPCRdata$X_cov[["wash_idx"]])
    
    #   wash_idx <- unk %>%
    #     group_by(qpcr_sample_idx) %>%
    #     summarise(wash_idx=mean(wash_idx_obs)) %>%
    #     pull("wash_idx")
    # 
    stan_qPCR_data <- c(stan_qPCR_data,
                        list(wash_idx = wash_idx)) # design matrix for covariates (right now, just the wash effect)
  }
  if("station_depth_idx"%in%names(unk)){
    qPCRdata$X_cov[["station_depth_idx"]]
    
    X_station_depth_obs <- formatted_qPCR_data$X_cov[["station_depth_idx"]]
    N_station_depth <- ncol(X_station_depth_obs)
    stan_qPCR_data <- c(stan_qPCR_data,
                        list(X_station_depth_obs = X_station_depth_obs,
                             N_station_depth = N_station_depth))
  }
  
  return(stan_qPCR_data)
}

# Metabarcoding to Stan ------------------------------------------------------------------------------
# inputs: cleaned metabarcoding field samples, clean mock communities
format_metabarcoding_data <- function(input_metabarcoding_data, input_mock_comm_data){ #nested within level2 units, e.g., technical replicates of biological replicates (NA if absent)

  Observation <- input_metabarcoding_data
  
  Mock <- input_mock_comm_data  
  
  # index species to a common standard 
  sp_list <- data.frame(
    species = c(Mock$species, Observation$species) %>% unique(),
    species_idx = NA)
  sp_list$species_idx <- match(sp_list$species, unique(sp_list$species)) 
  
  Observation <- Observation %>% 
    left_join(sp_list,by=join_by(species))
  
  Mock_unique <- Mock %>% distinct(Sample) %>% 
                  arrange(Sample) %>% mutate(mockID = 1:nrow(.))
  Mock <- left_join(Mock,Mock_unique,by=join_by(Sample)) %>% ungroup() %>% 
           left_join(sp_list,by=join_by(species))
    # # make a combined mockID
    # arrange(Mock_name,Mock_type,Rep) %>% 
    # group_by(Mock_name,Mock_type,Rep) %>% 
    # mutate(mockID=cur_group_id()) %>% 
    # ungroup()
    # 
  return(
    metabarcoding_data <- list(
      Observation = Observation,
      Mock = Mock,
      N_pcr_mock = 43, 
      NSpecies = nrow(sp_list),
      sp_list = sp_list
    ))
}

# additive log-ratio transform of a matrix; defined here, used in makeDesign()
alrTransform <- function(formatted_mock){
  
  # number of species
  nspp <- length(unique(formatted_mock$species))
  
  # wide form to long form, and fill in 1e-09 for "zeroes" to enable the log ratio calculation
  p_mock <- formatted_mock %>%
    select(mockID,species_idx,b_proportion) %>% 
    distinct(mockID,species_idx,b_proportion) %>%  
    pivot_wider(names_from = species_idx, names_sort=T,names_prefix="alr_",values_from = b_proportion, values_fill = 1e-9) %>% 
    ungroup() %>% 
    arrange(mockID)

  p_mock <- left_join(formatted_mock %>% distinct(mockID,Rep) %>%  select(mockID),
                      p_mock,by=join_by(mockID))
  
  colnames(p_mock)[2:ncol(p_mock)] <- paste0("alr_", 1:nspp)
  
  # calculate the ALRs
  p_mock <- alr(p_mock[,2:ncol(p_mock)]) %>% as.matrix() %>% as.data.frame()
  
  p_mock[,nspp] <- 0  #add zero expressly for reference species
  names(p_mock)[nspp] <- paste0("alr_", nspp)
  
  p_mock <-  cbind(formatted_mock %>% select(mockID) %>% distinct(),
                   p_mock) %>% ungroup()
  
  return(p_mock)
}
  
## Make design matrices for Stan from the metabarcoding field samples

prepare_stan_data_MB <- function(obs, #obs is the list output from format_metabarcoding_data
                       qPCR_tube_obs, # qpcr field samples formatted
                       N_pcr_cycles){ #N_pcr_cycles is the number of PCR cycles in your experimental/enviro samples; currently a single value, could be made into a vector if this number varies
  require(MCMCpack)
  require(compositions)
  require(rstan)
  require(dplyr)
    
  MOCK <- obs$Mock
  OBSERVED <- obs$Observation
  
  rep_level_mock <- "Rep"  #name of the column at the lowest level of replication
  rep_level_samp <- "tech" #name of the column at the lowest level of replication
  
  # ALR transform using the function above
  p_mock_all <- alrTransform(MOCK)
  
  MOCK <- MOCK %>% 
    dplyr::select(species_idx,Mock_name=Sample,Mock_type,mockID,  #unique biological samples
                  all_of(rep_level_mock),  #lowest level of replication
                  Nreads) %>% 
    ungroup() %>%
    pivot_wider(names_from = species_idx, names_prefix="sp_",names_sort=T,values_from = Nreads, values_fill = 0)
  
  N_pcr_mock <- rep(obs$N_pcr_mock, nrow(p_mock_all)) #assumes all have the same Npcr
  
  p_samp_all <- OBSERVED %>% 
    ungroup() %>% 
    dplyr::select(species_idx,
                  tubeID,  #unique biological samples
                  all_of(rep_level_samp),  #lowest level of replication
                  nReads) %>%
    ungroup() %>% 
    pivot_wider(names_from = species_idx, names_prefix="sp_",names_sort=T,values_from = nReads, values_fill = 0)
  
  # HERE IS WHERE YOU CHECK TO MAKE SURE THE SAMPLES IN THE METABARCODING ARE PRESENT IN THE QPCR DATA.
  # OTHERWISE YOU HAVE TO EXCLUDE THEM IN THE JOINT MODEL.
  
  p_samp_all <- p_samp_all %>%  filter(tubeID %in% tube_dat$tubeID)
  # removed 3 samples
  
  N_pcr_samp <- rep(N_pcr_cycles, nrow(p_samp_all))
  
  # make link to qPCR data here.
  p_samp_all <- qPCR_tube_obs %>% 
    distinct(tubeID,tube_idx,station_depth_idx,station_idx) %>%
    left_join(p_samp_all,by=join_by(tubeID)) %>% 
    filter(!is.na(sp_1)) 
  
  tube_link_idx = p_samp_all$tube_idx
  
  # Make design matrices for mocks and samps
  
  NOM <- as.name(colnames(p_mock_all)[1])
  formula_a <- eval(NOM) ~ N_pcr_mock -1
  model_frame <- model.frame(formula_a, p_mock_all)
  model_vector_a_mock <- model.matrix(formula_a, model_frame) %>% as.numeric()

  N_obs_mock <- nrow(p_mock_all)
  
  # unknown communities
  # species compositions (betas)
  
  NOM <- as.name(colnames(p_samp_all)[1])    
  p_samp_all$tubeID <- as.factor(p_samp_all$tubeID) 
  N_S = length(unique(p_samp_all$tubeID))
  p_samp_all[rep_level_samp] <- as.factor(unlist(p_samp_all[rep_level_samp]))
  if(N_S == 1){
    formula_b <- eval(NOM) ~ 1  
  } else {
    formula_b <- eval(NOM) ~ 0+tubeID
  }
  
  model_frame <- model.frame(formula_b, p_samp_all)
  model_matrix_b_samp <- model.matrix(formula_b, model_frame)
  
  # efficiencies (alphas)
  formula_a <- eval(NOM) ~ N_pcr_samp -1
  model_frame <- model.frame(formula_a, p_samp_all)
  model_vector_a_samp <- model.matrix(formula_a, model_frame) %>% as.numeric()
  
  #counters 
  # N_obs_samp_small <- nrow(model_matrix_b_samp_small)
  N_obs_samp <- nrow(p_samp_all)
  N_b_samp_col <- ncol(model_matrix_b_samp)
  
  # Make index for most abundant species in each metabarcoding sample.
  which.max <- function(a){return(which(a == max(a)))}
  
  ref_sp_idx <- p_samp_all %>% dplyr::select(contains("sp")) %>%
    apply(.,1,which.max) %>% as.data.frame()
  colnames(ref_sp_idx) <- "ref_sp_idx"
  ref_sp_idx <- bind_cols(p_samp_all,ref_sp_idx)  
  
  #### Make Stan objects
  stan_data <- list(
    sp_list = obs$sp_list,
    N_species = ncol(p_samp_all %>% dplyr::select(contains("sp"))),   # Number of species in data
    N_obs_mb_samp = nrow(p_samp_all), # Number of observed community samples and tech replicates ; this will be Ncreek * Nt * Nbiol * Ntech * 2 [for upstream/downstream observations]
    N_obs_mock = nrow(p_mock_all), # Number of observed mock samples, including tech replicates
    
    # Observed data of community matrices
    sample_data_labeled = p_samp_all ,
    sample_data = p_samp_all %>% dplyr::select(contains("sp")),
    ref_sp_idx = ref_sp_idx$ref_sp_idx,
    tube_link_idx = tube_link_idx,

    mock_data_labeled   = MOCK ,
    mock_data   = MOCK %>% dplyr::select(contains("sp")),

    # True proportions for mock community
    #mock_true_prop = p_mock_all %>% dplyr::select(contains("sp")),
    alr_mock_true_prop = p_mock_all %>% dplyr::select(contains("alr")),
    
    # Dirichlet alpha0
    log_dm_alpha0_mock = log(1000),
    model_vector_a_samp = model_vector_a_samp,
    
    # Design matrices: mock community samples
    model_vector_a_mock = as.array(model_vector_a_mock),
    
    # Priors
    alpha_prior = c(0,0.01),  # normal prior
    # tau_prior = c(10,1000),   # gamma prior on eta_mock = ~0.01
    log_D_mu = c(0),
    log_D_scale = c(5)
  )
  return(stan_data)
}
  
# Setting initial values ------------------------------------------------------------------------------

make_stan_inits <- function(n.chain,
                         jointData,
                         log_D_link_sp_init_mean){
  
  sample_data <- jointData$sample_data
  ref_col = jointData$qpcr_mb_link_sp_idx
  
  p <-  ((sample_data) / rowSums(sample_data)) + 1e-5
  log_p <- log(p / rowSums(p))
  
  log_p_rel <- log_p - log_p[,ref_col]
  log_D_init <- log_p_rel + log_D_link_sp_init_mean
  
  N_obs_mb = jointData$N_obs_mb_samp
  N_species = jointData$N_species
  Nplates=jointData$Nplates
  N_station_depth=jointData$N_station_depth
  
  # populate the list of inits
  A <- list()
  for(i in 1:n.chain){
    A[[i]] <- list(
      log_D_link_sp_init_mean = log_D_link_sp_init_mean, # initial mean log conc of qpcr link species
      log_D_init = log_D_init,
      mean_hake  = rnorm(1,log_D_link_sp_init_mean,0.01),
      log_D_station_depth=rnorm(N_station_depth,0,0.01),
      alpha_raw=jitter(rep(0,N_species-1),factor=1),
      beta_std_curve_0=runif(Nplates,-2,2),
      beta_std_curve_1=runif(Nplates,-1.44,-1.3),
      phi_0 = runif(1,1.5,1.8),
      phi_1 = runif(1,1,1.1),
      gamma_1 = runif(1,-0.01,0)
      # tau=rgamma(1,10,1000),
      # log_dm_alpha0_mock=log(1000)
      #eta_mock_raw = matrix((rnorm(N_species-1)*N_obs_mock),N_obs_mock,N_species-1)
    )
  }  
  return(A)
}
