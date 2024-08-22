### THIS SCRIPT INCLUDES HELPER FUNCTIONS FOR THE QM/qPCR JOINT MODEL TO HELP TRANSLATE/FEED INTO STAN
### INCLUDING FUNCTIONS TO HELP WITH FORMATTING, DATA TRANSFORMS AND INITIAL VALUES

library(mgcv)

# code for parsing smoothers
source(here("code","smoothers.R"),local=TRUE)

### FORMAT THE METABARCODING DATA FOR INGESTION INTO STAN
format_metabarcoding_data <- function(input_metabarcoding_data, input_mock_comm_data,
                                      Level_1_treatment_envir, #e.g., unique sampling site
                                      Level_2_treatment_envir, #nested within level1 units, e.g., unique biological samples
                                      Level_3_treatment_envir,  #nested within level2 units, e.g., technical replicates of biological samples (NA if absent)
                                      Level_1_treatment_mock, #e.g., unique biol community
                                      Level_2_treatment_mock, #nested within level1 units, e.g., unique biological replicates
                                      Level_3_treatment_mock  #nested within level2 units, e.g., technical replicates of biological replicates (NA if absent)
                                      ){
  require(tidyverse)
  
  Observation <- input_metabarcoding_data
    
  Mock <- input_mock_comm_data %>% 
    replace(is.na(.), 0) %>% 
    filter(b_proportion > 0) %>% 
    filter(Nreads > 0) %>% #omit things that are absent from the mocks
    filter(!is.na(Nreads)) %>% 
    rename(level1 = as.name(Level_1_treatment_mock),
           level2 = as.name(Level_2_treatment_mock))
  
  if (is.na(Level_3_treatment_mock)) {Mock$level3 <- 1} else {Mock <- Mock %>% rename(level3 = as.name(Level_3_treatment_mock))}  
  
  Mock <- Mock %>% 
    unite(c(level1, level2, level3), col = "Sample", sep = ".", remove = F)
    
  # only keep species present in both mocks and observations
  keepSpecies <- intersect(Mock$species, Observation$species)
  Observation <- Observation %>% 
    filter(species %in% keepSpecies) 
  Mock <- Mock %>% filter(species %in% keepSpecies)
  
  # index species to a common standard 
  sp_list <- data.frame(
    species = c(Mock$species, Observation$species) %>% unique(),
    species_idx = NA)
  sp_list$species_idx <- match(sp_list$species, unique(sp_list$species)) 
  
  #reindex and renormalize to deal with omitted species
  Mock <- Mock %>% 
    left_join(sp_list) %>% 
    mutate(level1 = match(level1, unique(level1)),
           speciesname = species,
           species = species_idx) %>% 
    group_by(Sample) %>% 
    mutate(b_proportion = b_proportion/sum(b_proportion)) %>% 
    ungroup()
  
  Observation <- Observation %>% 
    rename(level1 = as.name(Level_1_treatment_envir),
           level2 = as.name(Level_2_treatment_envir))
  
  if (is.na(Level_3_treatment_envir)) {Observation$level3 <- 1} else {Observation <- Observation %>% rename(level3 = as.name(Level_3_treatment_envir))}  
    
    
  Observation <- Observation %>% 
    left_join(sp_list) %>% 
    mutate(speciesname = species,
           species = species_idx) %>% 
    # mutate(station = case_when(stationname == "Dn" ~ 1,
    #                            stationname == "Up" ~ 2)) %>% 
    group_by(level1) %>% 
    mutate(level2 = match(level2, unique(level2))) %>% #reindex, if necess
    mutate(level3 = match(level3, unique(level3))) %>% #reindex, if necess
    unite(c(level1, level2, level3), col = "Sample", sep = ".", remove = F)
  
  station_list <- data.frame(
    station = Observation$level1 %>% unique(),
    station_idx = 1:length(unique(Observation$level1)))
  
  
  #list object containing named elements Observation, Mock, and Npcr
  
  return(
  metabarcoding_data <- list(
    Observation = Observation,
    Mock = Mock,
    N_pcr_mock = 43, 
    NSpecies = nrow(sp_list),
    station_list = station_list,
    sp_list = sp_list
  ))
}

# prep qPCR data
# function takes two dataframe inputs, for qPCR unknown samples, and the qPCR standards
format_qPCR_data <- function(qPCR_unknowns, 
                             qPCR_standards,
                             unk_covariates=NULL,cov_type=NULL,
                             unk_smoothes=NULL,
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
      } else if(cov_type[i] =="factor"){
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
      if(offset_type =="log"){
        if(i == 1){
          X_offset <-qPCR_unknowns[,unk_offsets[i]] %>% log() 
        } else{
          X_offset <- cbind(X_offset, qPCR_unknowns[,unk_offsets[i]] %>% log())
        }
        colnames(X_offset)[i] = paste0("log_",unk_offsets[i])
      }else{print("Offset type not supported at present")}
      
    }
  }
  
  # pull covars (not sure how to generalize this yet)
  # X_cov <- map_df(X_cov,~bind_cols(as.numeric))
  # X_offset <- bind_cols(X_offset)
  
  #unknowns
  qPCR_unk <- qPCR_unknowns %>% 
    # pick columns we care about and rename
    select(qPCR, well,tubeID,type,task,IPC_Ct,inhibit.val,inhibit.bin,wash_idx_obs=wash_idx,
           dilution,depth_cat,Ct=hake_Ct,copies_ul=hake_copies_ul) %>% 
    mutate(z=ifelse(Ct=="Undetermined",0,1)) %>%
    mutate(Ct=str_replace_all(Ct,"Undetermined",'')) %>% 
    mutate(Ct=as.numeric(Ct) %>% round(2)) %>% 
    mutate(Ct = ifelse(is.na(Ct), 99, Ct)) %>%  # Stan doesn't like NAs
    filter(task=="UNKNOWN",type=="unknowns") %>%  #this shouldn't really filter anything out (i.e., the data should have been cleaned before this step)
    # INDEX OF UNIQUE PLATES
    mutate(plate_idx=match(qPCR,unique(qPCR))) %>% 
    # INDEX OF UNIQUE SAMPLES
    mutate(qpcr_sample_idx=match(tubeID,unique(tubeID))) %>% 
    # WASH AND DILUTION EFFECTS
    mutate(log_dilution=log(dilution))
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
                   FORM = FORM,
                   X_cov = X_cov,
                   X_offset = X_offset,
                   SM_FORM = SM_FORM,
                   SM = SM)
  return(qPCRdata)
}

# a last piece is we need a sample identifier across QM and qPCR data
# We need to link unique qPCR biological samples to unique QM samples
prepare_stan_qPCR_mb_join <- function(input_metabarcoding_data,unk_formatted,link_species="Merluccius productus"){
  
  mb_link_1 <- input_metabarcoding_data %>% 
    mutate(mb_link=match(Sample,unique(Sample))) %>% 
    distinct(Sample,.keep_all = T)
  
  mb_link_2 <- unk_formatted %>% 
    distinct(plate_idx,qpcr_sample_idx,.keep_all = T) %>% 
    left_join(mb_link_1,by=join_by(tubeID==Sample)) %>% 
    dplyr::select(qpcr_sample_idx,mb_link) %>% 
    arrange(qpcr_sample_idx) %>% 
    filter(!is.na(mb_link))
  
  # get the index of the right species
  qpcr_mb_link_sp_idx <- match(link_species,sort(unique(input_metabarcoding_data$species)))
  
  # return just the link vector
  return(list(
    # datout = mb_link_2,
    N_mb_link = nrow(mb_link_2), # length of the linking vector
    mb_link_sp_idx= qpcr_mb_link_sp_idx,
    mb_link_idx=mb_link_2$mb_link))
}

# make stan data for qPCR part
prepare_stan_data_qPCR <- function(qPCRdata){
  
  unk <- qPCRdata$qPCR_unk
  std <- qPCRdata$qPCR_std
  
  if("wash_idx_obs"%in%names(unk)){
    wash_idx <- unk %>% 
      group_by(qpcr_sample_idx) %>% 
      summarise(wash_idx=mean(wash_idx_obs)) %>% 
      pull("wash_idx")
  }
  log_dil <-unk %>% 
    group_by(qpcr_sample_idx) %>% 
    summarise(log_dil=mean(log_dilution)) %>% 
    pull(log_dil)
  
  stan_qPCR_data <- list(
    Nplates = length(unique(unk$qPCR)),
    Nobs_qpcr = nrow(unk),
    NSamples_qpcr = length(unique(unk$qpcr_sample_idx)),
    NstdSamples = nrow(std),
    plate_idx = unk$plate_idx,
    std_plate_idx=std$plate_idx,
    qpcr_sample_idx = unk$qpcr_sample_idx,
    y_unk = unk$Ct, # cycles, unknowns
    z_unk = unk$z, # did it amplify? unknowns
    y_std = std$Ct, # cycles, standards
    z_std = std$z, # did it amplify? standards
    known_concentration = std$copies_ul,# known copy number from standards
    stdCurvePrior_intercept = c(39, 3), #normal distr, mean and sd ; hyperpriors
    stdCurvePrior_slope = c(-3, 1), #normal distr, mean and sd ; hyperpriors
    # hard coded covariates and offsets- COULD GENERALIZE THIS LATER
    wash_idx = wash_idx, # design matrix for covariates (right now, just the wash effect)
    wash_prior = c(-1,1),
    log_dil = log_dil # dilution offset
    # COULD ADD SMOOTHS HERE EVENTUALLY
  )
  
  return(stan_qPCR_data)
  
}

### END qPCR part ###

  # additive log-ratio transform of a matrix; defined here, used in makeDesign()
alrTransform <- function(MOCK, 
                           Nlevels = 2 #if both biol and tech replication, levels = 3, otherwise levels = 2
                           ){
    require(tidyverse)
    require(compositions)
    
    MOCK <- MOCK %>% 
      unite(c("level1", "level2", "level3")[1:Nlevels-1], col = S, sep = "_", remove = F)
    
      p_mock <- MOCK %>% 
        dplyr::select(species, S, c("level1", "level2", "level3")[Nlevels], b_proportion) %>% 
        pivot_wider(names_from = species, values_from = b_proportion, values_fill = 1e-9) %>% 
        ungroup() 

    colnames(p_mock)[3:(length(unique(MOCK$species))+2)] <- paste0("alr_", 1:length(unique(MOCK$species)))
    
    p_mock <- alr(p_mock[,3:ncol(p_mock)]) %>% as.matrix() %>% as.data.frame()
    p_mock[,length(unique(MOCK$species))] <- 0  #add reference zero expressly
    names(p_mock)[length(unique(MOCK$species))] <- paste0("alr_", length(unique(MOCK$species)))
    
    p_mock <-  cbind(MOCK %>% dplyr::select(c("level1", "level2", "level3")[Nlevels], S) %>% distinct(),
                     p_mock) %>% ungroup()
    #names(p_mock)[1] <- "tech_rep" ##omit in favor of level-specific labeling? 
    
    
    return(p_mock)
}
  
  
  
makeDesign <- function(obs, #obs is a named list with elements Observation, Mock, N_pcr_mock, sp_list
                         N_pcr_cycles,
                         Nlevels_mock = 2, #levels of replication. Samples with either tech or biol replicates = 2; Samples with tech AND biol replicates = 3
                         Nlevels_samp = 2 #levels of replication. Samples with either tech or biol replicates = 2; Samples with tech AND biol replicates = 3
                         ){ #N_pcr_cycles is the number of PCR cycles in your experimental/enviro samples; currently a single value, could be made into a vector if this number varies
    #library(tidyverse)
    library(MCMCpack)
    library(compositions)
    library(rstan)
    library(dplyr)
    
    
    mock <- obs$Mock %>% 
      unite(c("level1", "level2", "level3")[1:Nlevels_mock-1], col = S, sep = "_", remove = F) #create identifier to distinguish the lowest level of replication (i.e., site-sample)
    observed <- obs$Observation %>% 
      arrange(species, Sample) %>% 
      unite(c("level1", "level2", "level3")[1:Nlevels_samp-1], col = S, sep = "_", remove = F) #create identifier to distinguish the lowest level of replication
    
    rep_level_mock <- c("level1", "level2", "level3")[Nlevels_mock]  #name of the column with highest level of replication
    rep_level_samp <- c("level1", "level2", "level3")[Nlevels_samp]  #name of the column with highest level of replication
    
    p_mock_all <- alrTransform(mock, Nlevels_mock)
    
    mock <- mock %>% 
      dplyr::select(species, 
                    S,  #unique biological samples
                    all_of(rep_level_mock),  #lowest level of replication
                    Nreads) %>% 
      ungroup() %>% 
      mutate(species = paste0("sp_", species)) %>% 
      pivot_wider(names_from = species, values_from = Nreads, values_fill = 0)
    N_pcr_mock <- rep(obs$N_pcr_mock, nrow(p_mock_all)) #assumes all have the same Npcr
    
    
    p_samp_all <- observed %>% 
      ungroup() %>% 
      dplyr::select(species, 
                    S,  #unique biological samples
                    all_of(rep_level_samp),  #lowest level of replication
                    Nreads) %>% 
      mutate(species = paste0("sp_", species)) %>% 
      #arrange(species) %>% 
      pivot_wider(names_from = species, values_from = Nreads, values_fill = 0)
    N_pcr_samp <- rep(N_pcr_cycles, nrow(p_samp_all))

    ########################################################################
    #### Create data frames that can be read into Stan model
    ########################################################################
    
    NOM <- as.name(colnames(p_mock_all)[1])
    formula_a <- eval(NOM) ~ N_pcr_mock -1
    model_frame <- model.frame(formula_a, p_mock_all)
    model_vector_a_mock <- model.matrix(formula_a, model_frame) %>% as.numeric()
    N_pcr_mock_small <- cbind(N_pcr_mock, p_mock_all) %>%  slice(match(unique(p_mock_all$S), p_mock_all$S)) %>% pull(N_pcr_mock)
    formula_b <- eval(NOM) ~ N_pcr_mock_small -1
    model_frame <- model.frame(formula_b, p_mock_all%>% slice(match(unique(p_mock_all$S), p_mock_all$S)))
    model_vector_a_mock_small <- model.matrix(formula_b, model_frame) %>% as.numeric()
    N_obs_mock       <- nrow(p_mock_all)
    
    # unknown communities second
    # species compositions (betas)
    
    NOM <- as.name(colnames(p_samp_all)[1])    
    p_samp_all$S <- as.factor(p_samp_all$S) 
    N_S = length(unique(p_samp_all$S))
    p_samp_all[rep_level_samp] <- as.factor(unlist(p_samp_all[rep_level_samp]))
    if(N_S == 1){
      formula_b <- eval(NOM) ~ 1  
    } else {
      formula_b <- eval(NOM) ~ 0+S
    }
    
    model_frame <- model.frame(formula_b, p_samp_all)
    model_matrix_b_samp <- model.matrix(formula_b, model_frame)
    
    # choose a single representative for each station to make predictions to
    model_frame <- model.frame(formula_b, p_samp_all[match(unique(p_samp_all$S), p_samp_all$S),])
    model_matrix_b_samp_small <- model.matrix(formula_b, model_frame)
    
    # efficiencies (alpha)
    formula_a <- eval(NOM) ~ N_pcr_samp -1
    model_frame <- model.frame(formula_a, p_samp_all)
    model_vector_a_samp <- model.matrix(formula_a, model_frame) %>% as.numeric()
    N_pcr_samp_small <- cbind(N_pcr_samp, p_samp_all) %>% slice(match(unique(p_samp_all$S), p_samp_all$S)) %>% pull(N_pcr_samp)
    formula_b <- eval(NOM) ~ N_pcr_samp_small -1
    
    model_frame <- model.frame(formula_b, p_samp_all %>% slice(match(unique(p_samp_all$S), p_samp_all$S)))
    model_vector_a_samp_small <- model.matrix(formula_b, model_frame) %>% as.numeric()
    
    #counters 
    N_obs_samp_small <- nrow(model_matrix_b_samp_small)
    N_obs_samp <- nrow(p_samp_all)
    N_b_samp_col <- ncol(model_matrix_b_samp)
    
    
    #### Make Stan objects
    
    stan_data <- list(
      N_species = ncol(p_samp_all)-2,   # Number of species in data
      N_obs_mb_samp = nrow(p_samp_all), # Number of observed community samples and tech replicates ; this will be Ncreek * Nt * Nbiol * Ntech * 2 [for upstream/downstream observations]
      N_obs_mock = nrow(p_mock_all), # Number of observed mock samples, including tech replicates
      N_obs_mb_samp_small = nrow(p_samp_all[match(unique(p_samp_all$S), p_samp_all$S),]), # Number of unique observed community samples ; this will be Ncreek * Nt * Nbiol * 2 [for upstream/downstream observations]
      
      # Observed data of community matrices
      sample_data = p_samp_all %>% dplyr::select(contains("sp")),
      # sample_vector = p_samp_all$S,
      mock_data   = mock %>% dplyr::select(contains("sp")),
      # sp_list = obs$sp_list,
      
      # True proportions for mock community
      #mock_true_prop = p_mock_all %>% dplyr::select(contains("sp")),
      alr_mock_true_prop = p_mock_all %>% dplyr::select(contains("alr")),
      
      # vectors of PCR numbers
      # N_pcr_samp = N_pcr_samp,
      # N_pcr_mock = N_pcr_mock,
      
      # Design matrices: field samples
      N_b_samp_col = N_b_samp_col,
      model_matrix_b_samp = model_matrix_b_samp,
      model_matrix_b_samp_small = as.array(model_matrix_b_samp_small),
      model_vector_a_samp = model_vector_a_samp,
      model_vector_a_samp_small = as.array(model_vector_a_samp_small),
      
      # Design matrices: mock community samples
      model_vector_a_mock = as.array(model_vector_a_mock),
      
      # Priors
      alpha_prior = c(0,0.5),  # normal prior
      # beta_prior = c(0,5),    # normal prior
      tau_prior = c(1,2)   # gamma prior
    )
    
    return(stan_data)
    
}

  #example
  #stan_metabarcoding_data <- makeDesign(metabarcoding_data, N_pcr_cycles = 43)    
  
  
###########################################
### Setting Initial Values
stan_init_f1 <- function(n.chain,N_obs_mb,N_species,Nplates){
  set.seed(78345)
  A <- list()
  for(i in 1:n.chain){
    A[[i]] <- list(
      log_D_raw=matrix(data=rnorm(N_obs_mb*N_species,mean=5,sd=2),nrow = N_obs_mb,ncol=N_species),
      alpha_raw=jitter(rep(0,N_species-1),factor=0.5),
      beta_std_curve_0=runif(Nplates,35,45),
      beta_std_curve_1=runif(Nplates,-1.5,-1.2)
    )
  }  
  return(A)
}


