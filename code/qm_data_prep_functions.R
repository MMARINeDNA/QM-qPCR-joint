# Metabarcoding calibration model function
# RPK Aug 2022; Revised ORL July 2024

################DEFINE SUB-FUNCTIONS

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
format_qPCR_data <- function(qPCR_unknowns, qPCR_standards){
  
  #unknowns
  qPCR_1<- qPCR_unknowns %>% 
    # pick columns we care about and rename
    select(qPCR, well,tubeID,type,task,IPC_Ct,inhibition_rate,Ct=hake_Ct,copies_ul=hake_copies_ul) %>% 
    # remove completely empty rows
    filter(!is.na(qPCR)) %>%
    mutate(z=ifelse(Ct=="Undetermined",0,1)) %>%
    mutate(Ct=str_replace_all(Ct,"Undetermined",'')) %>% 
    mutate(Ct=as.numeric(Ct) %>% round(2)) %>% 
    filter(task=="UNKNOWN",type=="unknowns")
  
  #standards
  qPCR_2 <- qPCR_standards %>% 
    # pick columns we care about and rename
    select(qPCR, well,tubeID,type,task=hake_task,IPC_Ct,inhibition_rate,Ct=hake_Ct,copies_ul=hake_copies_ul) %>% 
    mutate(z=ifelse(Ct=="Undetermined",0,1)) %>%
    mutate(Ct=str_replace_all(Ct,"Undetermined",'')) %>% 
    mutate(Ct=as.numeric(Ct) %>% round(2))%>% 
    filter(task=="STANDARD") %>% 
  # HAD AN ISSUE WITH NON-UNIQUE SAMPLE NAMES BETWEEN STDS AND UNKS BECAUSE OF '5' BEING USED AS SAMPLE ID FOR STANDARDS WITH CONC. OF 5 COPIES
    mutate(tubeID=ifelse(tubeID=="5","5C",tubeID))
  
  # bind standards and unknowns
  qPCRdata <- qPCR_1 %>% 
    bind_rows(qPCR_2) %>% 
    mutate(Ct = ifelse(is.na(Ct), 99, Ct)) %>%  # Stan doesn't like NAs
    mutate(plate_idx=match(qPCR,unique(qPCR))) %>% 
    unite(c(qPCR,tubeID), col = "plateSample", remove = F) %>% 
    mutate(plateSample_idx = match(plateSample, unique(plateSample))) %>% 
    group_by(plateSample) %>% 
    add_tally(Ct==99) %>% 
    # filter(n < 3) %>% #do away with examples of three non-detections; we have no basis for modeling these.
    dplyr::select(-n) %>%
    ungroup() %>% 
    mutate(plateSample_idx = match(plateSample, unique(plateSample)))
  
  return(qPCRdata)
  
}

# a last piece is we need a sample identifier across QM and qPCR data
# We need to link unique qPCR biological samples to unique QM samples
prepare_stan_qPCR_mb_join <- function(input_metabarcoding_data,qPCR_unknowns, qPCR_standards,link_species="Merluccius productus"){

  qpcr_formatted <- format_qPCR_data(qPCR_unknowns,qPCR_standards)
  
  mb_link_1 <- input_metabarcoding_data %>% 
    mutate(mb_link=match(Sample,unique(Sample))) %>% 
    distinct(Sample,.keep_all = T)
  
  mb_link_2 <- qpcr_formatted %>% 
    distinct(plateSample,plateSample_idx,.keep_all = T) %>% 
    left_join(mb_link_1,by=join_by(tubeID==Sample)) %>% 
    dplyr::select(plateSample_idx,mb_link) %>% 
    arrange(plateSample_idx) %>% 
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
  
  type <- qPCRdata %>% distinct(plateSample, task) %>% pull(task)
  
  stan_qPCR_data <- list(
    Nplates = length(unique(qPCRdata$qPCR)),
    Nobs_qpcr = nrow(qPCRdata),
    NSamples_qpcr = length(unique(qPCRdata$plateSample)),
    NstdSamples = qPCRdata %>% filter(task == "STANDARD") %>% distinct(plateSample) %>% nrow(),
    plate_idx = qPCRdata$plate_idx, 
    std_idx =  which(type=="STANDARD"),
    unkn_idx = which(type == "UNKNOWN"),
    plateSample_idx = qPCRdata$plateSample_idx,
    y = qPCRdata$Ct,
    z = qPCRdata$z,
    known_concentration = qPCRdata %>% filter(task == "STANDARD") %>% distinct(plateSample,.keep_all=T) %>% pull(copies_ul),
    stdCurvePrior_intercept = c(39, 3), #normal distr, mean and sd ; hyperpriors
    stdCurvePrior_slope = c(-3, 1) #normal distr, mean and sd ; hyperpriors
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
stan_init_f1 <- function(n.chain,N_obs_mb,N_species){
  set.seed(78345)
  A <- list()
  for(i in 1:n.chain){
    A[[i]] <- list(
      log_D_raw=matrix(data=rnorm(N_obs_mb*N_species,mean=5,sd=2),nrow = N_obs_mb,ncol=N_species),
      alpha_raw=jitter(rep(0,N_species-1),factor=0.5)
    )
  }  
  return(A)
}


