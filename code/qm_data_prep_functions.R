# Metabarcoding calibration model function
# RPK Aug 2022; Revised March 2023

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
  #require(rlang)

  # level1 <- rlang::enquo(Level_1_treatment) #e.g., unique sampling site
  # level2 <- rlang::enquo(Level_2_treatment) #nested within level1 units, e.g., unique biological samples
  # level3 <- rlang::enquo(Level_3_treatment) #nested within level2 units, e.g., technical replicates of biological samples
  
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
    select(qPCR, well,type,task,IPC_Ct,inhibition_rate,Ct=hake_Ct,copies_ul=hake_copies_ul) %>% 
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
  
  #QC filter (old?)
  # qPCRdata <- qPCRdata %>%
  #   group_by(Plate, plateSample) %>% 
  #   mutate(sdCt = sd(Ct)) %>% 
  #   filter(sdCt < 1.5 | sdCt > 20) %>% #really bad triplicates; don't throw out instances where 1 or 2 out of 3 was a detection (which, because Ct = 99 for nondetect, inflates the sd)
  #   ungroup() %>% 
  #   mutate(plateSample_idx = match(plateSample, unique(plateSample))) #reindex
  
  
  return(qPCRdata)
  
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

# a last piece is we need a sample identifier across MB and qPCR data
prepare_qPCR_mb_join <- function(input_metabarcoding_data,qPCR_unknowns){
  
  mb_samps<-unique(input_metabarcoding_data$Sample)
  qPCR_samps <- unique(qPCR_unknowns$tubeID)
  inner <- intersect(mb_samps,qPCR_samps)
  
  qPCR_mb_link <- which(qPCR_unknowns$tubeID%in%mb_samps)
    
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
      formula_b <- eval(NOM) ~ S
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
###########################################
QM_likelihood <- function(stanmodelname, stan_metabarcoding_data){
    M <- stan_model(stanmodelname)
    
    stanOpt <- optimizing(M, data=stan_metabarcoding_data, iter=30000,draws=0,
                          verbose=T,  
                          tol_param=1e-40,
                          algorithm="LBFGS",
                          hessian = TRUE)
    
    MLest <- stanOpt$par[grep("int_samp_small", names(stanOpt$par))] %>%
      matrix(ncol = stan_metabarcoding_data$N_species) %>% 
      as.data.frame()
      names(MLest) <- stan_metabarcoding_data$sp_list$species
      rownames(MLest) <- unique(stan_metabarcoding_data$sample_vector) %>% sort()
    ML_a <- stanOpt$par[grep("alpha\\[", names(stanOpt$par))]  
    ML_a <- data.frame("alpha_est" = ML_a, 
               "species" = stan_metabarcoding_data$sp_list$species)
    
    
    return(list(
      ML_modelfit = stanOpt,
      ML_estimates = MLest,
      ML_alpha_est = ML_a
    ))
    
  }
  
  #example; note this can fail stochastically; run it several times if need be
  # ML_out <- QM_likelihood(here("quant_metabar_rosetta_noSampleEta.stan"), stan_metabarcoding_data)
###########################################
###########################################
  
###########################################
###########################################
  
QM_bayes <- function(stanmodelname, stan_metabarcoding_data, NCHAINS = 3, WARMUP = 500, ITER = 1500){
  require(tidyverse)
  require(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  
  # stan_pars <- c( 
  #   "alpha",
  #   "beta",
  #   "eta_mock",
  #   "tau",
  #   "mu_samp",
  #   "mu_mock",
  #   "int_samp_small"
  # )
  
    stanMod = stan(file = stanmodelname ,data = stan_metabarcoding_data,
                   verbose = FALSE, chains = NCHAINS, thin = 1,
                   warmup = WARMUP, iter = ITER,
                   control = list(adapt_init_buffer = 175,
                                  max_treedepth=12,
                                  stepsize=0.01,
                                  adapt_delta=0.7,
                                  metric="diag_e"),
                   # pars = stan_pars,
                   refresh = 10,
                   boost_lib = NULL 
    )  
    
    
    mean_est <- summary(stanMod, par = "int_samp_small")$summary[,1] %>%
      matrix(ncol = stan_metabarcoding_data$N_species, byrow = TRUE) %>% 
      as.data.frame()
      names(mean_est) <- stan_metabarcoding_data$sp_list$species
      rownames(mean_est) <- unique(stan_metabarcoding_data$sample_vector) %>% sort()
    ci25_est <- summary(stanMod, par = "int_samp_small")$summary[,5] %>%
      matrix(ncol = stan_metabarcoding_data$N_species, byrow = TRUE) %>% 
      as.data.frame()
      names(ci25_est) <- stan_metabarcoding_data$sp_list$species
      rownames(ci25_est) <- unique(stan_metabarcoding_data$sample_vector) %>% sort()
    
    ci75_est <- summary(stanMod, par = "int_samp_small")$summary[,7] %>%
      matrix(ncol = stan_metabarcoding_data$N_species, byrow = TRUE) %>% 
      as.data.frame()
      names(ci75_est) <- stan_metabarcoding_data$sp_list$species
      rownames(ci75_est) <- unique(stan_metabarcoding_data$sample_vector) %>% sort()
    
    mean_a_est <- summary(stanMod, par = "alpha")$summary[,1]
    mean_a_est <- data.frame("alpha_est" = mean_a_est, 
                       "species" = stan_metabarcoding_data$sp_list$species)
    
    #unique(meta.samples$sample)
    
    return(list(
      Bayes_modelfit = stanMod,
      Bayes_estimates = mean_est,
      Bayes_25ci = ci25_est,
      Bayes_75ci = ci75_est,
      Bayes_alpha_est = mean_a_est
    ))
    
  }
  
##example
  #QM_bayes_out <- QM_bayes(here("quant_metabar_rosetta_noSampleEta.stan"), stan_metabarcoding_data)
  #QM_bayes_out$Bayes_estimates
  ###########################################
  ###########################################
  
  
  
  ########################      
  
  #all-in-one function that calls those above:
  #inputs == input_csv, stanmodelname
  #output == fitted stan model object
  
  run_QM_model <- function(input_metabarcoding_RDS, 
                                      input_mock_comm_RDS, 
                                      N_pcr_cycles = 43,
                                      stanmodel, 
                                      method = "ML"  # "ML" or "Bayes"
                                      ){
    require(dplyr)
    require(here)
    
    
    metabarcoding_data <- format_metabarcoding_data(input_metabarcoding_RDS,
                                                    input_mock_comm_RDS)
    
    stan_metabarcoding_data <- makeDesign(metabarcoding_data, N_pcr_cycles = N_pcr_cycles) 
    
    if (method == "ML"){
      ML_out <- QM_likelihood(stanmodel, stan_metabarcoding_data)
    } 
    
    if (method == "Bayes"){
      QM_bayes_out <- QM_bayes(stanmodel, stan_metabarcoding_data)
    } 
    
    if (!method %in% c("ML", "Bayes")) {
      print("Pick a valid method: `ML' or `Bayes'")
    }
    
    
    if (method == "ML"){
      ML_out <- QM_likelihood(stanmodel, stan_metabarcoding_data)
    } 
    
    if (method == "Bayes"){
      QM_bayes_out <- QM_bayes(stanmodel, stan_metabarcoding_data)
    } 
    
    if (!method %in% c("ML", "Bayes")) {
      print("Pick a valid method: `ML' or `Bayes'")
    }
    
    
    if (method == "ML"){
      fit_out <- list(metabarcoding_data = metabarcoding_data,
                      ML_out = ML_out)
    } 
    
    if (method == "Bayes"){
      fit_out <- list(metabarcoding_data = metabarcoding_data,
                      QM_bayes_out = QM_bayes_out)
    }
    
    return(fit_out)
  }
  
  ## ML example
  # library(here)
  # ML_out <- run_QM_model(input_metabarcoding_RDS = here("../example_data/padden_sample.rds"),
  #                        input_mock_comm_RDS = here("/Volumes/GoogleDrive/My Drive/RPKDesktop/github_repos/quantitative_salmon_culverts/Input/mockcommunity/MiFishmockdata_noZ.RDS"),
  #                        stanmodel = here("quant_metabar_rosetta_noSampleEta.stan"),
  #                        method = "ML"
  #                )
  # # ML_out$ML_out$ML_estimates %>% 
  #   rownames_to_column("sample") %>% 
  #   pivot_longer(-sample, names_to = "species") %>% 
  #   separate(col = sample, into = c("time", "creek", "station", "biol"), remove = FALSE) %>% 
  #   filter(value > 0.001) %>% 
  #   ggplot(aes(x = sample, fill = species, y = value)) +
  #     geom_col() +
  #     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      
  ##bayesian example
  # library(here)
  # Bayes_out <- run_QM_model(input_metabarcoding_RDS = here("../example_data/padden_sample.rds"),
  #                           input_mock_comm_RDS = here("/Volumes/GoogleDrive/My Drive/RPKDesktop/github_repos/quantitative_salmon_culverts/Input/mockcommunity/MiFishmockdata_noZ.RDS"),
  #                           stanmodel = here("quant_metabar_rosetta_noSampleEta.stan"),
  #                           method = "Bayes"
  # )
  # Bayes_out$QM_bayes_out$Bayes_estimates %>%
  #   rownames_to_column("sample") %>%
  #   pivot_longer(-sample, names_to = "species") %>%
  #   separate(col = sample, into = c("time", "creek", "station", "biol"), remove = FALSE) %>%
  #   filter(value > 0.001) %>%
  #   ggplot(aes(x = sample, fill = species, y = value)) +
  #   geom_col() +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  # saveRDS(Bayes_out, "mod_fit_aug21.rds")
  # 
  # 
  # 
  # 
  # 
