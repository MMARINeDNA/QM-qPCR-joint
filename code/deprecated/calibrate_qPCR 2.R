# qPCR calibration model function
# RPK Aug 2022
# added EA Aug 2022

################DEFINE SUB-FUNCTIONS

    input_qPCR_data <- function(datafilepath, metafilepath, plate_num) {
      #function to take .csv file straight from instrument and make .csv file using only final
      #input = directory of files from instrument, number within files, and metadata saying which samples to use from which plate
      #output = csv file with final samples, Ct values for each marker
      #dependencies = listed packages
      
      require(tidyverse)
      require(readxl)
      # first just get file open - they are listed by date so plate number is not in order
        qPCRfiles <- list.files(path = datafilepath, pattern = "*data", recursive = T, full.names = T)
        qPCRmeta <- read_excel(path=metafilepath, sheet = "info_samples")
        qPCRmeta$use <- as.numeric(qPCRmeta$use)
        adj_vol_filtered <- readRDS(here("Output","qpcr","backpack","adj_vol_filtered.RDS"))
        
        whichmarker <- str_detect(qPCRfiles[1], "CUT")
        marker <- ifelse(whichmarker, "cutt", "coho")
        fn3 <- ifelse(whichmarker, 3,2)
        
        filename = qPCRfiles[plate_num]
        filename = gsub(".+quantitative","",filename) #RPK added to deal w changing paths on different computers
        filename2 = unlist(strsplit(filename, "_"))[4]
        filename3 = unlist(strsplit(filename2,"-"))[fn3]
        filename4 = str_sub(filename3, start=1, end =3)
        plate <- read_excel(qPCRfiles[plate_num], sheet = "Results", skip=6)
        colnames(plate)[7] <- "Ct"  
        colnames(plate)[2] <- "Sample"
        colnames(plate)[3] <- "Assay"
        plate$Plate <- as.numeric(filename4) # now right the actual plate number that it is (i.e., what Megan called it)
        
        # now only look at assay data and Ct values 
        sampdata <- plate %>% 
          filter(Assay != "IPC") %>% 
          dplyr::select(c("Plate", "Well","Sample", "Ct")) 
        
        # now use metadata file to only select samples that are not inhibited (many samples were run 2-3x)
        touse <- qPCRmeta %>% 
          filter(qPCR_no == filename2) %>% 
          filter(target == marker) %>%  # all the inhibition is on cutthroat plates
          filter(use == 1) %>% # sample by sample, plate by plate (only use good SC plates, and only use the samples we THOUGHT were not inhibited using the 3 ct threshold)
          filter(! str_detect(sample, "NTC")) %>% 
          filter(! str_detect(sample, "EB")) 
        
        dilution <- touse %>% 
          dplyr::select(c(sample,dilution)) %>% 
          rename(Sample = sample) 
        
        sampdata <- sampdata %>% 
          filter(Sample %in% touse$sample) %>% 
          mutate(Ct = na_if(Ct,"Undetermined")) %>% 
          mutate(Ct = as.numeric(Ct)) %>% 
          left_join(dilution, by="Sample") %>% 
          left_join(adj_vol_filtered, by = "Sample")
        return(sampdata)
    }
    
    # #example:
    # library(here)
    # qPCRfiles <- here("Input","qpcr","Results","COHO")
    # qPCRmeta <- here("Input","qpcr","qPCR_samples_ALL.xlsx")
    # test <- input_qPCR_data(qPCRfiles,qPCRmeta, 1)
    
format_qPCR_data <- function(input_csv){
  #function to format input data for stan-code calibration
  #input = data frame with required columns
  #output = named list
  #dependencies = listed packages
  
  require(tidyverse)

        qPCRdata <- read.csv(input_csv) 
        
        
        qPCRdata <- qPCRdata %>% 
          mutate(Task=ifelse(str_detect(Sample,"S"), "STANDARD", "UNKNOWN")) %>% 
          mutate(z = ifelse(is.na(Ct), 0, 1),
                 Ct = ifelse(is.na(Ct), 99, Ct)) %>% 
          mutate(Ct = as.numeric(Ct),
                 Ct = round(Ct, 2)) %>%  #do away with ridiculous decimal places 
          mutate(conc = NA) %>% 
          mutate(Quantity = 0) %>% 
          mutate(., Quantity = case_when(Sample == "St1" ~ 100000,
                                         Sample == "St2" ~ 10000,
                                         Sample == "St3" ~ 1000,
                                         Sample == "St4" ~ 100,
                                         Sample == "St5" ~ 10,
                                         Sample == "St6" ~ 5,
                                         Sample == "St7" ~ 3,
                                         Sample == "St8" ~ 1,
                                         TRUE ~ Quantity)) %>% 
          mutate(PlateName = Plate) %>% #retain original plate name
          mutate(Plate = match(Plate, unique(Plate))) #index plates sequentially
        
        
        qPCRdata$conc[qPCRdata$Task == "STANDARD"] <- qPCRdata$Quantity[qPCRdata$Task == "STANDARD"] %>% as.numeric()
        qPCRdata$Sample[qPCRdata$Task == "STANDARD"] <- qPCRdata$Quantity[qPCRdata$Task == "STANDARD"] 
        
        qPCRdata <- qPCRdata %>% 
          unite(c(Plate,Sample), col = "plateSample", remove = F) %>% 
          mutate(plateSample_idx = match(plateSample, unique(plateSample))) %>% 
          group_by(plateSample) %>% 
          add_tally(Ct == 99) %>% 
          filter(n < 3) %>% #do away with examples of three non-detections; we have no basis for modeling these. 
          dplyr::select(-n) %>% 
          mutate(plateSample_idx = match(plateSample, unique(plateSample))) #reindex
        
        #QC filter
        qPCRdata <- qPCRdata %>%
          group_by(Plate, plateSample) %>% 
          mutate(sdCt = sd(Ct)) %>% 
          filter(sdCt < 1.5 | sdCt > 20) %>% #really bad triplicates; don't throw out instances where 1 or 2 out of 3 was a detection (which, because Ct = 99 for nondetect, inflates the sd)
          ungroup() %>% 
          mutate(plateSample_idx = match(plateSample, unique(plateSample))) #reindex
  
  
        return(qPCRdata)
  
}

  #example:
  # library(here)
  # q <- format_qPCR_data(here("Output","qpcr","coho_final.csv"))



prepare_stan_data_qPCR <- function(qPCRdataname){
  
  require(tidyverse)
  
  qPCRdata <- qPCRdataname
  
  type <- qPCRdata %>% dplyr::select(plateSample, Task) %>% distinct() %>% pull(Task)
  
  stan_qPCR_data <- list(
    Nplates = length(unique(qPCRdata$Plate)),
    Nobs = nrow(qPCRdata),
    NSamples = length(unique(qPCRdata$plateSample)),
    NstdSamples = qPCRdata %>% filter(Task == "STANDARD") %>% dplyr::select(plateSample) %>% distinct() %>% nrow(),
    plate_idx = qPCRdata$Plate, 
    std_idx =  which(type=="STANDARD"),
    unkn_idx = which(type != "STANDARD"),
    plateSample_idx = qPCRdata$plateSample_idx, 
    y = qPCRdata$Ct,
    z = qPCRdata$z,
    known_concentration = qPCRdata %>% filter(Task == "STANDARD") %>% dplyr::select(plateSample_idx, conc) %>% distinct() %>% pull(conc),
    stdCurvePrior_intercept = c(39, 3), #normal distr, mean and sd ; hyperpriors
    stdCurvePrior_slope = c(-3, 1) #normal distr, mean and sd ; hyperpriors
  )
  
  return(stan_qPCR_data)
  
}
    
    #example
    # stan_qPCR_data <- prepare_stan_data_qPCR(q)


calibrate_qPCR <- function(stanmodelname, dataname, NCHAINS = 3, WARMUP = 500, ITER = 2500){
  #actual calibration in stan
  #input = named list
  #output = posterior samples for many params
  #dependencies = listed packages
  
  require(tidyverse)
  require(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  
  
  qMod = stan(file = stanmodelname, data = dataname,
              verbose = FALSE, chains = NCHAINS, thin = 1,
              warmup = WARMUP, iter = ITER,
              control = list(adapt_init_buffer = 175,
                             max_treedepth=12,
                             stepsize=0.01,
                             adapt_delta=0.7,
                             metric="diag_e"),
              refresh = 10,
              boost_lib = NULL
  )
  
  return(qMod)
  
}

#example:
# library(here)
# qMod <- calibrate_qPCR(stanmodelname = here("Scripts", "qm_qpcr_model", "qPCR_calibration_enchilada.stan"),
#                        dataname = stan_qPCR_data)


########################      

#all-in-one function that calls those above:
#inputs == input_csv, stanmodelname
#output == fitted stan model object

run_qPCR_model <- function(input_csv_name, stanmodel){
  require(dplyr)
  
  qMod <- input_csv_name %>% 
    format_qPCR_data() %>% 
    prepare_stan_data_qPCR() %>% 
    calibrate_qPCR(stanmodel, dataname = .)

  qPCRdata <- input_csv_name %>% 
    format_qPCR_data()
  
  results_qPCR <- qPCRdata %>% 
    filter(Task != "STANDARD") %>% 
    separate(Sample, into = c("time","creek","station","biorep"), remove = F) %>% 
    dplyr::select(Plate, time, creek, station, biorep, dilution, Adj_Vol) %>% 
    filter(creek != "EB") %>% 
    distinct() %>% 
        mutate(mean_concentration_est = dilution*(100/(Adj_Vol/2))*10^(summary(qMod, par = "envir_concentration")$summary[,1]),
               ci25_concentration_est = dilution*(100/(Adj_Vol/2))*10^(summary(qMod, par = "envir_concentration")$summary[,5]),
               ci75_concentration_est = dilution*(100/(Adj_Vol/2))*10^(summary(qMod, par = "envir_concentration")$summary[,7])) %>% 
    left_join(qPCRdata %>% dplyr::select(time, creek, station, biorep, Ct))
  
  
  return(
    list("qMod" = qMod, 
         "results_qPCR" = results_qPCR)
  )
  
}
      
      ## example
      # library(here)
      # qMod_out <- run_qPCR_model(here("Output/qpcr/cut_data_for_stan.csv"),
      #                here("Scripts/functions/qPCR_calibration_enchilada.stan")
      #                )
      # qMod_out$results_qPCR %>% 
      #   filter(creek == "4Pad") %>% 
      #   ggplot(aes(x = time, y = log(mean_concentration_est), color = station)) +
      #     geom_point() 
          
      
      