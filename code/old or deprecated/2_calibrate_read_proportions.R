library(here)
library(tidyverse)
source(here("code/1_LoadPackagesData.R"))

keep <- intersect(unique(mfu$BestTaxon), template_props$Species) %>% intersect(fish)
keep <- c(keep, "Zz_Engraulis mordax")
mfu$BestTaxon[mfu$BestTaxon == "Engraulis mordax"] <- "Zz_Engraulis mordax"
mfu$BestTaxon[mfu$BestTaxon == "Merluccius"] <- "Merluccius productus"

###NOTE, if calibrating, change BestTaxon to Species, and:
mfu <- mfu %>%
  rename(Species = BestTaxon) %>% 
  filter(Species %in% keep) #just take fish species in mock
#mfu$Species[which(mfu$Species == "Engraulis mordax")] <- "Zz_Engraulis mordax" #set as ref species
# mfu <- mfu %>%
#   filter(Species %in% template_props$Species) #filter for species for which we have mock data
# mfu <- mfu[,-which(colSums(mfu[,2:ncol(mfu)])==0)] #filter now-empty samples
mfu <- mfu %>%
  select(sample, Species, nReads) %>%
  drop_na() %>%
  # unite(station, depth, col = "station") %>%
  group_by(Species, sample) %>%
  summarise(nReads = sum(nReads)) %>% #NOTE will need to fix this; sums across dilutions, replicates, etc.
  pivot_wider(names_from = sample, values_from = nReads, values_fill = 0)

#just keep species in mocks that are also in the environmental data
template_props <- template_props %>% 
  filter(Species %in% keep)

#Fish
mock <- mock_reads %>% 
  filter(Primer == "MFU" & 
           Species %in% template_props$Species &
           Taq == "NPHF" &
           BSA == 1 &
           TD == 0) %>% 
  select(-Prop) %>% 
  mutate(Sample = substr(Sample, 1, 8))
true_prop <- template_props %>% 
  filter(Marker == "MV1" & Species %in% mock$Species) %>%  #note: using MV1 as true proportion, given lack of mismatches, as in Megan's paper
  mutate(template_prop = template_prop / sum(template_prop))
mock <- mock %>% 
  left_join(true_prop %>% select(Species, template_prop)) %>% 
  rename(b_proportion = template_prop,
         Nreads = Reads,
         species = Species) %>% 
  arrange(species)
  
  
  
mfu <- mfu %>% 
  pivot_longer(cols = -Species, names_to = "Sample", values_to = "Nreads") %>% 
  # filter(!grepl("positive",Sample)) %>% 
  # filter(!grepl("NTC",Sample)) %>% 
  # separate(Sample, into = c("Primer", "Project", "Sample", "Dilution"), sep = "\\.") %>% 
  # group_by(Sample, Species) %>% 
  # summarise(Nreads = sum(Nreads)) %>%  #sum across dilutions
  group_by(Sample) %>%
  mutate(totalReads = sum(Nreads)) %>% 
  filter(totalReads > 1000) %>% #impose total read-depth min
  ungroup() %>% 
  select(-totalReads) %>% 
  mutate(biol = 1, tech = 1) %>% 
  rename(species = Species) %>% 
  mutate(station = match(Sample, unique(Sample))) #let station = unique sampling location, treating different depths in the same x-y plane as independent


##The below bit does the actual calibration, both to est alpha and to est corrected proportions in all samples

stan_metabarcoding_data <- format_metabarcoding_data(mfu, mock,
                                                  Level_1_treatment_envir <- "Sample",
                                                  Level_2_treatment_envir <- "tech",
                                                  Level_3_treatment_envir <- NA,
                                                  Level_1_treatment_mock <- "Sample",
                                                  Level_2_treatment_mock <- "Rep",
                                                  Level_3_treatment_mock <- NA) %>% 
                           makeDesign(N_pcr_cycles = 43, Nlevels_mock = 2, Nlevels_samp = 2)
  
QM_bayes_out <- QM_bayes(here("code/quant_metabar_no_overdispersion.stan"),
                         stan_metabarcoding_data,
                         WARMUP = 500,
                         ITER = 2000)
saveRDS(QM_bayes_out, file = here("model_output/QM_bayes_MFU.RDS"))
#QM_bayes_out <- readRDS(here("model_output/QM_bayes_MFU.RDS"))
#QM_bayes_out$Bayes_estimates
QM_bayes_out$Bayes_alpha_est
#summary(QM_bayes_out$Bayes_modelfit)$summary[,"Rhat"] %>% sort() %>% tail()
rstan::traceplot(QM_bayes_out$Bayes_modelfit, par = "alpha")
      
      # head(QM_bayes_out$Bayes_estimates)
      
      
  
      
      
      