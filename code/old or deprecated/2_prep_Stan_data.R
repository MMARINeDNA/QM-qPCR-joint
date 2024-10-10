library(here)
library(tidyverse)

keep <- intersect(unique(mfu$BestTaxon), template_props$Species) %>% intersect(fish)
keep <- c(keep, "Zz_Engraulis mordax")
mfu$BestTaxon[mfu$BestTaxon == "Engraulis mordax"] <- "Zz_Engraulis mordax"
mfu$BestTaxon[mfu$BestTaxon == "Merluccius"] <- "Merluccius productus"

###NOTE, if calibrating, change BestTaxon to Species, and:
mfu <- mfu %>%
  rename(Species = BestTaxon) %>% 
  filter(Species %in% keep) #just take fish species in mock
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

# finally, data in Stan-input format (list)
# metabarcoding
stan_metabarcoding_data <- format_metabarcoding_data(mfu, mock,
                                                     Level_1_treatment_envir <- "Sample",
                                                     Level_2_treatment_envir <- "tech",
                                                     Level_3_treatment_envir <- NA,
                                                     Level_1_treatment_mock <- "Sample",
                                                     Level_2_treatment_mock <- "Rep",
                                                     Level_3_treatment_mock <- NA) %>% 
  makeDesign(N_pcr_cycles = 43, Nlevels_mock = 2, Nlevels_samp = 2)

#qPCR
stan_qPCR_data <- format_qPCR_data(qPCR_unk,qPCR_std) %>% 
  prepare_stan_data_qPCR()

stan_data_joint <- c(stan_qPCR_data,stan_metabarcoding_data)
