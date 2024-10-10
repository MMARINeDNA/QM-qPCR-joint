#Load packages and data for MURI Amp Bias Analysis

####Packages and Functions
library(here)
suppressMessages(library(tidyverse))
  options(tidyverse.quiet = TRUE)
  options(dplyr.summarise.inform = FALSE)
  options(dplyr.left_join.inform = FALSE)
suppressMessages(library(rstan))
  options(mc.cores = parallel::detectCores())
  rstan_options(threads_per_chain = 4)
  rstan_options(auto_write = TRUE)
suppressMessages(library(compositions))
suppressMessages(library(MCMCpack))


logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
} #from https://gist.github.com/aufrank/83572

select <- dplyr::select
source(here("code/calibrate_metabarcoding.R"))
      
####Data
      
#hake sample metadata
meta <- read.csv(here("data/metadata/Hake_2019_metadata.csv"))

#species and all assigned metabarcoding data
sp <- read.table(here("data/mocks/species_list.txt"), sep = "\t") %>% 
  rename("Species" = "V1")
q <- read.csv(here("data/mocks/mock_paper_ASV_assignments.csv"), row.names = 1) %>% 
  filter(Species %in% sp$Species) %>% 
  filter(!is.na(Species)) %>% 
  select(Sample, Reads, Marker, Species) %>% 
  mutate(BSA = ifelse(str_detect(Sample, "BSA"), 1, 0),
         TD = ifelse(str_detect(Sample, "TD"), 1, 0),
         SKEW = ifelse(str_detect(Sample, "SKEW"), 1, 0)) %>% 
  mutate(Sample = str_replace_all(Sample, "BSA\\.", ""),
         Sample = str_replace_all(Sample, "TD\\.", ""),
         Sample = str_replace_all(Sample, "SKEW\\.", "")) %>% 
  separate(Sample, into = c("Marker", "Taq", "rep")) %>% 
  group_by(Marker, Taq, rep, Species, BSA, TD, SKEW) %>% 
  summarise(Reads = sum(Reads)) %>%  #collapse ASVs into taxa
  ungroup()
q$Species[which(q$Species == "Engraulis mordax")] <- "Zz_Engraulis mordax"

#we can't estimate species that are only rarely observed, so we are going to drop those. 
allspecies <- unique(q$Species)
lose <- q %>% group_by(Species) %>% summarise(s = sum(Reads)) %>% arrange(s) %>% head(2) %>% pull(Species)
keep <- allspecies[-which(allspecies %in% lose)]
# keep <- c(keep, "Zz_Engraulis mordax") #to make Engraulis the reference species

q <- q %>% 
  select(-any_of(lose))

mismatch <- read.csv(here("data/mocks/mock_paper_mismatches_v3.csv"))

fish <- mismatch$Species[mismatch$Class=="Actinopteri"]
fish[which(fish == "Engraulis mordax")] <- "Zz_Engraulis mordax"

ddPCR <- read.csv(here("data/mocks/mock_paper_ddPCR_v3.csv")) %>% 
  #filter(Species %in% keep) %>% 
  filter(!is.na(Conc_copies_ul)) %>% #quick/dirty; omit dropouts
  group_by(Species, Marker, Qubit_0.05) %>% 
  summarise(conc = mean(Conc_copies_ul)) %>% 
  ungroup() %>% 
  rowwise(conc = conc/Qubit_0.05) %>%  #correct for different templates
  ungroup()

#read counts from mocks
mock_reads <- read.csv(here("data/mocks/mock_paper_df_v2.csv"))
  mock_reads$Species[which(mock_reads$Species == "Engraulis mordax")] <- "Zz_Engraulis mordax" #set as ref species
#mock community gDNA proportions, by group and by species
mock_props <- mismatch %>% 
  select(Species, Class) %>% 
  left_join(tibble(Class = c("Chondrichthyes", "Actinopteri","Mammalia", "Cephalopoda", "Malacostraca"),
                   Group = c("Fish", "Fish", "Mammal", "Ceph", "Fish"),
                   GroupProportion = c(.65, .65, .05, .3, .65))) %>% 
  group_by(Group) %>% 
  mutate(Species_in_Group = n()) %>% 
  ungroup() %>% 
  mutate(SpeciesProportion = GroupProportion/Species_in_Group) %>% 
  select(-c(Species_in_Group, GroupProportion)) 
# filter(Group != "Krill") %>%  #renormalize species proportions without krill, if desired
# mutate(SpeciesProportion = SpeciesProportion/sum(SpeciesProportion))

#calculate template (mtDNA) proportions, given ddPCR data
template_props <- ddPCR %>% 
  left_join(mock_props) %>% 
  select(Species, conc, Marker, SpeciesProportion) %>% 
  group_by(Marker) %>% 
  mutate(genomic_prop = SpeciesProportion/sum(SpeciesProportion)) %>% 
  mutate(template_prop = (conc*SpeciesProportion)/sum(conc*SpeciesProportion)) %>% 
  select(-c(SpeciesProportion, conc))
template_props$Species[which(template_props$Species == "Engraulis mordax")] <- "Zz_Engraulis mordax" #set as ref species

## Annotation Database; NOTE other users will have a different shortcut to this file, but should have access through MURI Project_2022-2027 on Google Drive
db <- read.csv(here('data','MFU_database.csv'), row.names = 1)
      
## Metabarcoding Read Data from Environmental Samples
# MiFish
mfu1 <- read.csv(here("data/fish/MURI304_MFU_ASV_table.csv")) %>% 
  left_join(db) %>%  
  filter(!is.na(BestTaxon)) %>% 
  separate(Sample_name, into = c("Primer", "Project", "sample", "Dilution", "Well"), sep = c("-|_"), remove = F) %>% 
  mutate(Rep = "1") %>% 
  relocate(c("Primer", "Project", "sample", "Dilution", "Rep", "Well"))
# mfu2 <- read.csv(here("data/fish/MURI305_MFU_ASV_table.csv")) %>% 
#   filter(str_detect(Sample_name, "MFU-[0-9]")) %>% 
#   left_join(db) %>%  
#   filter(!is.na(BestTaxon))
mfu3 <- read.csv(here("data/fish/MURI313_MFU_ASV_table.csv")) %>% 
  left_join(db) %>%  
  filter(!is.na(BestTaxon)) %>% 
  separate(Sample_name, into = c("Primer", "Project", "sample", "Rep", "Dilution", "Well"), sep = c("-|_"), remove = F) %>% 
  relocate(c("Primer", "Project", "sample", "Dilution", "Rep", "Well"))
mfu4 <- read.csv(here("data/fish/MURI314_MFU_ASV_table.csv")) %>% 
  left_join(db) %>%  
  filter(!is.na(BestTaxon)) %>% 
  separate(Sample_name, into = c("Primer", "Project", "sample", "Dilution", "Rep", "Well"), sep = c("-|_"), remove = F)
mfu5 <- read.csv(here("data/fish/MURI315_MFU_ASV_table.csv")) %>% 
  left_join(db) %>%  
  filter(!is.na(BestTaxon)) %>% 
  separate(Sample_name, into = c("Primer", "Project", "sample", "Dilution", "Rep", "Well"), sep = c("-|_"), remove = F)
mfu6 <- read.csv(here("data/fish/MURI316_MFU_ASV_table.csv")) %>% 
  left_join(db) %>%  
  filter(!is.na(BestTaxon)) %>% 
  separate(Sample_name, into = c("Primer", "Project", "sample", "Dilution", "Rep", "Well"), sep = c("-|_"), remove = F)
mfu7 <- read.csv(here("data/fish/MURI317_MFU_ASV_table.csv")) %>% 
  left_join(db) %>%  
  filter(!is.na(BestTaxon)) %>% 
  separate(Sample_name, into = c("Primer", "Project", "sample", "Dilution", "Rep", "Well"), sep = c("-|_"), remove = F)
mfu8 <- read.csv(here("data/fish/MURI318_MFU_ASV_table.csv")) %>% 
  left_join(db) %>%  
  filter(!is.na(BestTaxon)) %>% 
  separate(Sample_name, into = c("Primer", "Project", "sample", "Dilution", "Rep", "Well"), sep = c("-|_"), remove = F)

      
mfu <- mfu1 %>% 
  bind_rows(mfu3,mfu4,mfu5,mfu6,mfu7,mfu8) %>% 
  group_by(Primer, Project, sample, Dilution, Rep, BestTaxon) %>% 
  summarise(nReads = sum(nReads)) %>% 
  pivot_wider(id_cols = c("Primer","Project","sample","Dilution","Rep"), names_from = "BestTaxon", values_from = "nReads", values_fill = 0) %>% 
  pivot_longer(cols = -c("Primer","Project","sample","Dilution","Rep"), names_to = "BestTaxon", values_to = "nReads")

mfu <- mfu %>% 
  ungroup() %>% 
  filter(!grepl("positive",Project)) %>%
  filter(!grepl("NTC",Project)) %>%
  mutate(sample = as.numeric(sample)) %>% 
  # separate(Sample_name, into = c("Primer", "Project", "sample", "Dilution", "Rep", "Well"), sep = c("-|_")) %>% 
  left_join(meta %>% select(sample, station, depth, lat, lon))
      #note: some of the dilution - rep columns seem to be switched. hence the merge will complain; need to sort this out.  