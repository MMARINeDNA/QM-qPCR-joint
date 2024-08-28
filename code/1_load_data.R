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
      
### Data
#hake sample metadata 

meta <- read_csv(here("data","metadata","Hake_2019_metadata.csv")) %>% 
  # keep only the relevant and unique columns from this data 
  select(station,Niskin,year,month,day,transect,lat,lon,utm.lat,utm.lon,bottom.depth.consensus,transect.dist.km) %>%
  distinct()

qPCR.sample.id <- read_csv(here('Data','hake_qPCR','Hake eDNA 2019 qPCR results 2023-02-10 sample details.csv'),
                           col_select = all_of(c("Tube #", "CTD cast","Niskin","depth","drop.sample","field.negative.type","water.filtered.L")))
# dat.station.id <- read_csv(here('Data','CTD_hake_data_10-2019.csv'))

### SAMPLE IDs ###
qPCR.sample.id <- qPCR.sample.id %>%  
  rename(tubeID=`Tube #`,
         station=`CTD cast`,
         volume=water.filtered.L) %>%
  distinct() %>% 
  mutate(depth=ifelse(Niskin=="sfc","0",depth)) %>%
  mutate(depth=ifelse(depth=="sfc","0",depth)) %>%
  mutate(depth=as.numeric(depth)) %>% 
  # empty stations or extraneous rows
  filter(!(station=="N/A"|station=="-")) %>%
  left_join(meta,by = join_by(station,Niskin)) %>% 
  #dplyr::select(-Zymo) %>% 
  distinct()

#### QPCR DATA ####
qPCR_unk <- read_csv(here('data','hake_qPCR','Hake eDNA 2019 qPCR results 2021-01-04 results.csv'),
                     col_types = 'ccccccdccdcccccclll') %>% 
  rename(tubeID=sample) %>% 
  left_join(.,qPCR.sample.id,by=join_by(tubeID)) %>% 
  # fix some columns from chr to numeric
  mutate(IPC_Ct=str_replace_all(IPC_Ct,"Undetermined","") %>% as.numeric) %>% 
  mutate(hake_copies_ul=str_replace_all(hake_copies_ul,",",""),
         eulachon_copies_ul=str_replace_all(eulachon_copies_ul,",",""),
         lamprey_copies_ul=str_replace_all(lamprey_copies_ul,",","")) %>% 
  mutate(across(contains("copies_ul"),~as.numeric(.))) 
  # this leaves 9208 rows

# Get rid of zymo filtered samples  FIX USEFUL HERE
qPCR_unk <- qPCR_unk %>% filter(is.na(Zymo)) %>% 
                mutate(useful = ifelse(is.na(useful),"missing",useful)) %>% 
                filter(useful %in% c("missing","YES"))
  # 7412 rows remain.

# Specify an inhibition limit for retaining samples.
INHIBIT.LIMIT <- 0.5

# Get rid of samples with dilution == 1 if a dilution series was run on a sample and those that were inhibited
dat.ntc <- qPCR_unk %>% filter(type=="ntc") %>% 
  mutate(IPC_Ct = as.numeric(as.character(IPC_Ct))) %>%
  group_by(qPCR) %>% 
  dplyr::summarise(mean.ntc = mean(IPC_Ct),sd.ntc=sd(IPC_Ct))

# This gets rid of inhibited samples.
qPCR_unk <- left_join(qPCR_unk,dat.ntc,by=join_by(qPCR)) %>% 
  # mutate(mean.ntc = as.numeric(as.character(mean.ntc))) %>%
  mutate(inhibit.val = IPC_Ct-mean.ntc,
         inhibit.bin=ifelse(inhibit.val < INHIBIT.LIMIT ,0,1)) %>%
  filter(inhibit.bin ==0) %>% 
  select(-inhibition_rate)
# This leaves 7266 rows of data.

### Add wash covariate
### CHECK THE SAMPLES THAT HAD TROUBLE WITH WASHING (drop.sample == "30EtOH" or "30EtOHpaired")
dat.wash <- qPCR_unk %>% filter(drop.sample %in% c("30EtOH","30EtOHpaired")) %>% mutate(status="washed")

# find unique depth-station combinations among these stations.
uni.wash <- dat.wash %>% 
  group_by(station,depth,Niskin,drop.sample) %>% 
  summarise(N=n()) %>%
  mutate(status="washed") %>% 
  ungroup()

# Find the paired, but unwashed samples from the remaining samples.
pairs.wash <- qPCR_unk %>% filter(!drop.sample %in% c("30EtOH","30EtOHpaired"))
pairs.wash <- uni.wash %>% 
  dplyr::select(station,depth) %>% 
  left_join(pairs.wash,by = join_by(station, depth)) %>%
  filter(!is.na(Niskin)) %>%
  mutate(status="unwashed")
#99 of these

dat.wash.all <- bind_rows(dat.wash,pairs.wash) %>% arrange(station,depth)
dat.wash.summary <- dat.wash.all %>% group_by(station,depth,Niskin,status) %>% summarise(N=n()) %>% 
              arrange(station,depth,status) %>% as_tibble() 

# There are 27 paired samples with which to estimate the effect of the 30% EtOH treatment
dat.wash.summary %>% count(status)

# add indicator for membership in 30EtOH club and associated pairs
dat.wash.all <- dat.wash.all %>% 
  # if wash.indicator is 0 = normal sample. 1 = washed with 30% etoh. 2= pair of washed with 30% EtOH sample
  mutate(wash.indicator = ifelse(status == "washed",1,2)) %>%
  # for STAN, just keep track as a binary variable whether each sample was washed or not
  # 1=washed, 0=unwashed
  mutate(wash_idx=as.numeric(wash.indicator==1)) %>% 
  dplyr::select(-status)

# find samples that were washed with 30% EtOH, exclude them from dat.samp,
# then add them back in with needed indicator variables
qPCR_unk <- qPCR_unk %>% 
  mutate(wash.indicator=0,wash_idx=0) %>% 
  filter(!tubeID %in% unique(dat.wash.all$tubeID)) %>% 
  bind_rows(.,dat.wash.all)
# Still at 7266 rows of data.

# Filter out various controls, field negatives, ntc, etc.
qPCR_unk <- qPCR_unk %>%
          filter(type == "unknowns") %>% # this gets rid of 646 rows.
          filter(!is.na(utm.lat)) # gets rid of 9 replicates (3 unique tubes) filtered from other hake projects
# Still at 6614 rows of data.

# Classify each listed depth into one of a few categories.
qPCR_unk <- qPCR_unk %>% mutate(depth_cat=case_when(depth < 25 ~ 0,
                                                    depth ==25 ~ 25,  
                                                    depth > 25  & depth <= 60  ~ 50,
                                                    depth > 60  & depth <= 100 ~ 100,
                                                    depth > 119 & depth <= 150 ~ 150,
                                                    depth > 151 & depth <= 200 ~ 200,
                                                    depth > 240 & depth <= 350 ~ 300,
                                                    depth > 400 & depth <= 500 ~ 500))

# Only keep observations with a dilution of 1 or 0.2.
qPCR_unk <- qPCR_unk %>% filter(dilution %in% c(0.2,1))
# Down to 5394 rows

# Add volume for offset
qPCR_unk <- qPCR_unk %>% mutate(volume_offset = volume / 2.5)

# Make a summary file for each depth-location combination.  Call this station_dat
station_dat <- qPCR_unk %>% 
                dplyr::select(year,tubeID, station,lat,lon,utm.lat,utm.lon,depth,depth_cat) %>% 
                distinct() %>% group_by(year,station, lat, lon, utm.lat,utm.lon,depth,depth_cat) %>%
                count() %>% rename(n_tube_station_depth = n) %>% ungroup() %>% 
                mutate(station_depth_idx = 1:nrow(.))
stat2 <- qPCR_unk %>% 
                dplyr::select(year,tubeID, station,lat,lon,utm.lat,utm.lon) %>% 
                distinct() %>% group_by(year,station, lat, lon, utm.lat,utm.lon) %>% 
                count() %>% rename(n_tube_station = n) %>% ungroup() %>% 
                mutate(station_idx = 1:nrow(.))
station_dat <- left_join(station_dat,stat2)

# Merge back into the qPCR_unk
qPCR_unk <- qPCR_unk %>% left_join(.,station_dat)


### RANDOM EFFECT DESIGN MATRICES
# Make a matrix for random effect associated with each station-depth combination at observation level 
  form <- "year ~ 0 + factor(tubeID): factor(n_tube_station_depth)"
  model_frame   <- model.frame(form, qPCR_unk)  
  A <- model.matrix(as.formula(form), model_frame)
    # get rid of factor levels (columns) that == 0
    X_bio_rep_obs <- A[,which(colSums(A)>0)]
    # Only keep columns that have at least 2 replicates
    THESE <- substr(colnames(X_bio_rep_obs),nchar(colnames(X_bio_rep_obs)),nchar(colnames(X_bio_rep_obs)))
    THESE <- as.numeric(THESE)
    X_bio_rep_obs <- X_bio_rep_obs[,which(THESE>1)]
    N_bio_rep_RE <- ncol(X_bio_rep_obs)

# Make a matrix for random effect associated with each station-depth combination at tubeID level 
  #reduce qPCR to just a single value for each tubeID
    tube_dat <-   qPCR_unk %>% distinct(tubeID,station_depth_idx,station_idx,n_tube_station_depth,depth_cat)
    form <- "depth_cat ~ 0 + factor(tubeID): factor(n_tube_station_depth)"
    model_frame   <- model.frame(form, tube_dat)  
    A <- model.matrix(as.formula(form), model_frame)
    # get rid of factor levels (columns) that == 0
    X_bio_rep_tube <- A[,which(colSums(A)>0)]
    # Only keep columns that have at least 2 replicates
    THESE <- substr(colnames(X_bio_rep_tube),nchar(colnames(X_bio_rep_tube)),nchar(colnames(X_bio_rep_tube)))
    THESE <- as.numeric(THESE)
    X_bio_rep_tube <- X_bio_rep_tube[,which(THESE>1)]

    # Check to make sure the order is correct
    identical(colnames(X_bio_rep_tube),colnames(X_bio_rep_obs))
    
### ONE MORE FIXED EFFECT DESIGN MATRIX AT THE TUBE LEVEL.
    form <- "depth_cat ~ 0 + factor(station_depth_idx)"
    model_frame   <- model.frame(form, tube_dat)  
    X_station_depth_tube <- model.matrix(as.formula(form), model_frame)
    
### Make random effect that sums to zero for the tubes.
    tube_dat <- tube_dat %>% group_by(station_depth_idx) %>% mutate(rep_id = rep(1:n()))
    tube_trim <- tube_dat %>% ungroup() %>% filter(rep_id != n_tube_station_depth) %>% mutate(bio_rep_idx1= 1:nrow(.))    
    N_bio_rep1 <- nrow(tube_trim)
    tube_dat <-  tube_dat %>% ungroup() %>%  mutate(bio_rep_idx2= 1:nrow(.))    
    
    
    
    
    
    
    
    for(i in 1:N_bio_rep2){
      temp = 0
      
      
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
# Make a new metadata file that has all of the requisite stuff.
META <- qPCR_unk %>% dplyr::select(tubeID, station,lat,lon,depth,depth_cat,wash_idx) %>% distinct()
  
####
qPCR_std <- read_csv(here('data','hake_qPCR','Hake eDNA 2019 qPCR results 2020-01-04 standards.csv')) %>% 
  rename(tubeID=sample)
 
#### METABARCODING DATA ####

#species and all assigned metabarcoding data
sp <- read.table(here("data/mocks/species_list.txt"), sep = "\t") %>% 
  rename("Species" = "V1")
q <- read.csv(here("data/mocks/mock_paper_ASV_assignments.csv"), row.names = 1) %>% 
  # filter to species we know we want
  filter(Species %in% sp$Species) %>% 
  filter(!is.na(Species)) %>% 
  select(Sample, Reads, Marker, Species) %>% 
  # separate into mock community types
  mutate(BSA = ifelse(str_detect(Sample, "BSA"), 1, 0),
         TD = ifelse(str_detect(Sample, "TD"), 1, 0),
         SKEW = ifelse(str_detect(Sample, "SKEW"), 1, 0)) %>% 
  mutate(Sample = str_replace_all(Sample, "BSA\\.", ""),
         Sample = str_replace_all(Sample, "TD\\.", ""),
         Sample = str_replace_all(Sample, "SKEW\\.", "")) %>% 
  separate(Sample, into = c("Marker", "Taq", "rep")) %>% 
  # now regroup and sum total reads by taxa
  group_by(Marker, Taq, rep, Species, BSA, TD, SKEW) %>% 
  summarise(Reads = sum(Reads)) %>%  #collapse ASVs into taxa
  ungroup()

# rename Engraulis (anchovy) to be at the "end" of an alphabetical list; this will be our metabarcoding reference species
q$Species[which(q$Species == "Engraulis mordax")] <- "Zz_Engraulis mordax"

#we can't estimate species that are only rarely observed, so we are going to drop those. 
# in this case, there are two species will substantially less total reads that any others (<1000 reads across all samples)
allspecies <- unique(q$Species)
lose <- q %>% group_by(Species) %>% summarise(s = sum(Reads)) %>% arrange(s) %>% head(2) %>% pull(Species)
keep <- allspecies[-which(allspecies %in% lose)]
# keep <- c(keep, "Zz_Engraulis mordax") #to make Engraulis the reference species

q <- q %>% 
  filter(Species %in% keep)

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
  filter(Project == "52193") %>%  #just keep hake-cruise samples
  #mutate(sample = as.numeric(sample)) %>% 
  # separate(Sample_name, into = c("Primer", "Project", "sample", "Dilution", "Rep", "Well"), sep = c("-|_")) %>% 
  left_join(.,META %>% rename(sample = tubeID))

## PRE-PROCESSING DATA FOR STAN ##

keep <- intersect(unique(mfu$BestTaxon), template_props$Species) %>% intersect(fish)
keep <- c(keep, "Zz_Engraulis mordax")
keep <- keep[-which(keep == "Oncorhynchus tshawytscha")] #remove for the moment, pending better tissue sample
mfu$BestTaxon[mfu$BestTaxon == "Engraulis mordax"] <- "Zz_Engraulis mordax"
mfu$BestTaxon[mfu$BestTaxon == "Merluccius"] <- "Merluccius productus"
mfu$BestTaxon[mfu$BestTaxon == "Sardinops"] <- "Sardinops sagax"
mfu$BestTaxon[mfu$BestTaxon == "Clupea"] <- "Clupea pallasii"
mfu$BestTaxon[mfu$BestTaxon == "Clupeidae"] <- "Clupea pallasii"
mfu$BestTaxon[mfu$BestTaxon == "Trachurus"] <- "Trachurus symmetricus"

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

#samples w no hake MB reads  
no_hake_MB <- colnames(mfu)[2:ncol(mfu)][which(mfu[which(mfu$Species == "Merluccius productus"), 2:ncol(mfu)] == 0)] %>% as.numeric()

mfu <- mfu %>% 
  pivot_longer(cols = -Species, names_to = "Sample", values_to = "Nreads") %>% 
  group_by(Sample) %>%
  mutate(totalReads = sum(Nreads)) %>% 
  filter(totalReads > 1000) %>% #impose total read-depth min
  ungroup() %>% 
  select(-totalReads) %>% 
  filter(!Sample %in% no_hake_MB) %>%  #filter out samples w no hake MB reads
  mutate(biol = 1, tech = 1) %>% 
  rename(species = Species) %>% 
  mutate(station = match(Sample, unique(Sample))) #let station = unique sampling location, treating different depths in the same x-y plane as independent

