# Load data for QM-qPCR joint model

# Packages ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
options(tidyverse.quiet = TRUE, dplyr.summarise.inform = FALSE,dplyr.left_join.inform = FALSE)
suppressWarnings(suppressPackageStartupMessages(library(rstan)))
options(mc.cores = parallel::detectCores())
rstan_options(threads_per_chain = 4)
rstan_options(auto_write = TRUE)
suppressWarnings(suppressPackageStartupMessages(library(compositions)))
suppressWarnings(suppressPackageStartupMessages(library(MCMCpack)))

# Functions ------------------------------------------------------------------------------
logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
} #from https://gist.github.com/aufrank/83572

select <- dplyr::select

## Raw Data ------------------------------------------------------------------------------

#hake sample attributes
atts <- read_csv(here("data","metadata","Hake_2019_metadata.csv"),col_types = cols()) %>% 
  # keep only the relevant and unique columns from this data 
  select(station,Niskin,year,month,day,transect,lat,lon,utm.lat,utm.lon,bottom.depth.consensus,transect.dist.km) %>%
  distinct()

qPCR.sample.id <- suppressMessages(read_csv(here('data','hake_qPCR','Hake eDNA 2019 qPCR results 2023-02-10 sample details.csv'),
                           col_select = all_of(c("Tube #", "CTD cast","Niskin","depth","drop.sample","field.negative.type","water.filtered.L")),
                           col_types=cols()))

### SAMPLE IDs ###
qPCR.sample.id <- qPCR.sample.id %>%  
  rename(tubeID=`Tube #`,
         station=`CTD cast`,
         volume=water.filtered.L) %>%
  distinct() %>% 
  # empty stations or extraneous rows
  filter(!(station=="N/A"|station=="-")) %>%
  # organize and clean sample depth attribute
  mutate(depth=ifelse(Niskin=="sfc","0",depth)) %>%
  mutate(depth=ifelse(depth=="sfc","0",depth)) %>%
  mutate(depth=ifelse(depth=="300/150","300",depth)) %>%
  mutate(depth=as.numeric(depth)) %>% 
  left_join(atts,by = join_by(station,Niskin)) %>%
  distinct()

# Load qPCR data ------------------------------------------------------------------------------
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

# Get rid of zymo filtered samples
qPCR_unk <- qPCR_unk %>% 
  filter(is.na(Zymo)) %>% 
  mutate(useful = ifelse(is.na(useful),"missing",useful)) %>% 
  filter(useful %in% c("missing","YES"))
# 7412 rows remain.


## Add qPCR covariates -------------------------------------------------------------------------
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


### Add wash covariate --------------------------------------------------------------------------
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
#96 of these

dat.wash.all <- bind_rows(dat.wash,pairs.wash) %>% arrange(station,depth)
dat.wash.summary <- dat.wash.all %>% group_by(station,depth,Niskin,status) %>% summarise(N=n()) %>% 
  arrange(station,depth,status) %>% as_tibble() 

# There are 27 paired samples with which to estimate the effect of the 30% EtOH treatment
# to check: dat.wash.summary %>% count(status)

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

### Filter out various controls -----------------------------------------------------------------
# Filter out various controls, field negatives, ntc, etc.
qPCR_unk <- qPCR_unk %>%
  filter(type == "unknowns") %>% # this gets rid of 646 rows.
  filter(!is.na(utm.lat)) # gets rid of 9 replicates (3 unique tubes) filtered from other hake projects
# Now 6614 rows of data.


### Classify each listed depth into one of a few categories-----------------------------------------------------------------
qPCR_unk <- qPCR_unk %>% mutate(depth_cat=case_when(depth < 25 ~ 0,
                                                    depth ==25 ~ 25,  
                                                    depth > 25  & depth <= 60  ~ 50,
                                                    depth > 60  & depth <= 100 ~ 100,
                                                    depth > 119 & depth <= 150 ~ 150,
                                                    depth > 151 & depth <= 200 ~ 200,
                                                    depth > 240 & depth <= 350 ~ 300,
                                                    depth > 400 & depth <= 500 ~ 500))

### Filter sample dilution ----------------------------------------------------------------------
# Only keep observations with a dilution of 1 or 0.2.
qPCR_unk <- qPCR_unk %>% filter(dilution %in% c(0.2,1))
# Down to 5394 rows

### Add volume offset ---------------------------------------------------------------------------
qPCR_unk <- qPCR_unk %>% mutate(volume_offset = volume / 2.5)

## Summary qPCR files ---------------------------------------------------------------------------
# Make a summary file for each depth-location combination.  Call this station_dat
station_dat_depth <- qPCR_unk %>% 
  dplyr::select(year,tubeID, station,lat,lon,utm.lat,utm.lon,depth,depth_cat) %>% 
  distinct() %>% group_by(year,station, lat, lon, utm.lat,utm.lon,depth,depth_cat) %>%
  count() %>% rename(n_tube_station_depth = n) %>% ungroup() %>% 
  mutate(station_depth_idx = 1:nrow(.))
# Same, except for location only (i.e., depth-integrated)
station_dat_flat <- qPCR_unk %>% 
  dplyr::select(year,tubeID, station,lat,lon,utm.lat,utm.lon) %>% 
  distinct() %>% 
  group_by(year,station, lat, lon, utm.lat,utm.lon) %>% 
  count() %>% 
  rename(n_tube_station = n) %>% ungroup() %>% 
  mutate(station_idx = 1:nrow(.))

station_dat <- left_join(station_dat_depth,station_dat_flat,by = join_by(year, station, lat, lon, utm.lat, utm.lon))

## Merging and indexing ------------------------------------------------------------------------
# Merge back into the qPCR_unk
qPCR_unk <- qPCR_unk %>% left_join(.,station_dat,by = join_by(station, depth, year, lat, lon, utm.lat, utm.lon, depth_cat))

# Make tube_idx
tube_dat <-   qPCR_unk %>% distinct(tubeID,station_depth_idx,station_idx,n_tube_station_depth,depth_cat)
tube_dat <- tube_dat %>% mutate(tube_idx = row_number())

# Merge back into the qPCR_unk
qPCR_unk <- qPCR_unk %>% left_join(tube_dat,by = join_by(tubeID, depth_cat, n_tube_station_depth, station_depth_idx, station_idx))

# Make a metadata file that has all of the requisite stuff post-filtering.
META <- qPCR_unk %>% dplyr::select(tubeID, station,lat,lon,depth,depth_cat,wash_idx) %>% distinct()
write_rds(META,here('data','metadata','Hake_qPCR_META_after_data_prep.rds'))

#### Standards don't need much formatting
qPCR_std <- read_csv(here('data','hake_qPCR','Hake eDNA 2019 qPCR results 2020-01-04 standards.csv'),col_types=cols()) %>% 
  rename(tubeID=sample)

## Load metabarcoding data ---------------------------------------------------------------------

# Annotation database
db <- read.csv(here('data','MFU_database.csv'), row.names = 1)

# two tiny consolidations in the database for species we care about
db <- db %>% 
  mutate(BestTaxon=if_else(Tax_list=='Microstomus pacificus, Myzopsetta proboscidea',
                           'Microstomus pacificus',BestTaxon)) %>% 
  mutate(BestTaxon=if_else(Tax_list=='Hozukius emblemarius, Sebastes entomelas, Sebastes matsubarae',
                           'Sebastes entomelas',BestTaxon))


## Load sequencing data batches------------------------------------------------------------------

# function to load a run (because this will be the same for each different run, we can apply a consistent function to load them)
load_asv_table <- function(seq_run_number){
  # filepath
  fp <- here("data","fish",paste0("MURI",seq_run_number,"_MFU_ASV_table.csv"))
  # load
  read_csv(fp,col_types=cols())%>% 
    left_join(db,by = join_by(Hash)) %>%  
    filter(!is.na(BestTaxon)) %>% 
    mutate(Sample_name=str_replace(Sample_name,"-1-","-")) %>%
    mutate(Sample_name=str_replace(Sample_name,"-1_","_")) %>% 
    separate(Sample_name, into = c("Primer", "Project", "sample", "Dilution", "Well"), sep = c("-|_"), remove = F) %>% 
    mutate(Rep = "1") %>% 
    relocate(c("Primer", "Project", "sample", "Dilution", "Rep", "Well"))
}

# apply the above function
mfu <- purrr::map(c(304,313:318),load_asv_table)%>% 
  
  # row-bind all batches of samples
  list_rbind() %>% 
  
  # rename sample to tubeID
  rename(tubeID=sample) %>% 
  
  # summarise reads by taxon
  group_by(Primer, Project, tubeID, Dilution, Rep, BestTaxon) %>% 
  summarise(nReads = sum(nReads)) %>%
  ungroup() %>% 
  
  # fill in explicit zeroes for missing taxon/sample combinations
  complete(nesting(Primer,Project,tubeID,Dilution,Rep),BestTaxon,fill=list(nReads=0)) %>% 
  
  ungroup() %>% 
  filter(!grepl("positive",Project)) %>%
  filter(!grepl("NTC",Project)) %>%
  filter(Project == "52193") %>%  #just keep hake-cruise samples
  #mutate(sample = as.numeric(sample)) %>% 
  # separate(Sample_name, into = c("Primer", "Project", "sample", "Dilution", "Rep", "Well"), sep = c("-|_")) %>% 
  left_join(META,by=join_by(tubeID))

## Seq_depth filter factor ---------------------------------------------------------------------
# Which samples have sequencing depth >1000 reads?
samp_over_1000 <-mfu %>% 
  group_by(tubeID) %>% 
  summarise(sum_R=sum(nReads)) %>% 
  filter(sum_R>1000) %>% 
  pull(tubeID)

# Load mock data ------------------------------------------------------------------------------
#Collect all the mock rds files
mock_list <- list.files(here('data/mocks'),recursive = T,pattern = '.rds')

# as with the MFU dataset, define a function to consistently load mock data
load_mock_file <- function(fp){
  df <- read_rds(here('data','mocks',fp)) %>%
    
    # remove species index because we don't need it
    select(-sp_idx) %>% 
    
    # pivot just the reads from wide to long form
    pivot_longer(cols = where(is.integer), names_to = 'Sample',values_to = 'Reads') %>% 
    
    # rename the columns
    set_names(c('Species','b_proportion','Sample','Reads'))
}

# load and combine mocks
mock_comb <- purrr::map(mock_list,load_mock_file) %>% 
  list_rbind() %>% 
  arrange(Species)

# Use Sample string to extract attributes we want (this shouldn't remove any data)
mock <- mock_comb %>%
  mutate(Sample = gsub('^M', 'Mock', Sample) %>%
           gsub('\\.e\\.', '_even_', .) %>%
           gsub('\\.s\\.', '_skew_', .) %>%
           gsub('_S[0-9]+$', '', .)) %>% 
  separate(Sample, into = c("Mock_samp", "mock_type", "Rep"), sep = c("_"), remove = T) %>% 
  relocate(c('Mock_samp','mock_type','Rep','Species','Reads','b_proportion')) %>% 
  setNames(c('Sample','Mock_type','Rep','species','Nreads','b_proportion')) %>% 
  mutate(Primer='MFU')

# look at them?
# x <- mock_comb %>% group_by(Sample) %>% mutate(prop=Reads/sum(Reads)) %>% ungroup()
y <- mock %>% unite(sn,Sample,Mock_type,Rep) %>% ungroup()
y %>% ggplot(aes(sn,b_proportion,fill=species))+geom_col(position='stack')+theme(axis.text.x=element_text(angle=45))
x %>% ggplot(aes(Sample,prop,fill=Species))+geom_col(position='stack')+theme(axis.text.x=element_text(angle=45))

########################################################
mock <- mock %>% rename(Sample_short=Sample) %>% mutate(Sample = paste0(Sample_short,"_",Mock_type))
########################################################

# Species in the mocks
species_mock_list <- mock %>% 
  # Add "Zz" to make hake the QM reference species
  mutate(species=ifelse(species=="Merluccius productus","Zz_Merluccius productus",species)) %>% 
  # pull out all unique speices
  pull(species) %>% 
  unique() %>% 
  sort()
  

## List of species to keep ---------------------------------------------------------------------------

keep <- c('Clupea pallasii',
          #'Diogenichthys atlanticus', #we don't have lantern fish anywhere else other than in run 305 which I don't know why is turned off
          'Engraulis mordax',
          'Hippoglossus stenolepis',
          'Leuroglossus stilbius',
          'Zz_Merluccius productus', # QM reference species that we recoded to be last in the list
          'Microstomus pacificus',
          'Sardinops sagax',
          'Scomber japonicus',
          'Sebastes entomelas',
          'Stenobrachius leucopsarus',
          'Tactostoma macropus',
          'Tarletonbeania crenularis',
          'Thaleichthys pacificus',
          'Trachurus symmetricus')%>% sort()

#Species in mock discarded
cat('\n');cat('Species from mock left out::');cat('\n')
writeLines(unique(species_mock_list)[!species_mock_list%in%keep]);cat('\n')

# Filter the mocks to just the species we want
mock <- mock %>% 
  filter(species%in%keep) %>%
  # renormalize b_proportions after omitting species from Mocks

  group_by(Sample,Rep) %>% 
  mutate(b_proportion = b_proportion/sum(b_proportion)) %>% 
  ungroup()

## Some manual curation of MB data ---------------------------------------------------------------------

mfu <- mfu %>% 
  mutate(BestTaxon=case_when(
    BestTaxon == "Clupea"~ "Clupea pallasii",
    BestTaxon == "Clupeidae"~ "Clupea pallasii",
    # BestTaxon == "Diogenichthys"~ "Diogenichthys atlanticus" #we don't have lantern fish anywhere else other than in run 305 which I don't know why is turned off
    BestTaxon == "Hippoglossus"~ "Hippoglossus stenolepis",
    BestTaxon == "Leuroglossus"~ "Leuroglossus stilbius",
    BestTaxon == "Merluccius"~ "Zz_Merluccius productus",
    BestTaxon == "Sardinops"~ "Sardinops sagax",
    BestTaxon == "Scomber"~ "Scomber japonicus",
    BestTaxon == "Stenobrachius"~ "Stenobrachius leucopsarus",
    BestTaxon == "Tarletonbeania"~ "Tarletonbeania crenularis",
    BestTaxon == "Trachurus"~ "Trachurus symmetricus",
    BestTaxon == "Merluccius"~ "Zz_Merluccius productus",
    BestTaxon == "Merluccius productus"~ "Zz_Merluccius productus",
    TRUE ~ BestTaxon
  ))

## Filter for MB data ---------------------------------------------------------------------

# Filter to selected study species
mfu <- mfu %>%
  rename(Species = BestTaxon) %>% 
  filter(Species %in% keep)

cat('Metabarcoding summary: ');cat('\n')
mfu %>% group_by(Species) %>% 
  summarise(Sum_reads=sum(nReads)) %>% 
  mutate(prop=(Sum_reads/sum(Sum_reads))*100) %>% print()

# METADATA for metabarcoding
mfu_META <- mfu %>% 
  drop_na(tubeID, lat, Species, nReads) %>%
  dplyr::select(tubeID,station:wash_idx) %>%
  mutate(sampleidx = match(tubeID, unique(tubeID))) %>% 
  distinct()


# fill in zeroes for missing Species/tube combinations
mfu <- mfu %>%
  select(tubeID, Species, nReads) %>%
  drop_na() %>%
  group_by(Species, tubeID) %>%
  summarise(nReads = sum(nReads)) %>%
  ungroup() %>% 
  complete(Species,tubeID,fill=list(nReads=0))

cat('\n');cat('# of samples in MB samples: ')
cat(length(unique(mfu$tubeID)))

# Filter to samples with >1000 read depth
mfu <- mfu %>% 
  filter(tubeID%in%samp_over_1000) %>% 
  ungroup() %>%
  mutate(biol = 1, tech = 1) %>% 
  rename(species = Species)

cat('\n');cat('# of samples in MB samples with >1000 sequencing depth: ')
cat(length(unique(mfu$tubeID)))

# find samples w no hake MB reads  
no_hake_MB <- mfu %>% 
  filter(nReads==0,species=="Zz_Merluccius productus") %>% 
  pull(tubeID)

# find samples that have hake present
cat('\n');cat('# of samples that have hake present: ')
cat(length(unique(mfu$tubeID[!mfu$tubeID%in%no_hake_MB])))

# find samples with nothing BUT hake
no_others_MB <- mfu %>% 
  filter(species!="Zz_Merluccius productus") %>% 
  group_by(tubeID) %>% summarise(totreads=sum(nReads)) %>% 
  filter(totreads==0) %>% 
  pull(tubeID)

# find samples that have hake and another species present
cat('\n');cat('# of samples that have hake and another species present: ')
mfu %>% 
  filter(!tubeID%in%no_hake_MB) %>% 
  filter(!tubeID%in%no_others_MB) %>% 
  pull(tubeID) %>% unique() %>% length() %>% 
  cat()

# final filter- remove samples with no hake; we leave samples with only hake for now
mfu <- mfu %>% 
  filter(!(tubeID%in% no_hake_MB))


cat('\n');cat('Final # of metabarcoding samples loaded: ')
cat(length(unique(mfu$tubeID)))