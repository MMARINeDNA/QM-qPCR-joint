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
suppressWarnings(suppressPackageStartupMessages(library(MoMAColors)))

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

cat('\n');cat('# of qPCR reactions are: ')
cat(nrow(qPCR_unk))

cat('\n');cat('# of unique qPCR samples are: ')
cat(length(unique(qPCR_unk$tubeID)))

#### Standards don't need much formatting
qPCR_std <- read_csv(here('data','hake_qPCR','Hake eDNA 2019 qPCR results 2020-01-04 standards.csv'),col_types=cols()) %>% 
  rename(tubeID=sample)

## Load metabarcoding data ---------------------------------------------------------------------

# Annotation database
db <- read.csv(here('data','metabarcoding_db','MFU_database.csv'),row.names = 1)

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
  fp <- here('data','metabarcoding',paste0("MURI_",seq_run_number,"_MFU_ASV_table.csv"))
  # load
  read_csv(fp,col_types=cols())%>% 
    left_join(db,by = join_by(Hash)) %>%  
    filter(!is.na(BestTaxon)) %>% 
    mutate(Sample_name=str_replace(Sample_name,"-1-","-")) %>%
    mutate(Sample_name=str_replace(Sample_name,"-1_","_")) %>% 
    separate(Sample_name, into = c("Primer", "Project", "sample", "Dilution", "Well"), sep = c("-|_"), remove = F) %>% 
    mutate(Rep = "1") %>% 
    # mutate(Run = seq_run_number) %>% 
    relocate(c("Primer", "Project", "sample", "Dilution", "Rep", "Well"))
    # relocate(c("Primer",'Run', "Project", "sample", "Dilution", "Rep", "Well"))
}

# apply the above function
mfu <- purrr::map(c(304,313:318),load_asv_table)%>% 
  
  # row-bind all batches of samples
  list_rbind() %>% 
  
  # rename sample to tubeID
  rename(tubeID=sample) %>% 
  
  # summarise reads by taxon
  group_by(Primer, Project, tubeID, Dilution, Rep, BestTaxon) %>% 
  # group_by(Primer,Run, Project, tubeID, Dilution, Rep, BestTaxon) %>% 
  summarise(nReads = sum(nReads)) %>%
  ungroup() %>% 
  
  # fill in explicit zeroes for missing taxon/sample combinations
  complete(nesting(Primer,Project,tubeID,Dilution,Rep),BestTaxon,fill=list(nReads=0)) %>% 
  # complete(nesting(Primer,Run,Project,tubeID,Dilution,Rep),BestTaxon,fill=list(nReads=0)) %>% 
  
  ungroup() %>% 
  
  filter(Project == "52193" | grepl("positive",Project) | grepl("NTC",Project)) %>%  #just keep hake-cruise samples
  left_join(META,by=join_by(tubeID))

write.csv(mfu,here('data','mfu_raw0.csv'),row.names = F)

mfu <- mfu %>% 
  filter(!grepl("positive",Project))

## Seq_depth filter factor ---------------------------------------------------------------------
# Which samples have sequencing depth >1000 reads?
samp_over_1000 <-mfu %>% 
  group_by(tubeID) %>% 
  summarise(sum_R=sum(nReads)) %>% 
  filter(sum_R>1000) %>% 
  pull(tubeID)

n_samp_over_1000 <- mfu %>% 
  group_by(tubeID) %>% 
  summarise(sum_R=sum(nReads)) %>% 
  filter(sum_R>1000) %>% 
  nrow()

# Load mock data ------------------------------------------------------------------------------
#Collect all the Mock metabarcoding runs csv files
mock_list <- list.files(here('data','metabarcoding_mocks','asv_tables'),recursive = T,pattern = '.csv')

# Curate Mock samples and sample names
mock <- read.csv(here('data','metabarcoding_mocks','asv_tables',mock_list[1])) %>% 
  mutate(Sample_name=gsub('MFU.','',Sample_name)) %>%
  mutate(Sample_name=gsub('\\-','\\_',Sample_name)) %>%
  mutate(Sample_name=gsub('mock','M',Sample_name)) %>%
  mutate(Sample_name=gsub('skewed','s',Sample_name)) %>%
  mutate(Sample_name=gsub('skew','s',Sample_name)) %>%
  mutate(Sample_name=gsub('even','e',Sample_name)) %>%
  mutate(Sample_name=gsub('d10.','',Sample_name)) %>%
  mutate(Sample_name=gsub('d1.','',Sample_name)) %>%
  mutate(Sample_name=gsub('_S\\d+','',Sample_name)) %>%
  mutate(Sample_name=gsub('M_','M_1_',Sample_name)) %>%
  filter(grepl('^M\\_([1234])\\_([es])',Sample_name)) %>%
  rbind(.,
        read.csv(here('data','metabarcoding_mocks','asv_tables',mock_list[2])) %>% 
          mutate(Sample_name=gsub('MFU-CPS-','',Sample_name)) %>%
          mutate(Sample_name=gsub('\\-','\\_',Sample_name)) %>%
          mutate(Sample_name=gsub('mock','M_',Sample_name)) %>%
          mutate(Sample_name=gsub('skew','s',Sample_name)) %>%
          mutate(Sample_name=gsub('even','e',Sample_name)) %>%
          mutate(Sample_name=gsub('d10.','',Sample_name)) %>%
          mutate(Sample_name=gsub('d1.','',Sample_name)) %>%
          mutate(Sample_name=gsub('_S\\d+','',Sample_name)) %>%
          filter(grepl('^M\\_([1234])\\_([es])',Sample_name))
  )

# Annotate Mock ASVs
mock <- mock %>%
  left_join(.,db %>% select(Hash,BestTaxon),by='Hash') %>%
  mutate(BestTaxon=if_else(Hash=='b05eac0ae6bebbb0f133eb32789a54e3cb90ddbe','Cymatogaster aggregata',BestTaxon)) %>%
  mutate(BestTaxon=if_else(Hash=='72d93a313d0244edae1cc9229ebd2525ddaec2ed','Merluccius productus',BestTaxon)) %>%
  mutate(BestTaxon=if_else(Hash=='fe8e01e1080eac0ddd3ab3b57d749d9de0ee1fab','Leuroglossus stilbius',BestTaxon))

# List of species to keep 
keep <- c('Clupea pallasii',
          'Engraulis mordax',
          'Leuroglossus stilbius',
          'Merluccius productus',
          'Microstomus pacificus',
          'Sardinops sagax',
          'Scomber japonicus',
          'Sebastes entomelas',
          'Stenobrachius leucopsarus',
          'Tactostoma macropus',
          'Tarletonbeania crenularis',
          'Thaleichthys pacificus',
          'Trachurus symmetricus')%>% sort()

# Keep only the species of interest and collapse the reads of multiple hash of the same species
mock <- mock %>% filter(BestTaxon%in%keep) %>% 
  group_by(Sample_name,BestTaxon) %>% 
  summarise(nReads=sum(nReads)) %>% 
  arrange(Sample_name) 

# Pivot wider the replicates
mock <- mock %>% 
  mutate(Rep=substr(Sample_name,nchar(Sample_name),nchar(Sample_name))) %>% 
  mutate(Sample_name=substr(Sample_name,0,nchar(Sample_name)-2)) %>% 
  pivot_wider(
    names_from = Rep,      
    values_from = nReads,
    names_prefix = "Rep_") %>% 
  mutate(Sample_name=gsub('M_','Mock',Sample_name)) %>% 
  mutate(Sample_name=gsub('_e','_even',Sample_name)) %>% 
  mutate(Sample_name=gsub('_s','_skew',Sample_name)) %>% 
  as_tibble()

# Add initial proportion data
mock_ini_prop_1 <- 
  read.csv(here('data','metabarcoding_mocks','initial_proportions','M2_M3_M4_expected_props_corrected_oct_24.csv'))

mock_ini_prop_2 <- 
  read.csv(here('data','metabarcoding_mocks','initial_proportions','Mock_true_proportions_rpk.csv')) %>% 
  rename(BestTaxon='Species',mtDNA='template_prop',gDNA='genomic_prop') %>% 
  mutate(mock='mock1',even_skew='even') %>% 
  select(BestTaxon,gDNA,mtDNA,mock,even_skew) %>% 
  pivot_longer(cols = c(gDNA, mtDNA), 
               names_to = "prop_type",
               values_to = "Prop") %>% 
  arrange(mock,even_skew,prop_type)

mock_ini_prop_3 <-
  read.csv(here('data','metabarcoding_mocks','initial_proportions','mtDNA_props_skewed.csv')) %>% 
  filter(Class=='Actinopteri') %>% 
  rename(BestTaxon='Species') %>% 
  mutate(Prop=mtDNA_prop_skewed/sum(mtDNA_prop_skewed)) %>%
  mutate(mock='mock1',even_skew='skew',prop_type='mtDNA') %>% 
  arrange(mock,even_skew,prop_type) %>% 
  select(BestTaxon,mock,even_skew,prop_type,Prop)

# Combine the initial proportion data
mock_ini_prop <- 
  mock_ini_prop_1 %>% 
  rbind(.,mock_ini_prop_2) %>% 
  rbind(.,mock_ini_prop_3) %>% 
  arrange(mock,even_skew,prop_type)

# Filter and curate the initial proportion data
mock_ini_prop <- 
  mock_ini_prop %>% filter(prop_type=='mtDNA') %>% 
  filter(BestTaxon %in% keep) %>% 
  rename(species='BestTaxon',
         Mock_type = 'even_skew',
         b_proportion = 'Prop',
         Sample_short='mock') %>% 
  mutate(Primer='MFU',
         Sample = paste0(Sample_short,'_',Mock_type)) %>% 
  select(Sample_short,Mock_type,species,b_proportion,Primer,Sample) %>% 
  mutate(Sample_short=gsub('mock','Mock',Sample_short)) %>% 
  mutate(Sample=gsub('mock','Mock',Sample)) %>% 
  arrange(species,Sample) %>% 
  mutate(x=paste0(Sample,'_',species)) %>%
  as_tibble()

# Join the Mock MB reads with initial proportions
mock <-
  mock %>%
  mutate(x=paste0(Sample_name,'_',BestTaxon)) %>% 
  left_join(.,mock_ini_prop,by='x') %>% 
  filter(!is.na(b_proportion)) %>% 
  select(-BestTaxon,-x,-Sample_name) %>% 
  pivot_longer(
    cols = c(Rep_1,Rep_2,Rep_3), 
    names_to = "Rep",
    values_to = "Nreads") %>% 
  mutate(Rep=gsub('Rep_','',Rep)) %>% 
  arrange(species,Sample_short) %>% 
  select(Sample_short,Mock_type,Rep,species,Nreads,b_proportion,Primer,Sample) %>% 
  group_by(Sample,Rep) %>% 
  mutate(b_proportion = b_proportion/sum(b_proportion)) %>% 
  ungroup()

cat('\n');cat('Mock initial proportions are: ');cat('\n')
cat(mock %>% select(-Sample_short,-Mock_type,-Primer) %>% 
      group_by(species,Sample) %>% 
      summarise(b_proportion=mean(b_proportion)) %>% 
      mutate(b_proportion=round(b_proportion,3)) %>% 
      pivot_wider(names_from = Sample,
                  values_from = b_proportion) %>% 
      select(species,Mock1_even,Mock1_skew,Mock2_even,Mock2_skew, Mock3_even, Mock3_skew, Mock4_even,Mock4_skew) %>% 
      format(., justify = "left"), 
    sep = "\n")

y <- mock %>% unite(sn,Sample_short,Mock_type,Rep) %>% ungroup()
y %>% 
  ggplot(aes(sn,Nreads,fill=species))+
  geom_col(position='stack')+
  scale_fill_manual(values=moma.colors('Lupi',13))+
  ylim(0,220000)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45))

## Some manual curation of MB data ---------------------------------------------------------------------

mfu <- mfu %>% 
  mutate(BestTaxon=case_when(
    BestTaxon == "Clupea"~ "Clupea pallasii",
    BestTaxon == "Clupeidae"~ "Clupea pallasii",
    BestTaxon == "Hippoglossus"~ "Hippoglossus stenolepis",
    BestTaxon == "Leuroglossus"~ "Leuroglossus stilbius",
    BestTaxon == "Sardinops"~ "Sardinops sagax",
    BestTaxon == "Scomber"~ "Scomber japonicus",
    BestTaxon == "Stenobrachius"~ "Stenobrachius leucopsarus",
    BestTaxon == "Tarletonbeania"~ "Tarletonbeania crenularis",
    BestTaxon == "Trachurus"~ "Trachurus symmetricus",
    BestTaxon == "Merluccius"~ "Merluccius productus",
    BestTaxon == "Merluccius productus"~ "Merluccius productus",
    TRUE ~ BestTaxon
  ))

write.csv(mfu,here('data','mfu_raw1.csv'),row.names = F)

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
  filter(nReads==0,species=="Merluccius productus") %>% 
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

mfu <-
  mfu %>%
  mutate(species=if_else(species== "Merluccius productus","Zz_Merluccius productus",species)) %>% 
  arrange(species)
mock <- mock %>%
  mutate(species=if_else(species== "Merluccius productus","Zz_Merluccius productus",species)) %>% 
  arrange(species)


read.csv(here('data','mfu_raw1.csv')) %>% 
  filter(BestTaxon=='Merluccius productus') %>% 
  group_by(tubeID) %>% 
  summarise(nReads=sum(nReads)) %>% 
  ungroup() %>% 
  mutate(pres=if_else(nReads>0,1,0)) %>% 
  count(pres) %>% rename(number_of_samp='n') %>% 
  mutate(prop_presence=number_of_samp/sum(number_of_samp))

qPCR_unk %>% 
  mutate(pres=if_else(hake_Ct=='Undetermined',0,1)) %>% 
  count(pres) %>% rename(number_of_samp='n') %>% 
  mutate(prop_presence=number_of_samp/sum(number_of_samp))
