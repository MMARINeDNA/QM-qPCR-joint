# Plot raw MB reads data, for data included in the joint model
library(tidyverse)

mb_unk_reads_dat <- stan_data_joint %>%
  pluck('sample_data') %>% 
  mutate(rep=row_number()) %>% 
  pivot_longer(-rep,names_to = 'taxa',values_to='reads') %>% 
  # add total reads
  group_by(rep) %>% 
  mutate(sumreads=sum(reads),
         prop=reads/sumreads) %>% 
  ungroup()
mb_link_1 <- mfu %>% 
  mutate(mb_link=match(Sample,unique(Sample))) %>% 
  distinct(Sample,.keep_all = T)

mb_link_2 <- tube_dat %>% 
  distinct(tubeID,tube_idx,station_depth_idx,station_idx) %>%
  #distinct(plate_idx,qpcr_sample_idx,.keep_all = T) %>%
  left_join(mb_link_1,by=join_by(tubeID==Sample)) %>% 
  filter(!is.na(mb_link)) %>%
  arrange(tube_idx) 
# with spatial info
mb_join_key <- mb_link_2 %>% left_join(mfu_META,by=join_by(tubeID==sample)) %>% 
  dplyr::select(tubeID,tube_idx,lat,lon,depth,depth_cat,wash_idx) %>% 
  mutate(rep=row_number())

# join
mb_rawdat <- mb_unk_reads_dat %>% 
  left_join(mb_join_key,by=join_by(rep)) %>% 
  dplyr::select(-rep)
glimpse(mb_rawdat)

mb_rawdat %>%
  ggplot(aes(lon,lat))+
  geom_point()+
  facet_wrap(~factor(depth_cat),nrow=1)+
  coord_equal()

mb_rawdat %>%
  filter(taxa=="sp_6") %>% 
  ggplot(aes(lon,lat,size=reads))+
  geom_point()+
  facet_wrap(~factor(depth_cat),nrow=1)+
  coord_equal()+
  labs(title="Anchovy (raw reads)")


mb_rawdat %>%
  rename(Station=tubeID) %>% 
  filter(taxa=="sp_6") %>% 
  ggplot(aes(Station,reads))+
  geom_point()+
  facet_wrap(~factor(depth_cat),nrow=1)+
  labs(title="Anchovy (raw reads)")
  
test <- read_csv(here('data','metadata','metadata_Gled_10.2.24.csv'))
test2 <- mfu_META %>% mutate(sample_num=as.numeric(sample))
test3 <- test %>% rename(depth_Gled=depth) %>% left_join(test2,by=c("sample"="sample_num"))
glimpse(test)
glimpse(test3)
