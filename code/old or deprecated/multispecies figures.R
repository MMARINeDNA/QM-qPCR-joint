# FIGURE 2: MULTISPECIES SUMMARY
library(tidyverse)
library(corrplot)
library(sf)
library(here)
library(rnaturalearth)
# 
# op_win_cor <- op_dat %>% 
#   mutate(stock = factor(stock, level=op_riv_order)) %>%
#   filter(stock != "Mosquito Creek" & stock != "Moclips") %>%
#   dplyr::select(year, stock, total) %>%
#   na.omit() %>%
#   ungroup() %>% 
#   pivot_wider(.,id_cols=c("year"), names_from=c("stock"), values_from = c("total"))
# tmp <- op_riv_order[op_riv_order != "Mosquito Creek" & op_riv_order != "Moclips"]
# 
# corrplot(op_win_cor[,c("year",tmp)] %>% na.omit() %>% dplyr::select(-year) %>% cor(.), tl.col=
#            "black", tl.cex = 0.5)
# 
# op_win_cor <- op_dat %>% 
#   filter(stock != "Mosquito Creek" & stock != "Moclips") %>%
#   dplyr::select(year, stock, total) %>%
#   na.omit() %>%
#   ungroup() %>% 
#   pivot_wider(.,id_cols=c("year"), names_from=c("stock"), values_from = c("total"))
# 
# corrplot(op_win_cor %>% na.omit() %>% dplyr::select(-year) %>% cor(.), tl.col=
#            "black", addrect = 4, order="hclust", hclust.method="complete",
#          tl.cex = 0.5)


## LOAD DATA
# joint model output
fit <- read_rds(here('data','log_D_fitted.rds'))
# projected (gridded)
preds <- read_rds(here('data','all_smooth_predictions_13sp.rds'))

predgrid <- preds[[1]] %>% distinct(grID,depth_cat,x) %>% st_as_sf()

coaststates <- ne_states(country='United States of America',returnclass = 'sf') %>% 
  filter(name %in% c('California','Oregon','Washington','Nevada')) %>% 
  st_transform("+proj=utm +zone=10 +datum=WGS84 +units=km +no_defs +type=crs")

dat_bbox <- preds %>% pluck(1) %>% st_as_sf() %>% st_bbox()


preds_all <- preds %>% list_rbind()

# quantiles of predictions for each species
preds_summ <- preds_all %>% 
  dplyr::select(grID,depth_cat,est,species) %>% 
  group_by(species) %>% 
  mutate(quant=cume_dist(est)) %>% 
  ungroup()

# as a simple start, map the number of species above their X% quantile
# CPS first
spp_all <- unique(preds_all$species)
spp_cps <- c("Clupea pallasii","Engraulis mordax","Sardinops sagax","Scomber japonicus","Trachurus symmetricus","Thaleichthys pacificus")
spp_midwater <- c("Leuroglossus stilbius","Stenobrachius leucopsarus","Tactostoma macropus","Tarletonbeania crenularis","Merluccius productus","Sebastes entomelas")

#50% quantile, CPS

cps_quant50 <- preds_summ %>% 
  ungroup() %>% 
  mutate(is_abun=quant>=0.5) %>% 
  filter(is_abun,species%in%spp_cps) %>% 
  group_by(grID,depth_cat) %>% 
  summarise(num_spp=sum(is_abun)) %>% 
  ungroup() %>% 
  left_join(predgrid,by=c("grID","depth_cat")) %>% 
  st_as_sf()

cps_quant50_map <- ggplot()+
  geom_sf(data=coaststates,fill='gray80')+
  geom_sf(data=cps_quant50,aes(fill=factor(num_spp)),color=NA)+
  scale_fill_viridis_d(option = "D")+
  coord_sf(datum=NA)+
  facet_wrap(~depth_cat,nrow=1)+
  xlim(dat_bbox[1],dat_bbox[3])+ylim(dat_bbox[2],dat_bbox[4])+
  labs(fill="Num Species",title="50% quantile abundance or higher\nCoastal Pelagic Species")+
  theme(panel.background = element_rect(fill='white'))
cps_quant50_map

ggsave(here('plots','CPS 50 percent quantile.png'),cps_quant50_map,h=6,w=8,bg='white')

#80% quantile, CPS

cps_quant80 <- preds_summ %>% 
  ungroup() %>% 
  mutate(is_abun=quant>=0.8) %>% 
  filter(is_abun,species%in%spp_cps) %>% 
  group_by(grID,depth_cat) %>% 
  summarise(num_spp=sum(is_abun)) %>% 
  ungroup() %>% 
  left_join(predgrid,by=c("grID","depth_cat")) %>% 
  st_as_sf()

cps_quant80_map <- ggplot()+
  geom_sf(data=coaststates,fill='gray80')+
  geom_sf(data=cps_quant80,aes(fill=factor(num_spp)),color=NA)+
  scale_fill_viridis_d(option = "D")+
  coord_sf(datum=NA)+
  facet_wrap(~depth_cat,nrow=1)+
  xlim(dat_bbox[1],dat_bbox[3])+ylim(dat_bbox[2],dat_bbox[4])+
  labs(fill="Num Species",title="80% quantile abundance or higher\nCoastal Pelagic Species")+
  theme(panel.background = element_rect(fill='white'))
cps_quant80_map

ggsave(here('plots','CPS 80 percent quantile.png'),cps_quant80_map,h=6,w=8,bg='white')

#50% quantile, midwater and DVM species

midwater_quant50 <- preds_summ %>% 
  ungroup() %>% 
  mutate(is_abun=quant>=0.5) %>% 
  filter(is_abun,species%in%spp_midwater) %>% 
  group_by(grID,depth_cat) %>% 
  summarise(num_spp=sum(is_abun)) %>% 
  ungroup() %>% 
  left_join(predgrid,by=c("grID","depth_cat")) %>% 
  st_as_sf()

midwater_quant50_map <- ggplot()+
  geom_sf(data=coaststates,fill='gray80')+
  geom_sf(data=midwater_quant50,aes(fill=factor(num_spp)),color=NA)+
  scale_fill_viridis_d(option = "D")+
  coord_sf(datum=NA)+
  facet_wrap(~depth_cat,nrow=1)+
  xlim(dat_bbox[1],dat_bbox[3])+ylim(dat_bbox[2],dat_bbox[4])+
  labs(fill="Num Species",title="50% quantile abundance or higher\nMidwater Species")+
  theme(panel.background = element_rect(fill='white'))
midwater_quant50_map

ggsave(here('plots','midwater 50 percent quantile.png'),midwater_quant50_map,h=6,w=8,bg='white')

#80% quantile, midwater

midwater_quant80 <- preds_summ %>% 
  ungroup() %>% 
  mutate(is_abun=quant>=0.8) %>% 
  filter(is_abun,species%in%spp_midwater) %>% 
  group_by(grID,depth_cat) %>% 
  summarise(num_spp=sum(is_abun)) %>% 
  ungroup() %>% 
  left_join(predgrid,by=c("grID","depth_cat")) %>% 
  st_as_sf()

midwater_quant80_map <- ggplot()+
  geom_sf(data=coaststates,fill='gray80')+
  geom_sf(data=midwater_quant80,aes(fill=factor(num_spp)),color=NA)+
  scale_fill_viridis_d(option = "D")+
  coord_sf(datum=NA)+
  facet_wrap(~depth_cat,nrow=1)+
  xlim(dat_bbox[1],dat_bbox[3])+ylim(dat_bbox[2],dat_bbox[4])+
  labs(fill="Num Species",title="80% quantile abundance or higher\nMidwater Species")+
  theme(panel.background = element_rect(fill='white'))
midwater_quant80_map

ggsave(here('plots','midwater 80 percent quantile.png'),midwater_quant80_map,h=6,w=8,bg='white')

# Pairs plots (correlations)

# overall logD correlations, all species
logD_corr <- fit %>% 
  mutate(species=ifelse(species=="Zz_Merluccius productus","Merluccius productus",species)) %>% 
  dplyr::select(species,mean,tubeID) %>% 
  pivot_wider(names_from = "species",values_from="mean") %>% 
  dplyr::select(-tubeID) %>% cor()
corrplot(logD_corr,method="ellipse", type="full",
         # addCoef.col='grey80',
         tl.col="black", addrect = 3, order="hclust", hclust.method="ward.D2",
         tl.cex = 0.5)

png(file = here('plots','all species corrplot.png'),type = "cairo",width = 8,height=6,units='in',res=720)
corrplot(logD_corr,method="ellipse", type="full",
         # addCoef.col='grey80',
         tl.col="black", addrect = 3, order="hclust", hclust.method="ward.D2",
         tl.cex = 0.5)
dev.off()

# logD correlations by depth, all species
logD_corr_depth <- fit %>% 
  mutate(species=ifelse(species=="Zz_Merluccius productus","Merluccius productus",species)) %>% 
  dplyr::select(species,depth_cat,mean,tubeID) %>%
  pivot_wider(names_from = "species",values_from="mean") %>% 
  dplyr::select(-tubeID) %>% 
  group_by(depth_cat) %>% 
  nest() %>% 
  mutate(corrs=map(data,cor))

walk2(logD_corr_depth$depth_cat,logD_corr_depth$corrs,function(x,y){
  png(file = here('plots',paste0('all species corr_depth',x,'.png')),type = "cairo",width = 8,height=6,units='in',res=720)
  corrplot(y,method="ellipse", type="full",
           # addCoef.col='grey80',
           tl.col="black", addrect = 3, order="hclust", hclust.method="ward.D2",
           tl.cex = 0.5)
  dev.off()
})

png(file = here('plots','all species corrplot.png'),type = "cairo",width = 8,height=6,units='in',res=720)
corrplot(logD_corr,method="ellipse", type="full",
         # addCoef.col='grey80',
         tl.col="black", addrect = 3, order="hclust", hclust.method="ward.D2",
         tl.cex = 0.5)
dev.off()