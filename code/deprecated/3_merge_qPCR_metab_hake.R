library(here)
library(tidyverse)
source(here("code/1_LoadPackagesData.R"))
source(here("code/2_calibrate_read_proportions.R"))

hake <- read.csv(here("data/hake_qPCR/Hake eDNA 2019 qPCR results 2021-01-04 results.csv")) %>% 
  filter(!str_detect(sample, "ext")) %>% 
  mutate(hake_copies_ul = str_replace(hake_copies_ul, ",", "")) %>% 
  mutate(sample = str_replace(sample, "a", "")) %>% 
  mutate(sample = str_replace(sample, "b", "")) %>% 
  filter(hake_copies_ul != "") %>% 
  select(sample, hake_copies_ul) %>% 
  mutate(hake_copies_ul = as.numeric(hake_copies_ul),
         sample = as.numeric(sample)) %>% 
  filter(!is.na(hake_copies_ul)) %>% 
  filter(!is.na(sample)) %>% 
  left_join(meta, join_by(sample)) %>% 
  filter(!is.na(sampleID))

#QM_bayes_out <- readRDS(here("model_output/QM_bayes_MFU.RDS"))

##just using MEANS as point estimates
hake_means <- hake %>% 
  group_by(sample, depth, lat, lon) %>% 
  summarise(qPCR = mean(hake_copies_ul)) 
  # filter(depth == 0)
  
prop_means <- QM_bayes_out$Bayes_estimates %>% 
  rownames_to_column("sample")
  
fishColumns <- 2:ncol(prop_means) #NOTE annoying hard-coding to tell code which columns are fish data
abs_estimates <- prop_means %>% 
  mutate(sample = as.numeric(sample)) %>% 
  left_join(hake_means) %>% 
  drop_na() %>% 
  rowwise() %>% 
  mutate(totalMiFish = qPCR / `Merluccius productus`) %>% 
  mutate(across(all_of(fishColumns), ~ .x * totalMiFish))

write.csv(abs_estimates, file = here("model_output/abs_estimates_MFU.csv"))

# abs_estimates %>% 
#   filter(depth == 0) %>% 
#   ggplot(aes(x = lon, y = lat, size = `Sardinops sagax`)) +
#     geom_point() +
#   xlab("Longitude") + ylab("Latitude") + ggtitle("Mean Posterior Estimate \n(DNA copies/uL)")
# 
# abs_estimates %>% 
#   filter(depth == 0) %>% 
#   ggplot(aes(x = lon, y = lat, size = `Zz_Engraulis mordax`)) +
#   geom_point() +
#   xlab("Longitude") + ylab("Latitude") + ggtitle("Mean Posterior Estimate \n(DNA copies/uL)")
# 
# abs_estimates %>% 
#   filter(depth == 0) %>% 
#   ggplot(aes(x = lon, y = lat, size = `Trachurus symmetricus`)) +
#   geom_point() +
#   xlab("Longitude") + ylab("Latitude") + ggtitle("Mean Posterior Estimate \n(DNA copies/uL)")
# 
# abs_estimates %>% 
#   filter(depth == 0) %>% 
#   ggplot(aes(x = lon, y = lat, size = `Clupea pallasii`)) +
#   geom_point() +
#   xlab("Longitude") + ylab("Latitude") + ggtitle("Mean Posterior Estimate \n(DNA copies/uL)")
# 
# abs_estimates %>% 
#   filter(depth == 0) %>% 
#   ggplot(aes(x = lon, y = lat, size = log(`Thaleichthys pacificus`))) +
#   geom_point() +
#   xlab("Longitude") + ylab("Latitude") + ggtitle("Mean Posterior Estimate \n(DNA copies/uL)")
# 
(p_mifish_depth <- abs_estimates %>% 
  filter(depth %in% c(0, 50, 150, 300, 500)) %>% 
  select(depth, totalMiFish, lat, lon) %>% 
  #pivot_longer(-depth, names_to = "species", values_to = "conc") %>% 
  ggplot(aes(x = as.factor(depth), y = log(totalMiFish))) +
  geom_boxplot() +
  facet_grid(~as.factor(round(lat, 0))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))


(p_species_depth_latitude <- abs_estimates %>% 
  filter(depth %in% c(0, 50, 150, 300, 500)) %>% 
  select(2:9, lat) %>% 
  pivot_longer(-c(depth, lat), names_to = "species", values_to = "conc") %>% 
  ggplot(aes(x = as.factor(depth), y = log(conc), fill = species)) +
    geom_boxplot() +
    facet_grid(as.factor(round(lat, 0))~as.factor(species)) +
  ylim(c(-3, 12)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))

(p_species_depth <- abs_estimates %>% 
  filter(depth %in% c(0, 50, 150, 300, 500)) %>% 
  select(2:9, lat) %>% 
  pivot_longer(-c(depth, lat), names_to = "species", values_to = "conc") %>% 
  ggplot(aes(x = as.factor(depth), y = log(conc), fill = species)) +
  geom_boxplot() +
  facet_grid(~as.factor(species)) +
  ylim(c(-3, 12)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))

#write.csv(abs_estimates, file = here("model_output/abs_estimates_MFU.csv"))


#write out plots
ggsave(p_mifish_depth, file = here("plots/p_mifish_depth.pdf"), height = 8, width = 10)
ggsave(p_species_depth_latitude, file = here("plots/p_species_depth_latitude.pdf"), height = 8, width = 10)
ggsave(p_species_depth, file = here("plots/p_species_depth.pdf"), height = 8, width = 10)
          
          
          
          
#### Sampling Whole Posterior from QM (and normal approx for qPCR posteriors)
          
          # post.samples <- extract(QM_bayes_out$Bayes_modelfit, par = "int_samp_small")$int_samp_small #[iter, station, species]
          # metabarNames <- row.names(QM_bayes_out$Bayes_estimates) %>% as.numeric
          # 
          # hake_distr_apprx <- hake %>% 
          #   group_by(sample, depth, lat, lon) %>% 
          #   summarise(mean_qPCR = mean(hake_copies_ul),  ##Consider doing this on log scale, to keep positive
          #             sd_qPCR = sd(hake_copies_ul)) %>% 
          #   filter(sample %in% metabarNames)
          # qPCRNames <- hake_distr_apprx$sample
          # keep <- intersect(metabarNames, qPCRNames)
          # 
          # post.samples <- post.samples[,which(metabarNames %in% keep),] 
          # 
          # #check that things are lined up
          # table(metabarNames == hake_distr_apprx$sample)
          # 
          # #total Mifish copies = hake qPCR copies / metabarcoding proportions for `Merluccius productus`
          # #create posterior distribution for total copies for each sample in hand
          # #note minor issue: if we do the normal approx for the distribution of qPCR values, some posterior samples are negative on a linear scale, which doesn't make sense
          # totalCopies <- as.data.frame(matrix(NA, nrow = dim(post.samples)[1], ncol = dim(post.samples)[2]))
          # 
          # for (j in 1:ncol(totalCopies)){
          #   totalCopies[,j] <- rnorm(dim(post.samples)[1],
          #                            hake_distr_apprx$mean_qPCR[j], 
          #                            hake_distr_apprx$sd_qPCR[j]) / 
          #     post.samples[,j,4]
          # }
          #   
          # #now, given the above values for total MiFish copies, allocate those copies among the proportions of species as observed in the calibrated metabarcoding data
          # post.abundance <- as.data.frame(matrix(NA, nrow = dim(post.samples)[2] * dim(post.samples)[3], ncol = 8))
          # #e.g. for species 1, site 1
          # idx <- 1
          # for (i in 1:dim(post.samples)[2]){ #samples
          #   for (j in 1:dim(post.samples)[3]){ #species
          #     post.abundance[idx,1] <- i #sample index
          #     post.abundance[idx,2] <- j #species index
          #     post.abundance[idx,3:8] <- summary((totalCopies[,i] * post.samples[,i,j]))
          #     
          #     idx <- idx + 1
          #   }
          # }
          # names(post.abundance) <- c("Sample", "Species", "Min", "q.25", "q.50", "q.75", "Max")
          # 
          # head(post.abundance)
          