# Compare outputs between different joint model flavors
library(tidyverse)
dat1 <- read_rds(here('code','log_D_fitted.rds')) %>% 
  dplyr::select(variable,mean,se_mean,taxa:wash_idx) %>% 
  rename(dat1mean=mean,dat1se=se_mean)
glimpse(dat1)
dat2 <- read_rds(here('code',"2024 multispecies analyses & outputs (of 2019 data)",'log_D_fitted_Ole_all_hake_mocks_fixed.rds')) %>%  
  dplyr::select(variable,mean,se_mean,taxa:wash_idx) %>% 
  rename(dat2mean=mean,dat2se=se_mean)
glimpse(dat2)

datboth <- dat1 %>%
  left_join(dat2)
glimpse(datboth)

scatter_preds <- datboth %>% 
  ggplot(aes(dat1mean,dat2mean))+
  geom_point()+
  geom_abline(slope=1,intercept=0,linetype=2)+
  labs(x="Owen",y="Ole")
scatter_preds
# seems to be in the lower end (the zeroes) where the difference are, with 

scatter_resids <- datboth %>% 
  mutate(diff=dat1mean-dat2mean) %>% 
  ggplot(aes(dat1mean,diff))+
  geom_point()+
  labs(x='log conc.',y="Owen - Ole")
scatter_resids
# same pattern here- all of the difference is in the estimation of zeroes

# What about the alphas?
# mod1 <- read_rds(here('code','fit_dm.rds')) %>% extract()
# mod1$alpha %>% summary()
# 
# mod2 <- read_rds(here('code','2024 multispecies analyses & outputs (of 2019 data)','fit_dm_Ole_no_pure_hake_mocks_fixed.rds')) %>% extract()
# mod2$alpha %>% summary()