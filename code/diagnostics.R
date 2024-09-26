## Code to plot diagnostics for qPCR/QM joint model
## These are diagnostics for the fitted data, not for the Stan sampler

library(tidyverse)
library(rstan)
library(cowplot)

plot_theme <- theme_minimal()+theme(panel.border = element_rect(color="black",fill=NA),
                                    axis.text=element_text(size=8),plot.title = element_text(size=10))
theme_set(plot_theme)

# pull and organize the summary output from a fitted Stan model object
get_param_ests <- function(m){
  s <- summary(m)$summary 
  param_names <- rownames(s)
  out <- s %>% 
    as_tibble() %>% 
    mutate(variable=param_names) %>% 
    dplyr::select(variable,everything())
  out
}

# x <- get_param_ests(jointMod)

##### qPCR PLOTS ####

plot_obs_pred_qPCR <- function(s,stan_data,return_what="plot"){
  
  # qPCR standards
  qPCR_std_dat <- stan_data %>% pluck('qpcr_std')
  
  qPCR_std_bin_fit <- s %>% 
    filter(grepl('logit_theta_std',variable)) %>% 
    mutate(prob=plogis(mean))
  
  qPCR_std_pos_fit <- s %>% 
    filter(grepl("Ct_std",variable))
  
  qPCR_std_dat <- qPCR_std_dat %>% 
    mutate(prob_amp_fit=qPCR_std_bin_fit$prob,
           Ct_fit=qPCR_std_pos_fit$mean)
  
  # plots 
  qPCR_std_bin_plot <- qPCR_std_dat %>% 
    ggplot(aes(x=log10(copies_ul)))+
    geom_point(aes(y=z))+
    geom_line(aes(y=prob_amp_fit,color=qPCR))+
    labs(x="Log Copies",y="Probability of Amplification",title="qPCR Standards - Bernoulli")+
    guides(color='none')
  qPCR_std_pos_plot <- qPCR_std_dat %>% 
    filter(z==1) %>% 
    ggplot(aes(x=Ct,y=Ct_fit,col=qPCR))+
    geom_point()+
    geom_abline(slope=1,intercept=0,linetype=2)+
    labs(x="Observed Ct",y="Predicted Ct",title="qPCR Standards - Positive")+
    guides(color='none')
  
  # qPCR unknowns
  qPCR_unk_dat <- stan_data %>% pluck('qpcr_unk')
  
  qPCR_unk_bin_fit <- s %>% 
    filter(grepl('logit_theta_samp',variable)) %>% 
    mutate(prob=plogis(mean))
  
  qPCR_unk_pos_fit <- s %>% 
    filter(grepl('Ct',variable)&!grepl('Ct_std',variable))
  
  qPCR_unk_dat <- qPCR_unk_dat %>% 
    mutate(prob_amp_fit=qPCR_unk_bin_fit$prob,
           Ct_fit=qPCR_unk_pos_fit$mean)
  
  # plots
  qPCR_unk_bin_plot <- qPCR_unk_dat %>% 
    ggplot(aes(x=prob_amp_fit))+
    geom_point(aes(y=z))+
    stat_smooth(aes(y=z,x=prob_amp_fit)) +
    geom_abline(intercept=0,slope = 1,linetype="dashed",color="red") +
    scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1))+
    labs(x="Predicted Probability of Amplification",y="Observed Amplification",title="qPCR Unknowns - Bernoulli")+
    guides(color='none')
  qPCR_unk_pos_plot <- qPCR_unk_dat %>% 
    filter(z==1) %>% 
    filter(Ct_fit<50) %>% 
    ggplot(aes(Ct,Ct_fit,col=qPCR))+
    geom_point()+
    geom_abline(slope=1,intercept=0,linetype=2)+
    labs(x="Observed Ct",y="Predicted Ct",title="qPCR Unknowns - Positive")+
    guides(color='none')
  
  if(return_what=='data'){
    list(qPCR_unk=qPCR_unk_dat,
         qPCR_std=qPCR_std_dat)
  }else{
    # combined plots
    plot_grid(qPCR_std_bin_plot,qPCR_std_pos_plot,qPCR_unk_bin_plot,qPCR_unk_pos_plot,nrow=2)  
  }
  

}

#### alphas (metabarcoding efficiency) ####

plot_alphas <- function(s){
  
  # estimated alphas
  d <- s %>% 
    filter(grepl('alpha',variable),!grepl('alpha_raw',variable),!grepl('dm_alpha',variable))
  
  # plot
  log_eff_plot <- d %>% 
    # ggplot(aes(variable,mean,ymax=`97.5%`,ymin=`2.5%`))+
    ggplot(aes(variable,mean,ymax=mean+se_mean,ymin=mean-se_mean))+
    geom_pointrange()+
    geom_hline(yintercept=0)+
    labs(y="estimated log-efficiency",x="")
  # rel_eff_plot <- d %>% 
  #   mutate(releff=exp(mean),
  #          upper=exp(mean+se_mean),
  #          lower=exp(mean-se_mean)) %>% 
  #   # ggplot(aes(variable,mean,ymax=`97.5%`,ymin=`2.5%`))+
  #   ggplot(aes(variable,releff,ymax=upper,ymin=lower))+
  #   geom_pointrange()+
  #   geom_hline(yintercept=1)+
  #   labs(y="estimated relative efficiency",x="")
  # 
  # plot_grid(log_eff_plot,rel_eff_plot)
  log_eff_plot
}

#### predictions vs. observations, metabarcoding mock communities and field samples ####

plot_obs_pred_qm <- function(s,stan_data){
  
  # mock communities - reads
  mock_reads_dat <- stan_data %>%
    pluck('mock_data') %>% 
    mutate(rep=row_number()) %>% 
    pivot_longer(-rep,names_to = 'taxa',values_to='reads')
  
  # mock communities - proportions
  mock_alr_dat <- stan_data %>%
    pluck('alr_mock_true_prop') %>% 
    mutate(rep=row_number()) %>% 
    pivot_longer(-rep,names_to = 'taxa',values_to='alr')%>% 
    group_by(rep) %>% 
    mutate(prop=softmax(alr)) %>% 
    ungroup()
  
  # mock communities - fitted (estimated) values
  mock_fit <- s %>% 
    filter(grepl('logit_val_mock',variable)) %>% 
    # number of replicates (for doing the softmax) is equal to nrow(mock_dat)
    # you can also pull it from the names of the variables with regex
    mutate(rep=str_extract(variable,"(?<=\\[)\\d")) %>% 
    mutate(se_mean=ifelse(is.nan(se_mean),0,se_mean)) %>% 
    mutate(alr_upper=mean+se_mean,alr_lower=mean-se_mean) %>% 
    group_by(rep) %>% 
    mutate(prop=softmax(mean),prop_upper=softmax(mean+se_mean),prop_lower=softmax(mean-se_mean)) %>% 
    ungroup()
  
  # join obs/preds
  mock_alr_join <- mock_alr_dat %>% 
    mutate(alr_fit=mock_fit$mean,
           alr_upper=mock_fit$alr_upper,
           alr_lower=mock_fit$alr_lower,
           prop_fit=mock_fit$prop,
           prop_upper=mock_fit$prop_upper,
           prop_lower=mock_fit$prop_lower)
  
  mock_reads_join <- mock_reads_dat %>% 
    mutate(prop_fit=mock_fit$prop,
           prop_upper=mock_fit$prop_upper,
           prop_lower=mock_fit$prop_lower) %>% 
    # back out estimated reads from proportions
    group_by(rep) %>% 
    mutate(sumreads=sum(reads)) %>% 
    ungroup() %>% 
    mutate(reads_fit=prop_fit*sumreads,
           reads_upper=prop_upper*sumreads,
           reads_lower=prop_lower*sumreads)
  
  # mock communities - plots
  mock_alr_plot <- mock_alr_join %>% 
    ggplot(aes(alr_fit,alr,col=taxa,xmin=alr_lower,xmax=alr_upper))+
    geom_pointrange()+
    geom_abline(slope=1,intercept=0,linetype=2)+
    labs(x="ALR- Predicted",y="ALR - Observed",title="Mock Communities - ALR")+
    guides(color='none')
  
  mock_prop_plot <- mock_alr_join %>% 
    ggplot(aes(prop_fit,prop,col=taxa,xmin=prop_lower,xmax=prop_upper))+
    geom_pointrange()+
    geom_abline(slope=1,intercept=0,linetype=2)+
    labs(x="Proportion - Predicted",y="Proportion - Observed",title="Mock Communities - ALR")+
    guides(color='none')
  
  mock_reads_plot <- mock_reads_join %>% 
    ggplot(aes(reads_fit,reads,col=taxa,xmin=reads_lower,xmax=reads_upper))+
    geom_pointrange()+
    geom_abline(slope=1,intercept=0,linetype=2)+
    labs(x="Reads- Predicted",y="Reads - Observed",title="Mock Communities - Reads")+
    guides(color='none')
  
  # metabarcoding field samples - reads
  mb_unk_reads_dat <- stan_data %>%
    pluck('sample_data') %>% 
    mutate(rep=row_number()) %>% 
    pivot_longer(-rep,names_to = 'taxa',values_to='reads') %>% 
    # add total reads
    group_by(rep) %>% 
    mutate(sumreads=sum(reads),
           prop=reads/sumreads) %>% 
    ungroup()
  
  # metabarcoding field samples - fitted (estimated) reads
  mb_unk_reads_fit <- s %>% 
    filter(grepl('logit_val_samp',variable)) %>% 
    # extract a replicate number from the variable names using regex
    mutate(rep=str_extract(variable,"(?<=\\[)\\d+") %>% as.integer) %>% 
    mutate(se_mean=ifelse(is.nan(se_mean),0,se_mean)) %>% 
    group_by(rep) %>% 
    mutate(prop=softmax(mean),prop_upper=softmax(mean+se_mean),prop_lower=softmax(mean-se_mean)) %>% 
    ungroup()
  
  mb_unk_join <- mb_unk_reads_dat %>% 
    mutate(prop_fit=mb_unk_reads_fit$prop,
           prop_upper=mb_unk_reads_fit$prop_upper,
           prop_lower=mb_unk_reads_fit$prop_lower) %>% 
    # back out estimated reads from proportions
    mutate(reads_fit=prop_fit*sumreads,
           reads_upper=prop_upper*sumreads,
           reads_lower=prop_lower*sumreads)
  
  # metabarcoding field samples - plots
  
  mb_unk_prop_plot <- mb_unk_join %>%
    ggplot(aes(prop_fit,prop,col=taxa,xmin=prop_lower,xmax=prop_upper))+
    geom_pointrange()+
    geom_abline(slope=1,intercept=0,linetype=2)+
    labs(x="Proportion - Predicted",y="Proportion - Observed",title="Metabarcoding Field Samples - Proportions")+
    guides(color='none')
  mb_unk_reads_plot <- mb_unk_join %>% 
    # log 10 scaling:
    ggplot(aes(log10(reads_fit+1),log10(reads+1),col=taxa,xmin=log10(reads_lower+1),xmax=log10(reads_upper+1)))+
    # untransformed scaling:
    # ggplot(aes(reads_fit,reads,col=taxa,xmin=reads_lower,xmax=reads_upper))+
    geom_pointrange()+
    geom_abline(slope=1,intercept=0,linetype=2)+
    labs(x="Log10 Reads- Predicted",y="Log10 Reads - Observed",title="Metabarcoding Field Samples - Reads")
  
  # Combined plots
  mockp <- plot_grid(mock_alr_plot,mock_prop_plot,mock_reads_plot,nrow=1)
  mbunkp <- plot_grid(mb_unk_prop_plot,mb_unk_reads_plot,nrow=1)
  out <- plot_grid(mockp,mbunkp,nrow=2)
  
  out
}

#### distribution of log_D (estimated reads, i.e., our quantity of interest) ####

plot_est_reads <- function(s){
  
  mb_unk_reads_fit <- s %>% 
    filter(grepl('log_D',variable),
           !grepl('log_D_raw',variable),
           !grepl('log_D_station_depth',variable),
           !grepl('log_D_sp',variable),
           !grepl('log_D_mu',variable),
           !grepl('log_D_sigma',variable),
           !grepl('log_D_sigma',variable)) %>%
    # extract a replicate number from the variable names using regex
    mutate(rep=str_extract(variable,"(?<=\\[)\\d+") %>% as.integer)%>% 
    # extract a species number from the variable names using regex
    mutate(taxa=str_extract(variable,"(?<=\\,)\\d+") %>% as.integer)
  
  # plot 
  mb_unk_reads_plot <- mb_unk_reads_fit %>% 
    ggplot(aes(mean,fill=factor(taxa)))+
    geom_density()+
    facet_wrap(~taxa)+
    labs(x="Estimated Log Copies",y="Kernel Density",title="Field Samples - Distribution of Estimated Reads",fill="Species")
  
  mb_unk_reads_plot
}

#### Other variance parameters ####
