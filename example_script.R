### example of how to compare estimates of population trends while accounting for the precision of those trend estimates

library(tidyverse)
library(jagsUI)
library(ggmcmc)
library(tidybayes)



# reading in trends ---------------------------------------

tr = read.csv("data/Trend_Estimates.csv")

### funciton to transform %-change trends to log-scale geometric mean rates of change
log_trans <- function(x){
  log((x/100)+1)
}
## transforming the 95% CIs into an estimate of variance
log_var_func <- function(uci,lci){
  ((log_trans(uci)-log_trans(lci))/(1.96*2))^2
}

data_prep <- function(tr_dat,t1, t2){
  cls = paste(rep(c("annual_trend",
                    "lower",
                    "upper"),times = 2),rep(c(t1,t2),each = 3),sep = "_")
 
  tmp_dat = tr_dat[,c("species",cls)]
  names(tmp_dat) <- c("species","t1","t1lci","t1uci","t2","t2lci","t2uci")
  tmp_dat <- tmp_dat %>% drop_na() %>% 
    mutate(betahat1 = log_trans(t1),
           varhat1 = log_var_func(t1uci,t1lci),
           betahat2 = log_trans(t2),
           varhat2 = log_var_func(t2uci,t2lci))
  
  nspecies = nrow(tmp_dat)
  betahat = as.matrix(tmp_dat[,c("betahat1","betahat2")])
  varhat = as.matrix(tmp_dat[,c("varhat1","varhat1")])
  
  out_data <- list(nspecies = nspecies,
                    betahat = betahat,
                    varhat = varhat)
  return(list(jags_data = out_data,
                  full_data = tmp_dat))
}


comparisons = read.csv("data/Comparisons.csv")


# model -------------------------------------------------------------------

modl <- "

### model for comparing trend estimates from two (indexed by e) different analyses, for different species (indexed by s)
### Model is modified from the grassland bird model from Link and Barker 2010, pg 90, also described in Sauer and Link 2002, Ecology 83:1743-1751
### and includes covariance components and corelation coefficient aspects from Matzke et al. 2017, https://doi.org/10.1525/collabra.78

### in essence, this is a paired t-test style comparison, that accounts for the imprecision in each s*e trend estimate.

## varhat = estaimtes of the variance of the trend estimates
## betahat = estimates of the trends
## nspecies = number of species

model{

for (e in 1:2) {

for(s in 1:nspecies) {
  tau.betahat[s,e] <- 1/varhat[s,e] #transform variance to precision
	betahat[s,e] ~ dnorm(beta[s,e],tau.betahat[s,e] ) #betahat = data = trend estimates, tau.betahat = data = precision of trend estimate

} # end of s loop

mu[e] ~ dnorm(0,1) #mean trend for each survey
tau.beta[e] ~ dscaled.gamma(1,10) 
sd_beta[e] <- 1/sqrt(tau.beta[e])


	} #end of e loop (indexing two models being compared)


#### multivariate normal structure for the among-strata between survey variation in trends

for(s in 1:nspecies) {

beta[s,1:2] ~ dmnorm(mu[1:2],ISigma_cov[1:2,1:2])
}

#rho ~ dnorm(0,0.3)T(-1,1) #estimated Pearson correlation coefficient - prior truncated normal
rho ~ dunif(-1,1)
### sigma_cov is the covariance matrix
Sigma_cov[1,1] <- pow(sd_beta[1],2)
Sigma_cov[1,2] <- rho*sd_beta[1]*sd_beta[2]
Sigma_cov[2,1] <- rho*sd_beta[1]*sd_beta[2]
Sigma_cov[2,2] <- pow(sd_beta[2],2)
ISigma_cov[1:2,1:2] <- inverse(Sigma_cov[1:2,1:2])


for(s in 1:nspecies) {
dif[s] <- beta[s,1]-beta[s,2] # dif is a vector of the strata-specific trend differences after accounting for the imprecision of each estimate's trend and the group (survey/monitoring program) structure
} # end of second s-strata loop


m_dif <- mean(dif[])
mu_dif <- mu[1]-mu[2]
p_pos <- step(mu_dif)

} # end of model

   "  

trend_comp = "models/trend_comparison_correlation.txt"
cat(modl,file = trend_comp)   

betas_out <- vector(mode = "list",length = nrow(comparisons))
difs_out <- betas_out
summr <- betas_out
lines_out <- betas_out

for(ic in 2:3){
  
  trend1 = comparisons[ic,"t1"]
  trend2 = comparisons[ic,"t2"]

  tmp_data <- data_prep(tr,trend1,trend2)



params <- c("m_dif",
            "mu_dif",
            "beta",
            "dif",
            "rho",
            "sd_beta",
            "mu",
            "p_pos")

jags_data <- tmp_data$jags_data

burnInSteps = 2000            # Number of steps to "burn-in" the samplers. this is sufficient for testing, but you'll want to increase this
nChains = 3                   # Number of chains to run.
numSavedSteps=1000         # Total number of steps in each chain to save. this is sufficient for testing, but you'll want to increase this
thinSteps=40                   # Number of steps to "thin" (1=keep every step).
nIter = ceiling( ( (numSavedSteps * thinSteps )+burnInSteps)) # Steps per chain.



out = jagsUI(data = jags_data,
             parameters.to.save = params,
             n.chains = 3,
             n.burnin = burnInSteps,
             n.thin = thinSteps,
             n.iter = nIter,
             parallel = T,
             model.file = trend_comp)


torig = data.frame(p = c(1,2),
                   trend = c(trend1,trend2))
sporig = data.frame(species = tmp_data$full_data$species,
                    s = 1:jags_data$nspecies)

betasamples = gather_draws(out$samples,beta[s,p])
betas <- betasamples %>% group_by(s,p) %>% 
  summarise(mean = mean((exp(.value)-1)*100),
            lci = quantile((exp(.value)-1)*100,0.025),
            uci = quantile((exp(.value)-1)*100,0.975)) %>% 
  left_join(.,torig,by = "p") %>% 
  left_join(.,sporig,by = "s") %>% 
  ungroup() %>% 
  select(species,trend,mean,lci,uci) %>% 
  pivot_wider(names_from = trend,
              values_from = c(mean,lci,uci)) %>% 
  relocate(c(1,2,4,6,3,5,7)) %>% 
  mutate(comparison = paste(trend1,trend2,sep = "-"),
         version = "modeled")

orig <- tmp_data$full_data %>% 
  select(species,t1,t1lci,t1uci,
         t2,t2lci,t2uci) %>% 
  mutate(comparison = paste(trend1,trend2,sep = "-"),
         version = "original")
names(orig) <- names(betas)

betas <- bind_rows(betas,orig)




betas_out[[ic]] <- betas

rho_samples <- gather_draws(out$samples,rho,p_pos)
rho <- rho_samples %>% 
  summarise(mean = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975))




dif_samples <- gather_draws(out$samples,mu_dif,mu[p])


difs <- dif_samples %>% 
  summarise(mean = mean((exp(.value)-1)*100),
            lci = quantile((exp(.value)-1)*100,0.025),
            uci = quantile((exp(.value)-1)*100,0.975)) %>% 
  full_join(.,torig,by = "p") %>% 
  bind_rows(.,rho) %>% 
  mutate(comparison = paste(trend1,trend2,sep = "-"))

difs_out[[ic]] <- difs 
  
summr[[ic]] <- out$summary #table showing a summary of the posterior distribution and some basic convergence stats for all the monitored parameters

rho_samples <- rho_samples %>% filter(.variable == "rho") %>% 
  ungroup() %>% 
  select(.draw,.value) %>% 
  rename(rho = .value)


mum_samples <- gather_draws(out$samples,mu_dif) %>% 
  ungroup() %>% 
  select(.draw,.value) %>% 
  rename(mum = .value)


line_samples <- left_join(rho_samples,mum_samples,by = ".draw")

lines_out[[ic]] <- line_samples

}


all_difs <- NULL
beta_plot<- NULL
for(ic in 2:3){
  tmp <- difs_out[[ic]]
  all_difs <- bind_rows(all_difs,tmp)
  tmp2 = betas_out[[ic]]
  names(tmp2)[2:4] <- c("mean_Liard",
                   "lci_Liard",
                   "uci_Liard")
  beta_plot <- bind_rows(beta_plot,tmp2)
}
# plotting correlations ---------------------------------------------------

rhos <- all_difs %>% filter(.variable == "rho",
                            grepl(pattern = "BBS",x = comparison))

rho_plot <- ggplot(data = rhos,aes(x = comparison,y = mean))+
  geom_pointrange(aes(ymin = lci,ymax = uci))+
  xlab("")+
  ylab("Pearson Correlation Coefficient")+
  geom_hline(yintercept = 0,alpha = 0.2)+
  coord_flip()+
  theme_classic()

pdf(file = "figures/Correlation_Coefficients.pdf")
print(rho_plot)
dev.off()

# plotting differences ---------------------------------------------------
difs <- all_difs %>% filter(.variable == "mu_dif",
                            grepl(pattern = "BBS",x = comparison))

dif_plot <- ggplot(data = difs,aes(x = comparison,y = mean))+
  geom_pointrange(aes(ymin = lci,ymax = uci))+
  xlab("")+
  ylab("Modeled Mean Difference in Trends (%/year)")+
  geom_hline(yintercept = 0,alpha = 0.2)+
  coord_flip()+
  theme_classic()

pdf(file = "figures/Mean_differences_between_trends.pdf")
print(dif_plot)
dev.off()

nms = read.csv("data/Parameter_names.csv")
all_difs <- left_join(all_difs,nms,by = ".variable")
write.csv(all_difs,file = "parameters_estimates.csv",row.names = FALSE)




# plotting the betas ------------------------------------------------------
lbls <- c("Overall Liard Trend","After removing habitat effects")
lab_df <- data.frame(comparison_labs = lbls,
                     comparison = c("StudyArea-BBS","Regional-BBS"))

lab_df$comparison_labs <- factor(lab_df$comparison_labs,levels = lbls,
       ordered = TRUE)

beta_plot <- left_join(beta_plot,lab_df,by = "comparison")

bt = beta_plot %>% #filter(comparison == "StudyArea-BBS") %>% 
  mutate(sd_Liard = uci_Liard-lci_Liard,
         sd_BBS = uci_BBS-lci_BBS,
         prec = 1/(sd_Liard^2 + sd_BBS^2))
bp = ggplot(data = bt,aes(x = mean_BBS,y = mean_Liard))+
  geom_point(aes(colour = version,size = prec),alpha = 0.6)+
  geom_line(aes(group = species),alpha = 0.1,size = 0.3)+
  xlab("BBS Trend")+
  ylab("Liard Trend")+
  scale_colour_viridis_d(aesthetics = c("colour","fill"),
                         begin = 0.4,
                         end = 0.95,
                         option = "mako")+
  theme_classic()+
  facet_wrap(~comparison_labs,ncol = 2,nrow = 1)

print(bp)



pdf(file = "Figures/Liard_BBS_comparison.pdf"
  ,width = 8.5,
    height = 6)
print(bp)
dev.off()


# 
# ic = 1
# trend1 = comparisons[ic,"t1"]
# trend2 = comparisons[ic,"t2"]
# 
# plot_dat <- betas_out[[ic]]
# names(plot_dat) <- gsub(names(plot_dat),pattern = trend1,replacement = "1")
# names(plot_dat) <- gsub(names(plot_dat),pattern = trend2,replacement = "2")
# 
# prec_func = function(l,u){
#   1/(((u-l)/4)^2)
# }
# plot_dat <- plot_dat %>% 
#   mutate(size1 = prec_func(lci_1,uci_1),
#          size2 = prec_func(lci_2,uci_2),
#          sizem = (size1+size2)/2)
# 
# difs <- difs_out[[ic]]
# 
# # mns = unlist(c(plot_dat[which(plot_dat$version == "modeled"),c("mean_1")],plot_dat[which(plot_dat$version == "modeled"),c("mean_2")]))
# # 
# # s = min(mns)
# # e = max(mns)
# # 
# # lines_plot <- lines_out[[ic]] 
# # predl <- matrix(seq(from = s,to = e,length = 100),ncol = 100,nrow = nrow(lines_plot),byrow = TRUE)
# 
#  sl = as.numeric(difs[which(difs$.variable == "rho"),"mean"]) 
#  mu =  as.numeric(difs[which(difs$.variable == "mu_dif"),"mean"]) 
#   #sample_n(100)
#  line_dat1 <- plot_dat %>% select(species,mean_1,mean_2,version) %>% 
#    pivot_wider(values_from = c(mean_1,mean_2),
#                                       names_from = version)
# 
# trend_plot <- ggplot(data = plot_dat,aes(y = mean_1,x = mean_2,colour = version))+
#   geom_point(aes(size = sizem))+
#   #geom_errorbar(aes(ymin = lci_1,ymax = uci_1),alpha = 0.1)+
#   #geom_errorbarh(aes(xmin = lci_2,xmax = uci_2),alpha = 0.1)+
#   geom_abline(intercept = mu,slope = sl,alpha = 0.1)+
#   coord_cartesian(xlim = c(-10,10),ylim = c(-10,10))+
#   theme_minimal()+
#   xlab(trend2)+
#   ylab(trend1)
#   
# print(trend_plot)
# 









