library(bayesplot)
library(caret)
library(coda)
library(corrplot)
library(ggcorrplot)
library(grid)
library(gridExtra)
library(glmnet)
library(mvtnorm)
library(nimble)
library(tidyverse)
library(sf)
library(rgdal)
library(runjags)
library(viridis)

# Functions -----

ditch_the_axes <- theme(
  axis.text = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_blank())

# Load data -----
set.seed(432)

predictors <- read.csv("data/all_campaigns_predictors.csv") 
gas_vocs <- read.csv("data/gasoline_vocs.csv")

# Data manipulation -----

hexane_data <- gas_vocs %>% filter(type == "Hexane") %>% 
  mutate(log_concentration = log(concentration)) %>% 
  select(c("ID",  "X_km", "Y_km", "campaign", "log_concentration")) %>%  
  spread(key = "campaign", value = "log_concentration")  %>% 
  select( "ID", "X_km", "Y_km", "Dec 2005", "Apr 2006", "Aug 2006") %>% 
  na.omit()

voc_predictors <- merge(hexane_data,
                        predictors, by = "ID")

voc_predictors[is.na(voc_predictors)] <- 0

# Hexane predictors -----

hexane_predictors <- voc_predictors %>% select("ID", "X_km", "Y_km", 
                                                      "Dec 2005","Apr 2006", 
                                                      "Aug 2006", 
                                                      "Residential_1000m",
                                                      "Residential_100m",
                                                      "Residential_500m",
                                                      "Roads_100m",
                                                      "Roads_200m",
                                                      "Building_1000m",
                                                      "Pop_500m",
                                                      "Commercial_1000m",
                                                      "Resource.and.Industrial_1000m",
                                                      "Government.and.Institutional_200m",
                                                      "Government.and.Institutional_1000m",
                                                      "Waterbody_1000m",
                                                      "Av_NOx_50m",
                                                      "Av_NOx_500m",
                                                      "Av_NOx_1000m")

hexane_predictors[,c("X_km", "Y_km")] <- scale(hexane_predictors[,c("X_km", "Y_km")])
dist_matrix <- dist(hexane_data[c("X_km", "Y_km")]) %>% as.matrix()

campaigns_indx <- which(colnames(hexane_predictors) %in% c("ID","Dec 2005", "Apr 2006", "Aug 2006"))
hexane_predictors[,-campaigns_indx] <- scale(hexane_predictors[,-campaigns_indx])

# Define variables for program -----

Nsites <- nrow(final_predictors_vocs)

ones <- rep(1,Nsites) %>% as.data.frame()

D1 <- c(0,1,0)
D2 <- c(0,0,1)

covs <- data.frame(ones, hexane_predictors[,-c(1,4:6)])

D <- data.frame(D1, D2) %>% as.matrix()

# Model 1 ------

model1 <- nimbleCode(
  {
    
    for(i in 1:(Nsites-1)){
      for(j in (i+1):Nsites){
        for(k in 1:Ngroups){
          C[k,i,j] <- sigma2[k]*exp((-1)*dist_matrix[i,j]/phi[k])
          C[k,j,i] <- C[k,i,j];
        }
      }
    }
    
    for(i in 1:Nsites){
      for(k in 1:Ngroups){
        C[k,i,i] <- sigma2[k] + tau2
      }
    }
    
    for(i in 1:Ngroups){
      mu[,i] <-  rep(dummy[i,]%*%alpha[1:(Ngroups-1)], Nsites) + covs[,1:p]%*%beta[1:p] 
      y[,i] ~ dmnorm(mu[,i], cov = C[i,,])
      fitted[,i]~ dmnorm(mu[,i], cov = C[i,,])
    }
    
    
    # Priors
    for(k in 1:p){
      beta[k] ~ dnorm(0, sd = 10)
    }
    
    for(j in 1:(Ngroups-1)){
      alpha[j] ~ dnorm(0, sd = 10)
    }
    
    for(h in 1:Ngroups){
      sigma2[h] ~ dinvgamma(2,1)
      phi[h] ~ dexp(0.3)
      
    }
    
    tau2 ~ dinvgamma(2,1)
    
  })


spatialCodeinits <-list( "alpha" = rnorm(2),
                         "beta" = rnorm(ncol(covs)), 
                         "phi" = rexp(3,0.3), 
                         "tau2" = rinvgamma(1,2,1),
                         "sigma2" = rinvgamma(3,2,1))


model <- nimbleModel(model1, data = list(y = final_predictors_vocs[,4:6] %>% as.matrix(), 
                                                   covs = covs %>% as.matrix(),
                                                   dist_matrix = hexane_dist,
                                                   dummy = D),
                     constants = list( Nsites = Nsites, 
                                       Ngroups = 3,
                                       p = ncol(covs)), 
                     dimensions = list(mu = c(Nsites,3),
                                       C = c(3,Nsites, Nsites),
                                       fitted = c(nrow(final_predictors_vocs),3)),
                     inits = spatialCodeinits)

Cmodel <- compileNimble(model,showCompilerOutput = TRUE)

start_time <- Sys.time()

mcmcConf <- configureMCMC(Cmodel,
                          monitors = c("alpha","beta", 
                                       "tau2", "phi", 
                                       "sigma2", "fitted"))

modelMCMC <- buildMCMC(mcmcConf,
                       enableWAIC = TRUE)

CmodelMCMC <- compileNimble(modelMCMC, project = model, showCompilerOutput = TRUE)

mcmc.out <- runMCMC(CmodelMCMC, niter = 30000,
                    nburnin = 6000, thin = 14, nchains = 2,
                    summary = TRUE, inits = spatialCodeinits, WAIC = TRUE)

end_time <- Sys.time()
end_time - start_time

coda_sample_1 <- coda::mcmc(mcmc.out$samples$chain1) 
coda_sample_2 <- coda::mcmc(mcmc.out$samples$chain2) 
coda_sample <- coda::mcmc.list(coda_sample_1, coda_sample_2)

n <- ncol(coda_sample_1) #numero de parametros

Rhat_all <-  sapply(1:n, function(x) {
  rstan::Rhat(cbind(coda_sample_1[, x], coda_sample_2[, x]))
})

ess_all <- sapply(1:n, function(x) {
  rstan::ess_bulk(cbind(coda_sample_1[, x], coda_sample_2[, x]))
})

ess_all_tail <- sapply(1:n, function(x) {
  rstan::ess_tail(cbind(coda_sample_1[, x], coda_sample_2[, x]))
})

max(Rhat_all)
min(ess_all)
min(ess_all_tail)

# Model 2 -----

model2 <- nimbleCode(
  {
    
    for(i in 1:Ngroups){
      C[i,,] <-  diag(rep(tau2,Nsites))
      mu[,i] <-  rep(dummy[i,]%*%alpha[1:(Ngroups-1)], Nsites) + covs[,1:p]%*%beta[1:p] 
      y[,i] ~ dmnorm(mu[,i], cov = C[i,,])
      fitted[,i]~ dmnorm(mu[,i], cov = C[i,,])
    }
    
    
    # Priors
    for(k in 1:p){
      beta[k] ~ dnorm(0, sd = 10)
    }
    
    for(j in 1:(Ngroups-1)){
      alpha[j] ~ dnorm(0, sd = 10)
    }
    
    tau2 ~ dinvgamma(2,1)
    
  })


spatialCodeinits <-list( "alpha" = rnorm(2),
                         "beta" = rnorm(ncol(covs)), 
                         "tau2" = rinvgamma(1,2,1))


model <- nimbleModel(model2, data = list(y = final_predictors_vocs[,4:6] %>% as.matrix(), 
                                                    covs = covs %>% as.matrix(),
                                                    dummy = D),
                     constants = list( Nsites = Nsites, 
                                       Ngroups = 3,
                                       p = ncol(covs)), 
                     dimensions = list(mu = c(Nsites,3),
                                       C = c(3,Nsites, Nsites),
                                       fitted = c(nrow(final_predictors_vocs),3)),
                     inits = spatialCodeinits)

Cmodel <- compileNimble(model,showCompilerOutput = TRUE)

start_time <- Sys.time()

mcmcConf <- configureMCMC(Cmodel,
                          monitors = c("alpha","beta", 
                                       "tau2", "fitted"))

modelMCMC <- buildMCMC(mcmcConf,
                       enableWAIC = TRUE)

CmodelMCMC <- compileNimble(modelMCMC, project = model, showCompilerOutput = TRUE)

mcmc.out2 <- runMCMC(CmodelMCMC, niter = 30000,
                     nburnin = 6000, thin = 14, nchains = 2,
                     summary = TRUE, inits = spatialCodeinits, WAIC = TRUE)

end_time <- Sys.time()
end_time - start_time

# Model 3 -----

model3 <- nimbleCode(
  {
    
    for(i in 1:(Nsites-1)){
      for(j in (i+1):Nsites){
        for(k in 1:Ngroups){
          C[k,i,j] <- sigma2[k]*exp((-1)*dist_matrix[i,j]/phi[k])
          C[k,j,i] <- C[k,i,j];
        }
      }
    }
    
    for(i in 1:Nsites){
      for(k in 1:Ngroups){
        C[k,i,i] <- sigma2[k] + tau2
      }
    }
    
    for(i in 1:Ngroups){
      mu[,i] <-  covs[,1:p]%*%gamma[i,1:p] 
      y[,i] ~ dmnorm(mu[,i], cov =  C[i,,] )
      fitted[,i]~ dmnorm(mu[,i], cov = C[i,,])
    }
    
    # Priors
    for(k in 1:p){
      beta[k] ~ dnorm(0, sd = 10)
      for(j in 1:Ngroups){
        gamma[j,k] ~ dnorm(beta[k], var = psi2)
        
      }
    }
    
    
    for(h in 1:Ngroups){
      sigma2[h] ~ dinvgamma(2,1)
      phi[h] ~ dexp(0.3)
      
    }
    
    tau2 ~ dinvgamma(2,1)
    psi2 ~ dinvgamma(2,1)
    
  })

spatialCodeinits <-list("beta" = rnorm(ncol(covs)), 
                        "tau2" = rinvgamma(1,2,1))


model <- nimbleModel(model3, 
                     data = list(y = final_predictors_vocs[,4:6] %>% as.matrix(), 
                                 covs = covs %>% as.matrix(), 
                                 dist_matrix = hexane_dist),
                     constants = list(Nsites = nrow(final_predictors_vocs), 
                                      Ngroups = 3,
                                      p = ncol(covs)), 
                     dimensions = list(mu = c(nrow(final_predictors_vocs),3),
                                       C = c(3,nrow(final_predictors_vocs), 
                                             nrow(final_predictors_vocs)),
                                       fitted = c(nrow(final_predictors_vocs),3)),
                     inits = spatialCodeinits)

Cmodel <- compileNimble(model)

start_time <- Sys.time()

mcmcConf <- configureMCMC(Cmodel,
                          monitors = c("gamma","beta", 
                                       "tau2", "phi", 
                                       "sigma2", "fitted" ))


modelMCMC <- buildMCMC(mcmcConf,
                       enableWAIC = TRUE)

CmodelMCMC <- compileNimble(modelMCMC, project = model)

mcmc.out3 <- runMCMC(CmodelMCMC, niter = 30000,
                     nburnin = 6000, thin = 14, nchains = 2,
                     summary = TRUE, inits = spatialCodeinits, WAIC = TRUE)
end_time <- Sys.time()
end_time - start_time

coda_sample_1 <- coda::mcmc(mcmc.out3$samples$chain1) 
coda_sample_2 <- coda::mcmc(mcmc.out3$samples$chain2) 
coda_sample <- coda::mcmc.list(coda_sample_1, coda_sample_2)

n <- ncol(coda_sample_1) #numero de parametros

Rhat_all <-  sapply(1:n, function(x) {
  rstan::Rhat(cbind(coda_sample_1[, x], coda_sample_2[, x]))
})

ess_all <- sapply(1:n, function(x) {
  rstan::ess_bulk(cbind(coda_sample_1[, x], coda_sample_2[, x]))
})

ess_all_tail <- sapply(1:n, function(x) {
  rstan::ess_tail(cbind(coda_sample_1[, x], coda_sample_2[, x]))
})

max(Rhat_all)
min(ess_all)
min(ess_all_tail)

# Model 4 -----

model4 <- nimbleCode(
  {
    
    for(i in 1:Ngroups){
      mu[,i] <-  covs[,1:p]%*%gamma[i,1:p] 
      C[i,,] <-  diag(rep(tau2,Nsites))
      y[,i] ~ dmnorm(mu[,i], cov =  C[i,,] )
      fitted[,i]~ dmnorm(mu[,i], cov = C[i,,])
    }
    
    # Priors
    for(k in 1:p){
      beta[k] ~ dnorm(0, sd = 10)
      for(j in 1:Ngroups){
        gamma[j,k] ~ dnorm(beta[k], var = psi2)
      }
    }
    
    tau2 ~ dinvgamma(2,1)
    psi2 ~ dinvgamma(2,1)
    
  })

spatialCodeinits <-list("beta" = rnorm(ncol(covs)), "tau2" = rinvgamma(1,2,1))


model <- nimbleModel(lod_allcampaigns4, 
                     data = list(y = final_predictors_vocs[,4:6] %>% as.matrix(), 
                                 covs = covs %>% as.matrix()),
                     constants = list(Nsites = nrow(final_predictors_vocs), 
                                      Ngroups = 3,
                                      p = ncol(covs)), 
                     dimensions = list(mu = c(nrow(final_predictors_vocs),3),
                                       C = c(3,nrow(final_predictors_vocs), 
                                             nrow(final_predictors_vocs)),
                                       fitted = c(nrow(final_predictors_vocs),3)),
                     inits = spatialCodeinits)

Cmodel <- compileNimble(model)

start_time <- Sys.time()

mcmcConf <- configureMCMC(Cmodel,
                          monitors = c("gamma","beta", "tau2",  "fitted" ))


modelMCMC <- buildMCMC(mcmcConf,
                       enableWAIC = TRUE)

CmodelMCMC <- compileNimble(modelMCMC, project = model)

mcmc.out4 <- runMCMC(CmodelMCMC, niter = 30000,
                     nburnin = 6000, thin = 14, nchains = 2,
                     summary = TRUE, inits = spatialCodeinits, WAIC = TRUE)
end_time <- Sys.time()
end_time - start_time

# Load mcmc chain if you didn't load run the previous part -----

mcmc.out4 <- readRDS(file = "E:/20_sept_2021/Projects_data/spatialvocs/hexane_lur_M4.rds")

# Chains ------

coda_sample_1 <- coda::mcmc(mcmc.out4$samples$chain1) 
coda_sample_2 <- coda::mcmc(mcmc.out4$samples$chain2) 
coda_sample <- coda::mcmc.list(coda_sample_1, coda_sample_2)

n <- ncol(coda_sample_1) #numero de parametros

Rhat_all <-  sapply(1:n, function(x) {
  rstan::Rhat(cbind(coda_sample_1[, x], coda_sample_2[, x]))
})

ess_all <- sapply(1:n, function(x) {
  rstan::ess_bulk(cbind(coda_sample_1[, x], coda_sample_2[, x]))
})

ess_all_tail <- sapply(1:n, function(x) {
  rstan::ess_tail(cbind(coda_sample_1[, x], coda_sample_2[, x]))
})

max(Rhat_all)
min(ess_all)
min(ess_all_tail)

par(mfrow=c(2,2))

plot(mcmc.out4$samples$chain1[ , 'tau2'], type = 'l', xlab = 'iteration', 
     ylab = "tau2", bty="n")
lines(mcmc.out4$samples$chain2[, 'tau2'], col = "red")


par(mfrow=c(3,2))

plot(mcmc.out4$samples$chain1[ , 'beta[1]'], type = 'l', xlab = 'iteration', 
     ylab = "beta[1]", bty="n")
lines(mcmc.out4$samples$chain2[ , 'beta[1]'], col = "red")

plot(mcmc.out4$samples$chain1[ , 'beta[2]'], type = 'l', xlab = 'iteration', 
     ylab = "beta[2]", bty="n")
lines(mcmc.out4$samples$chain2[ , 'beta[2]'], col = "red")

plot(mcmc.out4$samples$chain1[ , 'beta[3]'], type = 'l', xlab = 'iteration', 
     ylab = "beta[3]", bty="n")
lines(mcmc.out4$samples$chain2[ , 'beta[3]'], col = "red")

plot(mcmc.out4$samples$chain1[ , 'beta[4]'], type = 'l', xlab = 'iteration', 
     ylab = "beta[4]", bty="n")
lines(mcmc.out4$samples$chain2[ , 'beta[4]'], col = "red")

# Coefficient Estimates -----

coeff_estimates <- mcmc.out4$summary$all.chains[grepl("gamma" , 
                                                      rownames(mcmc.out4$summary$all.chains)),
                                                c("Mean",  "95%CI_low",   "95%CI_upp")] %>% 
  round(2) %>% 
  as.data.frame()

coeff_estimates$CI <- paste0("(",coeff_estimates[,"95%CI_low"], ",", 
                             coeff_estimates[,"95%CI_upp"], ")")

coeff_estimates %>% select("Mean", "CI")

# Plot fitted values -----

summary_fitted_dec <- mcmc.out4$summary$all.chains[grepl("fitted.*, 1\\]" , rownames(mcmc.out4$summary$all.chains)),] %>% as.data.frame()
summary_fitted_apr <- mcmc.out4$summary$all.chains[grepl("fitted.*, 2\\]" , rownames(mcmc.out4$summary$all.chains)),] %>% as.data.frame()
summary_fitted_aug <- mcmc.out4$summary$all.chains[grepl("fitted.*, 3\\]" , rownames(mcmc.out4$summary$all.chains)),] %>% as.data.frame()

true_dec <- as.vector(final_predictors_vocs[,4])
true_apr <- as.vector(final_predictors_vocs[,5])
true_aug <- as.vector(final_predictors_vocs[,6])

names(true_dec)[1] <- "real"
names(true_apr)[1] <- "real"
names(true_aug)[1] <- "real"


fitted_real_dec <- cbind(summary_fitted_dec, true_dec)
fitted_real_apr <- cbind(summary_fitted_apr, true_apr)
fitted_real_aug <- cbind(summary_fitted_aug, true_aug)


plot_dec <- ggplot(data = fitted_real_dec) +
  geom_point(aes(x = true_dec, y = Mean)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.y.right = NULL) + xlim(c(0.5,3.5)) + ylim(c(0.5,3.5)) + 
  geom_abline(slope=1, intercept = 0) + 
  ylab(expression("Predicted log hexane "~"["*mu*g/m^3*"]")) + 
  xlab(expression("Observed log hexane "~"["*mu*g/m^3*"]")) +   
  ggtitle("December 2005")


plot_april <- ggplot(data = fitted_real_apr) +
  geom_point(aes(x = true_apr, y = Mean)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.y.right = NULL)  + xlim(c(0.5,3.5)) + ylim(c(0.5,3.5)) + 
  geom_abline(slope=1, intercept = 0) +
  ylab(expression("Predicted log hexane "~"["*mu*g/m^3*"]")) + 
  xlab(expression("Observed log hexane "~"["*mu*g/m^3*"]")) +   
  ggtitle("April 2006")


plot_august <- ggplot(data = fitted_real_aug) +
  geom_point(aes(x = true_aug, y = Mean)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.y.right = NULL) + xlim(c(-1,2)) + ylim(c(-1,2)) + 
  geom_abline(slope=1, intercept = 0) + 
  ylab(expression("Predicted log hexane "~"["*mu*g/m^3*"]")) + 
  xlab(expression("Observed log hexane "~"["*mu*g/m^3*"]")) +   
  ggtitle("August 2006")


grid.arrange(plot_dec, plot_april, plot_august, nrow = 1)


# Residual analysis ----

res_dec <- true_dec - summary_fitted_dec$Mean
res_apr <- true_apr - summary_fitted_apr$Mean
res_aug <- true_aug - summary_fitted_aug$Mean

par(mfrow = c(1,3))
plot( summary_fitted_dec$Mean,res_dec, ylab = "Residual December", xlab = "fitted")
abline(h = 0, col = 'red')

plot( summary_fitted_apr$Mean,res_apr, ylab = "Residual April", xlab = "fitted")
abline(h = 0, col = 'red')

plot(summary_fitted_aug$Mean, res_aug,  ylab = "Residual August", xlab = "fitted")
abline(h = 0, col = 'red')

# Plot grid ------


# Load grid data -----

predictors_grid <- read.csv("data/raw/voc_centroids_predictors.csv")

predictors_grid <- predictors_grid %>% select("id", "Easting", "Northing",  
                                              "Residential_1000m",
                                              "Residential_100m",
                                              "Residential_500m",
                                              "Roads_100m",
                                              "Roads_200m",
                                              "Building_1000m",
                                              "Pop_500m",
                                              "Commercial_1000m",
                                              "Resource.and.Industrial_1000m",
                                              "Government.and.Institutional_200m",
                                              "Government.and.Institutional_1000m",
                                              "Waterbody_1000m",
                                              "Av_NOx_50m",
                                              "Av_NOx_500m",
                                              "Av_NOx_1000m") %>% 
  mutate(X_km = Easting/1000, Y_km = Northing/1000)

dist_grid <- dist(predictors_grid[,c("X_km", "X_km")]) %>% as.matrix()

predictors_grid[is.na(predictors_grid)] <- 0


# Compute y_grid -----

# Obtain the unobserved measurements -----

Ngroups <- 3
N <- nrow(predictors_grid)

coda_sample_1 <- mcmc(mcmc.out4$samples$chain1) 
coda_sample_2 <- mcmc(mcmc.out4$samples$chain2) 
coda_sample <- mcmc.list(coda_sample_1, coda_sample_2)
samples_all_chains <- combine.mcmc(coda_sample)

niter <- nrow(samples_all_chains)

y_c_list <- list(matrix(nrow =  N, ncol = Ngroups))
y_c <- matrix(nrow = N, ncol = Ngroups)

mu <- vector(length = N)


tau_mcmc <- samples_all_chains[,"tau2"]


#gamma_mcmc <- array(NA, dim = c( Ngroups, niter, ncol(covs)))

func_gamma <- function(x){
  gamma_chain <- samples_all_chains %>% as.data.frame() %>% 
    select(starts_with(paste0("gamma[",x,",")))
  return(gamma_chain)
}


gamma_mcmc <- lapply(1:Ngroups, func_gamma) 


unlistaux <- function(lista, iter, param){
  vec_lista <- c()
  for(i in 1:param){
    for(j in 1:iter){
      for(k in 1:length(lista)){
        vec_lista <- c(vec_lista, lista[[k]][j, i])
      }
    }
  }
  return(vec_lista)
}


gamma_mcmc_array <- array(data = unlistaux(gamma_mcmc, niter, ncol(covs_no_scale)),
                          dim = c(Ngroups, niter,  ncol(covs_no_scale)))

#alpha_mcmc <- mcmc.out2$samples$chain1 %>% as.data.frame() %>% 
#  select(starts_with("alpha"))

covs_grid <- data.frame(rep(1, N), 
                        predictors_grid[,c("X_km", "Y_km",
                                           "Residential_1000m",
                                           "Residential_100m",
                                           "Residential_500m",
                                           "Roads_100m",
                                           "Roads_200m",
                                           "Building_1000m",
                                           "Pop_500m",
                                           "Commercial_1000m",
                                           "Resource.and.Industrial_1000m",
                                           "Government.and.Institutional_200m",
                                           "Government.and.Institutional_1000m",
                                           "Waterbody_1000m",
                                           "Av_NOx_50m",
                                           "Av_NOx_500m",
                                           "Av_NOx_1000m")]) %>% as.matrix()



covs_grid[,-1]  <-  scale(covs_grid[,-1] , center = colMeans(covs_no_scale[,-1]), scale = apply(covs_no_scale[,-1], 2, sd))

Sigma_sample <- array(dim = c(N,N,3))


myfnct <- function(j){
  for(i in 1:Ngroups){
    mu <-  covs_grid%*%gamma_mcmc_array[i,j,]
    y_c[,i] <- rmvnorm(1, mu, diag(rep(tau_mcmc[j], N)))
  }
  return(y_c)
}

j <- 1:niter
y_c_list <- mapply(myfnct, j, SIMPLIFY = F)

saveRDS(y_c_list, "data/processed/posteriors_maps/posterior_hexane_map.rds")

#y_c_list <- readRDS("data/processed/posteriors_maps/posterior_hexane_map.rds")
# 
# for(k in 1:nrow(beta_mcmc)){
#   for(i in 1:Ngroups){
#     # Sigma_sample[,,i] <-  sigma_mcmc[k,i]*exp((-1)*dist_grid/phi_mcmc[k,i])
#     # diag(Sigma_sample[,,i]) <- tau_mcmc[k] +  sigma_mcmc[k,i]
#     mu <- covs_grid%*%gamma_mcmc_array[i,k,]
#     y_c[,i] <- rmvnorm(1, mu, diag(rep(tau_mcmc[k], N)))
#   }
# 
#   print(k)
#   y_c_list[[k]] <- y_c
# 
# }


y_c_mean <- apply(simplify2array(y_c_list), 1:2, mean) %>% 
  as.data.frame()


y_c_025 <- apply(simplify2array(y_c_list), 1:2, function(x) quantile(x,probs = 0.025, na.rm = T)) %>% 
  as.data.frame() 


y_c_975 <- apply(simplify2array(y_c_list), 1:2, function(x) quantile(x,probs = 0.975, na.rm = T)) %>% 
  as.data.frame() 

y_c_sd <- apply(simplify2array(y_c_list), 1:2, function(x) sd(x,na.rm = TRUE)) %>% 
  as.data.frame() 


# Plot maps grid  -----

grid$means_dec <- y_c_mean$V1
grid$sd_dec <- y_c_sd$V1

grid$means_apr <- y_c_mean$V2
grid$sd_apr <- y_c_sd$V2

grid$means_aug <- y_c_mean$V3
grid$sd_aug <- y_c_sd$V3

grid$means_campaign <- rowMeans(y_c_mean)
grid$sd_campaigns <- rowMeans(y_c_sd)

grid_mean <- grid %>% 
  gather(key = "Campaign", 
         value = "Mean_hexane", -c(id, left, top, right, bottom, geometry, 
                                   sd_apr, sd_aug, sd_dec, sd_campaigns ))

grid_sd <- grid %>% 
  gather(key = "Campaign", 
         value = "Sd_hexane", -c(id, left, top, right, bottom, geometry, 
                                 means_apr, means_aug, means_dec, means_campaign))

grid_mean$Campaign <- factor(grid_mean$Campaign, 
                             levels = c("means_dec", "means_apr", "means_aug", 
                                        "means_campaign"), 
                             labels = c("December 2005", "April 2006",
                                        "August 2006", "Mean"))


grid_sd$Campaign <- factor(grid_sd$Campaign, 
                           levels = c("sd_dec", "sd_apr", "sd_aug", 
                                      "sd_campaigns"), 
                           labels = c("December 2005", "April 2006",
                                      "August 2006", "Mean"))


# st_write(grid_mean[grid_mean$Campaign == "April 2006", ], "output/hexane_mean_april.shp")
# st_write(grid_mean[grid_mean$Campaign == "December 2005", ], "output/hexane_mean_december.shp")
# st_write(grid_mean[grid_mean$Campaign == "August 2006", ], "output/hexane_mean_august.shp")



plot_mean <- ggplot(grid_mean) + 
  geom_sf(aes(fill = Mean_hexane), colour = NA) +  
  ditch_the_axes + 
  labs(fill = expression ("Log Hexane")) + 
  geom_point(data = hexane_data, 
             aes(x = X_km*1000, y = Y_km*1000), 
             color = "red", size = 1, alpha = 1/3) + 
  scale_color_viridis(discrete = FALSE, option = "D") +
  scale_fill_viridis(discrete = FALSE) +
  facet_wrap(~ Campaign, nrow = 1)

#ggsave("output/figures/mean_pred_hexane.png")

plot_sd <-  ggplot(grid_sd) + 
  geom_sf(aes(fill = Sd_hexane), colour = NA) +  
  ditch_the_axes + 
  labs(fill = expression ("Log Hexane")) + 
  geom_point(data = hexane_data, 
             aes(x = X_km*1000, y = Y_km*1000), 
             color = "red", size = 1, alpha = 1/3) + 
  scale_color_viridis(discrete = FALSE, option = "D") +
  scale_fill_viridis(discrete = FALSE) +
  facet_wrap(~ Campaign, nrow=1)

#ggsave("output/figures/sd_pred_hexane.png")


