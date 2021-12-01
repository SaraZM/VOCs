library(bayesplot)
library(coda)
library(grid)
library(gridExtra)
library(mvtnorm)
library(nimble)
library(sf)
library(tidyverse)

# Functions -----
set.seed(123)

ditch_the_axes <- theme(
  axis.text = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_blank())

# Load data -----

predictors <- read.csv("data/all_campaigns_predictors.csv") 
gas_vocs <- read.csv("data/gasoline_vocs.csv")

# Data processing -----

decane_data <- gas_vocs %>% filter(type == "nDecane") %>% 
  mutate(log_concentration = log(concentration)) %>% 
  select(c("ID",  "X_km", "Y_km", "campaign", "log_concentration")) %>%  
  spread(key = "campaign", value = "log_concentration") %>% 
  select( "ID", "X_km", "Y_km", "Dec 2005", "Apr 2006", "Aug 2006") %>% 
  na.omit()


voc_predictors <- merge(decane_data,
                           predictors, by = "ID")

# If the predictor is missing it means it is a zero
# Note: there are no missing values in the voc data

voc_predictors[is.na(voc_predictors)] <- 0

# Decane predictors ----

decane_predictors <- decane_predictors %>% select("ID", "X_km", "Y_km", 
                                                      "Dec 2005","Apr 2006", 
                                                      "Aug 2006", 
                                                      "Av_NOx_50m",
                                                      "Av_NOx_200m",
                                                      "Av_NOx_500m",
                                                      "Av_NOx_1000m",
                                                      "Building_50m",
                                                      "Building_100m",
                                                      "Building_200m",
                                                      "Building_500m",
                                                      "Commercial_1000m",
                                                      "Government.and.Institutional_200m",
                                                      "Government.and.Institutional_500m",
                                                      "Residential_50m",
                                                      "Residential_1000m",
                                                      "Roads_50m",
                                                      "Roads_100m",
                                                      "Roads_200m",
                                                      "Open.Area_100m",
                                                      "Pop_500m",
                                                      "Pop_1000m",
                                                      "Resource.and.Industrial_500m",
                                                      "Resource.and.Industrial_1000m",
                                                      "Waterbody_1000m")

decane_predictors[,c("X_km", "Y_km")] <- scale(decane_predictors[,c("X_km", "Y_km")])

campaigns_indx <- which(colnames(benzene_predictors) %in% c("ID","Dec 2005", "Apr 2006", "Aug 2006"))
decane_predictors[,-campaigns_indx] <- scale(decane_predictors[,-campaigns_indx])

# Calculate distance matrix in km

dist_matrix <- dist(decane_data[c("X_km", "Y_km")]) %>% as.matrix()

# Define variables for program -----

Nsites <- nrow(final_predictors_vocs)

ones <- rep(1,Nsites) %>% as.data.frame()
D1 <- c(0,1,0)
D2 <- c(0,0,1)

covs <- data.frame(ones, decane_predictors[,-campaigns_indx])

D <- data.frame(D1 , D2) %>% as.matrix()


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


spatialCodeinits <-list("alpha" = rnorm(2), "beta" = rnorm(ncol(covs)), 
                        "phi" = rexp(3,0.3), 
                        "tau2" = rinvgamma(1,2,1),
                        "sigma2" = rinvgamma(3,2,1))


model <- nimbleModel(model1, data = list(y = final_predictors_vocs[,4:6] %>% as.matrix(), 
                                                   covs = covs %>% as.matrix(),
                                                   dist_matrix = decane_dist,
                                                   dummy = D),
                     constants = list( Nsites = Nsites, 
                                       Ngroups = 3,
                                       p = ncol(covs)), 
                     dimensions = list(mu = c(Nsites,3),
                                       C = c(3,Nsites, Nsites),
                                       fitted = c(nrow(final_predictors_vocs),3)),
                     inits = spatialCodeinits)

Cmodel <- compileNimble(model)

start_time <- Sys.time()

mcmcConf <- configureMCMC(Cmodel,
                          monitors = c("alpha","beta", 
                                       "tau2", "phi", 
                                       "sigma2", "fitted"))

modelMCMC <- buildMCMC(mcmcConf,
                       enableWAIC = TRUE)

CmodelMCMC <- compileNimble(modelMCMC, project = model)

mcmc.out <- runMCMC(CmodelMCMC, niter = 20000,
                    nburnin = 1000, thin = 14, nchains = 2,
                    summary = TRUE, inits = spatialCodeinits, WAIC = TRUE)

end_time <- Sys.time()
end_time - start_time

# Check convergence 

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
      mu[,i] <- rep(dummy[i,]%*%alpha[1:(Ngroups-1)], Nsites) + covs[,1:p]%*%beta[1:p] 
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

spatialCodeinits <-list("beta" = rnorm(ncol(covs)), 
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

Cmodel <- compileNimble(model)

start_time <- Sys.time()

mcmcConf <- configureMCMC(Cmodel)


modelMCMC <- buildMCMC(mcmcConf,
                       enableWAIC = TRUE)

CmodelMCMC <- compileNimble(modelMCMC, project = model)

mcmc.out2 <- runMCMC(CmodelMCMC, niter = 20000,
                     nburnin = 2000, thin = 14, nchains = 2,
                     summary = TRUE, inits = spatialCodeinits, WAIC = TRUE)
end_time <- Sys.time()
end_time - start_time

# Check convergence 

coda_sample_1 <- coda::mcmc(mcmc.out2$samples$chain1) 
coda_sample_2 <- coda::mcmc(mcmc.out2$samples$chain2) 
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
      mu[,i] <- covs[,1:p]%*%gamma[i,1:p]
      y[,i] ~ dmnorm(mu[,i], cov = C[i,,])
      fitted[,i]~ dmnorm(mu[,i], cov = C[i,,])
    }
    
    # Priors
    for(k in 1:p){
      beta[k] ~ dnorm(0, sd = 10)
      for(j in 1:Ngroups){
        gamma[j,k] ~ dnorm(beta[k], var = psi2)
        
      }
    }
    
    
    for(i in 1:Ngroups){
      sigma2[i] ~ dinvgamma(2,1)
      phi[i] ~ dexp(0.1)
    }
    
    tau2 ~ dinvgamma(2,1)
    psi2 ~ dinvgamma(2,1)
    
  })

spatialCodeinits <-list("beta" = rnorm(ncol(covs)), 
                        "phi" = rexp(3,0.3), 
                        "tau2" = rinvgamma(1,2,1),
                        "sigma2" = rinvgamma(3,2,1))


model <- nimbleModel(model3, data = list(y = final_predictors_vocs[,4:6] %>% as.matrix(), 
                                                    covs = covs %>% as.matrix(),
                                                    dist_matrix = decane_dist),
                     constants = list( Nsites = Nsites, 
                                       Ngroups = 3,
                                       p = ncol(covs)), 
                     dimensions = list(mu = c(Nsites,3),
                                       C = c(3,Nsites, Nsites),
                                       fitted = c(nrow(final_predictors_vocs),3)),
                     inits = spatialCodeinits)

Cmodel <- compileNimble(model)

start_time <- Sys.time()

mcmcConf <- configureMCMC(Cmodel,
                          monitors = c("beta", "gamma", "psi2",
                                       "tau2", "phi", 
                                       "sigma2", "fitted" ))


modelMCMC <- buildMCMC(mcmcConf,
                       enableWAIC = TRUE)

CmodelMCMC <- compileNimble(modelMCMC, project = model)

mcmc.out3 <- runMCMC(CmodelMCMC, niter = 20000,
                     nburnin = 2000, thin = 14, nchains = 2,
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

spatialCodeinits <-list("beta" = rnorm(ncol(covs)),"tau2" = rinvgamma(1,2,1))


model <- nimbleModel(model4, 
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
                          monitors = c("beta", "gamma", 
                                       "tau2", "fitted" ))


modelMCMC <- buildMCMC(mcmcConf,
                       enableWAIC = TRUE)

CmodelMCMC <- compileNimble(modelMCMC, project = model)

mcmc.out4 <- runMCMC(CmodelMCMC, niter = 20000,
                     nburnin = 2000, thin = 5, nchains = 2,
                     summary = TRUE, inits = spatialCodeinits, WAIC = TRUE)
end_time <- Sys.time()
end_time - start_time

# Check convergence

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

# Print WAIC all models 

mcmc.out$WAIC
mcmc.out2$WAIC
mcmc.out3$WAIC
mcmc.out4$WAIC

# Plot chains best model ----

par(mfrow=c(2,1))

plot(mcmc.out$samples$chain1[ , 'alpha[1]'], type = 'l', xlab = 'iteration', 
     ylab = expression(alpha[1]), bty="n")
lines(mcmc.out$samples$chain2[ , 'alpha[1]'], col = "red")

plot(mcmc.out$samples$chain1[ , 'alpha[2]'], type = 'l', xlab = 'iteration', 
     ylab = expression(alpha[2]), bty="n")
lines(mcmc.out$samples$chain2[ , 'alpha[2]'], col = "red")

par(mfrow=c(2,2))

plot(mcmc.out$samples$chain1[ , 'beta[1]'], type = 'l', xlab = 'iteration', 
     ylab = expression(beta[1]), bty="n")
lines(mcmc.out$samples$chain2[ , 'beta[1]'], col = "red")

plot(mcmc.out$samples$chain1[ , 'beta[2]'], type = 'l', xlab = 'iteration', 
     ylab = expression(beta[2]), bty="n")
lines(mcmc.out$samples$chain2[ , 'beta[2]'], col = "red")

plot(mcmc.out$samples$chain1[ , 'beta[3]'], type = 'l', xlab = 'iteration', 
     ylab = expression(beta[3]), bty="n")
lines(mcmc.out$samples$chain2[ , 'beta[3]'], col = "red")

plot(mcmc.out$samples$chain1[ , 'beta[4]'], type = 'l', xlab = 'iteration', 
     ylab = expression(beta[4]), bty="n")
lines(mcmc.out$samples$chain2[ , 'beta[4]'], col = "red")

par(mfrow=c(2,2))

plot(mcmc.out$samples$chain1[ , 'sigma2[1]'], type = 'l', xlab = 'iteration', 
     ylab = expression(sigma[1]^2), bty="n")
lines(mcmc.out$samples$chain2[ , 'sigma2[1]'], col = "red")

plot(mcmc.out$samples$chain1[ , 'sigma2[2]'], type = 'l', xlab = 'iteration', 
     ylab = expression(sigma[2]^2), bty="n")
lines(mcmc.out$samples$chain2[ , 'sigma2[2]'], col = "red")

plot(mcmc.out$samples$chain1[ , 'sigma2[3]'], type = 'l', xlab = 'iteration', 
     ylab = expression(sigma[3]^2), bty="n")
lines(mcmc.out$samples$chain2[ , 'sigma2[3]'], col = "red")

plot(mcmc.out$samples$chain1[ , 'tau2'], type = 'l', xlab = 'iteration', 
     ylab = expression(tau^2), bty="n")
lines(mcmc.out$samples$chain2[ , 'tau2'], col = "red")

par(mfrow=c(2,2))

plot(mcmc.out$samples$chain1[ , 'phi[1]'], type = 'l', xlab = 'iteration', 
     ylab = expression(phi[1]), bty="n")
lines(mcmc.out$samples$chain2[ , 'phi[1]'], col = "red")

plot(mcmc.out$samples$chain1[ , 'phi[2]'], type = 'l', xlab = 'iteration', 
     ylab = expression(phi[2]), bty="n")
lines(mcmc.out$samples$chain2[ , 'phi[2]'], col = "red")

plot(mcmc.out$samples$chain1[ , 'phi[3]'], type = 'l', xlab = 'iteration', 
     ylab = expression(phi[3]), bty="n")
lines(mcmc.out$samples$chain2[ , 'phi[3]'], col = "red")



# Summary of the parameters -----
varnames <- c("alpha[1]", "alpha[2]", 
              "beta[1]", "tau2",
              "phi[1]", "phi[2]","phi[3]",
              "sigma2[1]", "sigma2[2]", "sigma2[3]")

summary_coef <- mcmc.out$summary$all.chains[varnames, ] %>% as.data.frame()


summary_coef$variable <- as.factor(rownames(summary_coef))


xtable(summary_coef[c( "tau2",
                       "phi[1]", "phi[2]","phi[3]",
                       "sigma2[1]", "sigma2[2]", "sigma2[3]"), 
                    c("Mean", "95%CI_low", "95%CI_upp")]) 


# Fitted values -----

summary_fitted_dec <- mcmc.out$summary$all.chains[grepl("fitted.*, 1\\]" , rownames(mcmc.out$summary$all.chains)),] %>% as.data.frame()
summary_fitted_apr <- mcmc.out$summary$all.chains[grepl("fitted.*, 2\\]" , rownames(mcmc.out$summary$all.chains)),] %>% as.data.frame()
summary_fitted_aug <- mcmc.out$summary$all.chains[grepl("fitted.*, 3\\]" , rownames(mcmc.out$summary$all.chains)),] %>% as.data.frame()

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
        axis.line.y.right = NULL) + xlim(c(-1.5,1.6)) + ylim(c(-1.5,1.6)) + 
  geom_abline(slope=1, intercept = 0) + ylab("fitted") + xlab("observed")+ 
  ggtitle("December 2005")


plot_april <- ggplot(data = fitted_real_apr) +
  geom_point(aes(x = true_apr, y = Mean)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.y.right = NULL)  + xlim(c(-0.5,2)) + ylim(c(-0.5,2)) + 
  geom_abline(slope=1, intercept = 0) + ylab("fitted") + xlab("observed")+ 
  ggtitle("April 2006")


plot_august <- ggplot(data = fitted_real_aug) +
  geom_point(aes(x = true_aug, y = Mean)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.y.right = NULL) + xlim(c(-2.2,1.5)) + ylim(c(-2.2,1.5)) + 
  geom_abline(slope=1, intercept = 0) + ylab("fitted") + xlab("observed")+ 
  ggtitle("August 2006")




grid.arrange(plot_dec, plot_april, plot_august, nrow = 1)

# Residual analysis ----

res_dec <- true_dec - summary_fitted_dec$Mean
res_apr <- true_apr - summary_fitted_apr$Mean
res_aug <- true_aug - summary_fitted_aug$Mean

par(mfrow = c(1,3))
plot(summary_fitted_dec$Mean, res_dec, ylab = "Residual December")
abline(h = 0, col = 'red')

plot(summary_fitted_apr$Mean, res_apr, ylab = "Residual April")
abline(h = 0, col = 'red')

plot(summary_fitted_aug$Mean, res_aug, ylab = "Residual August")
abline(h = 0, col = 'red')

scale_res_dec <- scale(res_dec)
scale_res_apr <- scale(res_apr)
scale_res_aug <- scale(res_aug)

par(mfrow = c(1,3))
plot(res_dec, ylab = "Scaled residual December")
abline(h = 0, col = 'red')

plot(res_apr, ylab = "Scaled residual April")
abline(h = 0, col = 'red')

plot(res_aug, ylab = "Scaled residual August")
abline(h = 0, col = 'red')

