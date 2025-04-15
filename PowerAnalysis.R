################
## Power analysis "Food effect on praziquantel efficacy for schistosomiasis"

library(parallel)
library(MASS)

nominal.alpha <- 0.05 # significance threshold
n.sim <- 1000 # number of data sets to simulate

# Study design options

# Total sample size
n <- seq(1000, 2000, 200)

# Prevalence
p <- seq(0.4, 0.8, 0.2)

## Carbohydrate meal as reference
## OR taken from Kabatende et al. 2022
# Matching the quantiles
# quantile(rnorm(1000,1.39,0.12),c(0.025,0.975))
ORprot <- 1.39 #sd 0.12

#quantile(rnorm(1000,1.52,0.19),c(0.025,0.975))
ORfat <- 1.52 #sd 0.19

# correlation of predictors (unrealistic to assume zero correlation)
r <- 0.25

# continuous drivers (age,school)
n.cont <- 2

# binary drivers (Schisto species = Smansoni or Shaematobium, Country = Malawi or Cote),
# both are split 50-50
bin.x.p <- c("Sman" = 0.5, "Cote" = 0.5)

# Total number of drivers
n.x <- n.cont + length(bin.x.p)

# Bonferroni correction (likely too conservative)
alpha <- c(nominal.alpha/n.x)

# Make table of all parameter and design combinations
par.tab <- expand.grid(n = n, p = p, or = ORprot, r = r, alpha = alpha, n.sim = n.sim)

## Run simulations
sim.res <- 
  sapply(1:nrow(par.tab), function(i) {
    
    # Print progress
    print(paste0(round(100*(i-1)/nrow(par.tab)), "% complete"))
    
    # Simulate Xs (the drivers)
    
    # Correlation among Xs
    sigma2 <- diag(n.x)
    sigma2[lower.tri(sigma2)] <- par.tab$r[i]
    sigma2[upper.tri(sigma2)] <- par.tab$r[i]
    
    # Simulate data and analysis n.sim times, returning the mean
    # number of drivers identified at P < alpha
    n.sim.out <-
      mclapply(1:par.tab$n.sim[i], function(k) {
        
        # Simulate Xs before dichotomising (so the binary Xs are derived from latent scales)
        X.raw <- mvrnorm(par.tab$n[i], mu = rep(0, n.x), Sigma = sigma2)
        
        # Dichotomise the binary Xs, but centre on zero by subtracting 0.5
        # and add a column of 1s for the intercept
        X <-
          cbind(1, 
                sapply(1:n.x, function(j) {
                  if(j <= length(bin.x.p)) {
                    X.raw[, j] <- as.integer(X.raw[, j] < qnorm(bin.x.p[j])) - 0.5
                  }
                  X.raw[, j]
                }))
        
        # vector of intercept (log odds) and log odds ratios
        b <- c(qlogis(par.tab$p[i]), log(rep(par.tab$or[i], n.x)))
        
        # Simulate failure to clear from linear predictor X %*% b
        response <- rbinom(par.tab$n[i], 1, plogis(X %*% b))
        
        # Put everything into a data frame
        dat <- data.frame(X[, -1], response = response)
        
        # fit GLM
        fit <- glm(response ~ X[, -1], family = binomial)
        res.tab <- coef(summary(fit))
        
        # How many drivers were detected (P < alpha)?
        sum(res.tab[-1, "Pr(>|z|)"] < par.tab$alpha[i])
        
      })
    
    # take mean number of drivers (Xs) detected across simulated data sets
    mean(unlist(n.sim.out))
    
  })

# Estimate power as the proportion of the n.x drivers with p < alpha
par.tab$prop.drivers <- sim.res/n.x
par.tab
