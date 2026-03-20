##install.packages("pwrss")
library(pwrss)

pwrss.z.logreg(p0 = 0.5, odds.ratio = 1/0.7, r2.other.x = 0,
               alpha = 0.05, power = 0.80,
               dist = "bernoulli")


pwrss.f.reg(r2 = 0.3, k = 5, power = 0.80, alpha = 0.05)

pwrss.t.reg(beta1 = .02, k = 5, r2 = .3, sdx = sqrt(.5*(1-.5)), sdy = .05,
            power = .80, alpha = 0.05, alternative = "not equal")


p <- 0.50
pwrss.t.reg(beta1 = 0.2, k = 3, r2 = 0.5, sdx = sqrt(p*(1-p)),
            power = .80, alpha = 0.05, alternative = "not equal")


pwrss.z.prop(p = 0.45, p0 = 0.50, arcsin.trans=T,
             alpha = 0.05, power = 0.80,
             alternative = "not equal")

pwr.2p.test(h = ES.h(p1 = 0.5, p2 = 0.45), sig.level = 0.05, power = .80)


library(pwr)
## Calculate Cohen's d
m1 = 2
m2 = 3.3
sd1 = .63
sd2 = .45
d <- abs(m1-m2)/((sd1^2+sd2^2)/2)^.5
d

pwr.t.test(d=.2, sig.level = 0.05, power = 0.8)

nominal.alpha <- 0.05
delta.AUC <- 0.51
target.power.AUC <- 0.8
n.AUC <- 
  power.anova.test(groups = 2, power = target.power.AUC, 
                   between.var = var(c(0, delta.AUC)),
                   within.var = 1, sig.level = nominal.alpha)$n


#######
## Assuming a logistic regression (binary outcome)
library(pwrss)
library(pwr)

## Adsorption - PK dynamics AUC change at least 167% Bagchus et al 2019
## Assumed within-subject coefficient of variation of 51% and 26.5%
## power of 0.8 and 0.05 significance
power.anova.test(groups = 2, power = 0.8, 
                 between.var = var(c(.265, .51)),
                 within.var = 1, sig.level = 0.05)$n

#d = 0.41 medium cohens d, f = 0.2, f2 = 0.04 
pwrss.f.ancova(f2=.04, n.levels = 4, n.covariates = 5, alpha = 0.05,
               power = .8)
## sample size of 199 (so 100 per group)

#d = 0.2 small cohens d, f = 0.1, f2 = 0.01
pwrss.f.ancova(f2=.1, n.levels = 4, n.covariates = 5, alpha = 0.05,
               power = .8)

## sample size of 787 (so 400 per group)


## Metabolism


##side effect occurrence - from Kabatende et al 17.7%, OR 1.39
##drug efficacy - from clark et al. drops from ~100% to 30-40%


# P0 is the base probability under null hypothesis
# odds.ratio = exp(beta1), with beta1 being the regression coeff for predictor X
# r2.other.x is the proportion of the variance of X explained by other predictors 
pwrss.z.logreg(p0 = 0.177, odds.ratio = 1.39, r2.other.x = 0.05,
               alpha = 0.05, power = 0.80,
               dist = "bernoulli")

## This gives a value of 1900
## A probability of the adverse effect at baseline happening with at least 17.7% prevalence
## Lowest CI, without treatment (taken from Kabatende et al. 2022)
## Odds ratio between groups of 1.39 (taken from Kabatende et al. 2022)
## case-control, so r2.other.x = 0 (or very low if we are not fully random i.e. 0.05)
## 95% confidence, 80% power



pwrss.z.2props(p2 = 0.25, p1 = 0.177,
             alpha = 0.05, power = 0.80,
             alternative = "not equal",
             arcsin.trans = FALSE)







###### FINAL SAMPLE SIZE CALCULATIONS FOR REVISED SUBMISSION
## 140 * 6
## 140 * 3
library(pwrss)
library(pwr)

### Change in prevalence
#Power calculation for two proportions 
# Calculate a difference between 2 proportions
pwr.2p.test(h=.4, sig.level = 0.05, power = 0.8)

## Just under 100 per group/arm for moderate cohen's d

### Change in side effects
# P0 is the base probability under null hypothesis
# odds.ratio = exp(beta1), with beta1 being the regression coeff for predictor X
# r2.other.x is the proportion of the variance of X explained by other predictors 
# Baseline of side effects 0.177, odds ratio 1.8, other covariates explain 0.05 of the variance
pwrss.z.logreg(p0 = 0.177, odds.ratio = 1.80, r2.other.x = 0.05,
               alpha = 0.05, power = 0.80,
               dist = "bernoulli")

## 561 in total for two arms (correcting for species)

### PK dynamics
## Adsorption - PK dynamics AUC change at least 167% Bagchus et al 2019
#d = 0.41 medium cohens d, f = 0.2, f2 = 0.04 
pwrss.f.ancova(f2=.04, n.levels = 2, n.covariates = 3, alpha = 0.05,
               power = .8)

## 199 in total for two arms, correcting for species


#### Evaluate micronutrient deficiencies
#### Will need to change the probabilities in the bernouilli 
#### (i.e. the balance between the group deficient vs normal)

### Change in prevalence
## deficient group can be 7.5% (or more) of the total sample
pwr.2p2n.test(h=.4, n1=53, sig.level = 0.05, power = 0.8)


### Change in side effects
# P0 is the base probability under null hypothesis
# odds.ratio = exp(beta1), with beta1 being the regression coeff for predictor X
# r2.other.x is the proportion of the variance of X explained by other predictors 
# Baseline of side effects 0.177, odds ratio 1.8, other covariates explain 0.05 of the variance
pwrss.z.logreg(p0 = 0.177, odds.ratio = 1.80, r2.other.x = 0.05,
               alpha = 0.05, power = 0.80,
               distribution = list(dist = "bernoulli", prob = 0.2))

## PK dynamics
## 14% gives us under the 420 mark that we need.
p = .14
pwrss.t.reg(beta1 = 0.4,
            k = 1, sdx = sqrt(p*(1-p)),
            alpha = 0.05, power = 0.80)
