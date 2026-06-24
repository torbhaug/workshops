# Setup
rm(list=ls())

# if not installed:
# install.packages("EMC2")
library(EMC2)


# Here we set a seed that's appropriate with parallel calculations
# For exact replication of manuscript results purposes.
RNGkind("L'Ecuyer-CMRG")
set.seed(123456)

# Load data
dat <- forstmann

# Examine data
head(dat)

# Plot density of RT distributions
plot_density(dat, factors = c("S", "E"), layout = c(2,3))

# Extract data for first subject for single-subject DDM
dat_single <- dat[dat$subjects == levels(dat$subjects)[1],]
dat_single$subjects <- factor(dat_single$subjects)


### Single Subject DDM ###
# Phase 1: Design Specification
# Specify DDM model with threshold varying by emphasis condition
### THIS REPLICATES CODE BLOCK 1 Single Subject DDM ###
design_DDMaE <- design(data = dat_single,model=DDM,constants=c(s=log(1)),
                       formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1))

### THIS REPLICATES CODE BLOCK 2 Single Subject DDM ###
mapped_pars(design_DDMaE)

# Phase 2: Prior specification
prior_mean=c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
           t0=log(.2),Z=qnorm(.5),sv=log(.5))


# Set prior standard deviation
prior_sd <- c(v_Sleft=1,v_Sright=1,a=.3,a_Eneutral=.3,a_Eaccuracy=.3,
           t0=.4,Z=1,sv=.4,SZ=1)

# Construct prior
prior_DDMaE <- prior(design_DDMaE,pmean = prior_mean, psd = prior_sd,
                     type = "single")


# Visualize prior
### THIS REPLICATES FIGURE 3 Single Subject Design ###
plot(prior_DDMaE, layout = c(2,4))

# Phase 3: Model Estimation
# Construct emc object
DDMaE <- make_emc(dat_single, design_DDMaE, prior = prior_DDMaE, type = "single")

# Fit DDM to first subject
### THIS REPLICATES CODE BLOCK 4 Single Subject DDM ###
DDMaE <- fit(DDMaE)

# Estimation checking
### THIS REPLICATES FIGURE 4 Single Subject DDM ###
check(DDMaE)

# Detailed parameter statistics
### THIS REPLICATES CODE BLOCK 6 Single Subject DDM ###
summary(DDMaE)

# Phase 4: Model Assessment
# Compare posterior to prior using plot
### THIS REPLICATES FIGURE 5 Single Subject DDM ###
plot_pars(DDMaE, layout = c(2,4))



# Generate posterior predictives
pp_DDMaE <- predict(DDMaE, n_cores = 12)

# Function to use to aggregate for levels of response correctness
acc_fun <- function(data) return(data$S == data$R)

# Plot model fit with CDFs
### THIS REPLICATES FIGURE 6 Single Subject DDM ###
plot_cdf(DDMaE, pp_DDMaE, factors = "E", functions = c(correct = acc_fun),
         layout = c(1,3), defective_factor = "correct",
         legendpos = c("bottomright", "right"))

### Hierarhical DDM ###
# Phase 1: Design Specification
# Create stimulus contrast matrix and create
# Design with both drift and threshold varying with emphasis
### THIS REPLICATES CODE BLOCK 1 Hierarchical DDM ###
Smat <- matrix(c(-1,1), nrow = 2,dimnames=list(NULL,"dif"))
design_DDMavE <- design(data = dat, model=DDM,
                        contrasts = list(S = Smat),
                        formula=list(v~S*E,a~E, t0~1, s~1, Z~1, sv~1),
                        constants=c(s=log(1), v_Eneutral = 0, v_Eaccuracy = 0))

# Inspect parameter mappings
### THIS REPLICATES CODE BLOCK 2 Hierarchical DDM ###
mapped_pars(design_DDMavE)

# Phase 2: Prior Specification
### THIS REPLICATES CODE BLOCK 3 Hierarchical DDM ###
prior_DDMavE <- prior(design_DDMavE, update = prior_DDMaE,
                    mu_mean=c(v=0,v_Sdif=2,'v_Sdif:Eneutral'=0.1,
                              'v_Sdif:Eaccuracy'=0.2),
                    mu_sd=c(v=1,v_Sdif=1,'v_Sdif:Eneutral'=1,
                            'v_Sdif:Eaccuracy'=1), type = "standard")


# Phase 3: Model Estimation
DDMavE <- make_emc(dat, design = design_DDMavE, prior = prior_DDMavE)
DDMavE <- fit(DDMavE, cores_per_chain = 4)

# Phase 4: Model Assessment
# Visualize design with evidence accumulation process
### THIS REPLICATES FIGURE 7 + CODE BLOCK 4 Hierarchical DDM ###
plot_design(DDMavE, factors = list(v = c("S", "E"), a = "E"),
            plot_factor = "S")

# Generate posterior predictives (these are used later in comparison with the LBA)
pp_DDMavE <- predict(DDMavE, n_cores = 12)


# Individual differences visualization
### THIS REPLICATES FIGURE 8 + CODE BLOCK 5 Hierarchical DDM ###
plot_pars(DDMavE, layout = c(1,4),
          use_par = c("a_Eneutral", "a_Eaccuracy",
                      "v_Sdif:Eneutral","v_Sdif:Eaccuracy"),
          all_subjects = TRUE)

# Parameter difference hypothesis tests
### THIS REPLICATES CODE BLOCK 6 Hierarchical DDM ###
vdiff <- function(p)diff(p[c("v_Sdif:Eneutral","v_Sdif:Eaccuracy")])
adiff <- function(p)diff(p[c("a_Eneutral","a_Eaccuracy")])
# Note this relies on density estimation, which under the hood varies between
# OS, so results might vary slightly between OS.
hypothesis(DDMavE,fun=vdiff)
hypothesis(DDMavE,fun=adiff)

### Hierarchical LBA ###

# Phase 1: Design Specification
# Match function for race models
### THIS REPLICATES CODE BLOCK 1 Hierarchical LBA ###
matchfun=function(data)data$S==data$lR

# Specifying design LBA
### THIS REPLICATES CODE BLOCK 2 Hierarchical LBA ###
E2 <- function(data) factor(data$E!="speed",labels=c("speed","nonspeed"))
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
design_LBABvE <- design(data = dat,model=LBA,matchfun=matchfun,
                       formula=list(v~lM*E2,sv~lM,B~E2+lR,A~1,t0~1),
                       contrasts=list(lM = ADmat),constants=c(sv=log(1)),
                       functions=list(E2=E2))

# Parameter mappings
### THIS REPLICATES CODE BLOCK 3 Hierarchical LBA ###
mapped_pars(design_LBABvE)

# Phase 2: Prior Specification
mu_mean <- c(v=1, v_E2nonspeed = -.2, v_lMd=1, "v_lMd:E2nonspeed"=.2,
             sv_lMd=log(1),B=log(1), B_E2nonspeed = log(1.5), B_lRright=0,
             A=log(0.25), t0=log(.2))
mu_sd <- c(v=1, v_lMd=0.5, "v_lMd:E2nonspeed"=0.5,
           sv_lMd=.5,B=0.3, B_E2nonspeed=0.3,B_lRright=0.3, A=0.4, t0=.5)

prior_LBABvE <- prior(design_LBABvE, type = 'standard',mu_mean=mu_mean,
                      mu_sd=mu_sd)

# Phase 3: Model Estimation
LBABvE <- make_emc(dat,design_LBABvE, prior=prior_LBABvE)
LBABvE <- fit(LBABvE, cores_per_chain = 4)


# Phase 4: Model Assessment
# Generate posterior predictives
pp_LBABvE <- predict(LBABvE, n_cores = 12)

# Compare LBA and DDM models using CDF
### THIS REPLICATES FIGURE 9 Hierarchical LBA ###
plot_cdf(LBABvE, list(LBA = pp_LBABvE, DDM = pp_DDMavE), factors = "E",
         functions = c(correct = acc_fun), layout = c(1,3),
         defective_factor ="correct", legendpos = c("right", "topleft"))

# Functions to calculate RT and accuracy differences
drt_ns <- function(data){
  all <- tapply(data$rt,data$E,mean)
  out <- all['neutral'] - all['speed']
  names(out) <- "Neutral - Speed"
  return(out)
}
drt_as <- function(data){
  all <- tapply(data$rt,data$E,mean)
  out <- all['accuracy'] - all['speed']
  names(out) <- "Accuracy - Speed"
  return(out)
}
derr_ns <- function(data){
  data$correct <- data$S == data$R
  all <- tapply(data$correct,data$E,mean)*100
  out <- all['neutral'] - all['speed']
  names(out) <- "Neutral - Speed"
  return(out)
}
derr_as <- function(data){
  data$correct <- data$S == data$R
  all <- tapply(data$correct,data$E,mean)*100
  out <- all['accuracy'] - all['speed']
  names(out) <- "Accuracy - Speed"
  return(out)
}

# Plot statistic comparisons between DDM and LBA
### THIS REPLICATES FIGURE 10 Hierarchical LBA ###
par(mfrow = c(2,2))
plot_stat(LBABvE, list(LBA = pp_LBABvE, DDM = pp_DDMavE), stat_fun = drt_ns,
          legendpos = c('topleft', 'topright'), main = "Neutral - Speed",
          layout=NULL,xlim=c(.08,.16),lwd=2,adjust=1,xlab = "RT (sec)",
          posterior_args = list(lwd = 2), prior_args = list(lwd = 2))
plot_stat(LBABvE, list(LBA = pp_LBABvE, DDM = pp_DDMavE), stat_fun = drt_as,
          legendpos = c('topleft', 'topright'), main = "Accuracy - Speed",
          layout=NULL,xlim=c(.08,.16),lwd=2,adjust=1,xlab = "RT (sec)",
          posterior_args = list(lwd = 2), prior_args = list(lwd = 2))
plot_stat(LBABvE, list(LBA = pp_LBABvE, DDM = pp_DDMavE), stat_fun = derr_ns,
          legendpos = c('topleft', 'topright'), main = "Neutral - Speed",
          layout=NULL,xlim=c(7,16),lwd=2,adjust=2, xlab = "Accuracy (%)",
          posterior_args = list(lwd = 2), prior_args = list(lwd = 2))
plot_stat(LBABvE, list(LBA = pp_LBABvE, DDM = pp_DDMavE), stat_fun = derr_as,
          legendpos = c('topleft', 'topright'), main = "Accuracy - Speed",
          layout=NULL,xlim=c(7,16),lwd=2,adjust=2, xlab="Accuracy (%)",
          posterior_args = list(lwd = 2), prior_args = list(lwd = 2))

# Model comparison
### THIS REPLICATES CODE BLOCK 4 Hierarchical LBA ###
compare(list(LBA = LBABvE, DDM = DDMavE))

# Hypothesis tests for LBA parameters
### THIS REPLICATES CODE BLOCK 4 Hierarchical LBA ###
b_diff <- function(p)diff(p[c("B","B_E2nonspeed")])
qual_diff <- function(p)diff(p[c("v_lMd","v_lMd:E2nonspeed")])
quant_diff <- function(p)diff(p[c("v","v_E2nonspeed")])
par(mfrow=c(1,3))
hypothesis(LBABvE,fun=b_diff)
hypothesis(LBABvE,fun=qual_diff)
hypothesis(LBABvE,fun=quant_diff)

# Parameter credible intervals
### THIS REPLICATES CODE BLOCK 5 Hierarchical LBA ###
credint(LBABvE, use_par = c("v_E2nonspeed", "v_lMd:E2nonspeed", "B_E2nonspeed"))

