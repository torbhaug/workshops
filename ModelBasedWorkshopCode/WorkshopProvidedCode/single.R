rm(list = ls())
library(EMC2)
library(cmdstanr)

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

n_reps <- 10

# Function to fit models and record time
fit_models <- function(rep, n_trials = 50) {
  # Initialize a data frame to store timing results for this repetition
  timing_results <- data.frame(
    repetition = rep,
    emc2_time_ddm5 = NA,
    stan_time_ddm5 = NA,
    emc2_time_ddm6 = NA,
    stan_time_ddm6 = NA,
    emc2_time_ddm7 = NA,
    stan_time_ddm7 = NA
  )
  ESS_results <- data.frame(
    repetition = rep,
    emc2_ESS_ddm5 = NA,
    stan_ESS_ddm5 = NA,
    emc2_ESS_ddm6 = NA,
    stan_ESS_ddm6 = NA,
    emc2_ESS_ddm7 = NA,
    stan_ESS_ddm7 = NA
  )

  # DDM5 Setup and Fitting
  design_DDM5 <- design(factors=list(subjects=1,C=c("HARD", "EASY")),
                        Rlevels = c(0, 1),
                        formula =list(v~C,a~1, t0~1, s~1, Z~1, sv~1),
                        constants=c(s=log(1)),
                        model = DDM)

  pmean5 <- sampled_pars(design_DDM5, doMap = F)
  pmean5[] <- c(1.5, 0.1, log(1), log(.435), qnorm(.5), log(1))
  DDM5dat <- make_data(pmean5, design_DDM5, n_trials = n_trials)

  # # DDM6 Setup and Fitting
  design_DDM6 <- design(factors=list(subjects=1,C=c("HARD", "EASY")),
                        Rlevels = c(0, 1),
                        formula =list(v~C,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
                        constants=c(s=log(1)),
                        model = DDM)
  #
  pmean6 <- sampled_pars(design_DDM6, doMap = F)
  pmean6[] <- c(1.5, 0.1, log(1), log(.435), qnorm(.5), log(1), qnorm(.25))
  DDM6dat <- make_data(pmean6, design_DDM6, n_trials = n_trials)

  # DDM7 Setup and Fitting
  design_DDM7 <- design(factors=list(subjects=1,C=c("HARD", "EASY")),
                        Rlevels = c(0, 1),
                        formula =list(v~C,a~1, t0~1, s~1, Z~1, sv~1, SZ~1, st0 ~ 1),
                        constants=c(s=log(1)),
                        model = DDM)

  pmean7 <- sampled_pars(design_DDM7, doMap = F)
  pmean7[] <- c(1.5, 0.1, log(1), log(.435), qnorm(.5), log(1), qnorm(.25), log(.183))
  DDM7dat <- make_data(pmean7, design_DDM7, n_trials = n_trials)

  # Fit EMC2 Models
  # DDM5
  emc2_time5 <- system.time({
    DDM5_emc2 <- make_emc(DDM5dat, design_DDM5, type = "single", n_chains = 4, compress = FALSE)
    DDM5_emc2 <- fit(DDM5_emc2, iter = 1000)
  })[3]
  emc2_ESS5 <- min(ess_summary(DDM5_emc2))

  # # DDM6
  # emc2_time6 <- system.time({
  #   DDM6_emc2 <- make_emc(DDM6dat, design_DDM6, type = "single", n_chains = 4, compress = FALSE)
  #   DDM6_emc2 <- fit(DDM6_emc2, iter = 1000)
  # })[3]
  # emc2_ESS6 <- min(ess_summary(DDM6_emc2))
  #
  # # DDM7
  # emc2_time7 <- system.time({
  #   DDM7_emc2 <- make_emc(DDM7dat, design_DDM7, type = "single", n_chains = 4, compress = FALSE)
  #   DDM7_emc2 <- fit(DDM7_emc2, iter = 1000)
  # })[3]
  # emc2_ESS7 <- min(ess_summary(DDM7_emc2))

  # Fit Stan Models
  # DDM5
  file5 <- file.path("ddm5.stan")
  stan_time5 <- system.time({
    mod5 <- cmdstan_model(file5,
                          force_recompile = TRUE)
    standata5 <- list(N=nrow(DDM5dat), cnd = ifelse(DDM5dat$C == "HARD", 1, 2),
                      Ncnds = length(unique(DDM5dat$C)), rt = DDM5dat$rt,
                      resp = as.numeric(as.character(DDM5dat$R)), parallel = 0)
    fit5 <- mod5$sample(data = standata5, seed = 123, chains = 4, parallel_chains = 4,
                        iter_warmup = 300, iter_sampling = 1000, refresh = 0)
  })[3]
  stan_ESS5 <- fit5$summary()[,'ess_bulk'] |> min()

  # DDM6
  standata6 <- list(N=nrow(DDM6dat), cnd = ifelse(DDM6dat$C == "HARD", 1, 2),
                    Ncnds = length(unique(DDM6dat$C)), rt = DDM6dat$rt,
                    resp = as.numeric(as.character(DDM6dat$R)), parallel = 0)

  file6 <- file.path("ddm6.stan")
  stan_time6 <- system.time({
    mod6 <- cmdstan_model(file6, force_recompile = TRUE)
    fit6 <- mod6$sample(data = standata6, seed = 123, chains = 4, parallel_chains = 4,
                        iter_warmup = 300, iter_sampling = 1000, refresh = 0)
  })[3]
  stan_ESS6 <- fit6$summary()[,'ess_bulk'] |> min()

  # DDM7
  standata7 <- list(N=nrow(DDM7dat), cnd = ifelse(DDM7dat$C == "HARD", 1, 2),
                    Ncnds = length(unique(DDM7dat$C)), rt = DDM7dat$rt,
                    resp = as.numeric(as.character(DDM7dat$R)), parallel = 0)

  file7 <- file.path("ddm7.stan")
  stan_time7 <- system.time({
    mod7 <- cmdstan_model(file7, force_recompile = TRUE)
    fit7 <- mod7$sample(data = standata7, seed = 123, chains = 4, parallel_chains = 4,
                        iter_warmup = 300, iter_sampling = 1000, refresh = 0)
  })[3]
  stan_ESS7 <- fit7$summary()[,'ess_bulk'] |> min()
  # Store results
  timing_results$emc2_time_ddm5 <- emc2_time5
  timing_results$stan_time_ddm5 <- stan_time5
  timing_results$emc2_time_ddm6 <- emc2_time6
  timing_results$stan_time_ddm6 <- stan_time6
  timing_results$emc2_time_ddm7 <- emc2_time7
  timing_results$stan_time_ddm7 <- stan_time7
  ESS_results$emc2_ESS_ddm5 <- emc2_ESS5
  ESS_results$stan_ESS_ddm5 <- stan_ESS5
  ESS_results$emc2_ESS_ddm6 <- emc2_ESS6
  ESS_results$stan_ESS_ddm6 <- stan_ESS6
  ESS_results$emc2_ESS_ddm7 <- emc2_ESS7
  ESS_results$stan_ESS_ddm7 <- stan_ESS7
  return(list(timing = timing_results, ESS = ESS_results))
}

# Perform repetitions and collect results
for(i in 1:n_reps){
  out[[i]] <- fit_models(i)
  gc()
}

timing_single <- do.call(rbind, lapply(out, \(x) x$timing))
ess_single <- do.call(rbind, lapply(out, \(x) x$ESS))

save(timing_single, file = "results/timing_single.RData")
save(ess_single, file = "results/ess_single")
