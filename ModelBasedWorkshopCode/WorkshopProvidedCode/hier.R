rm(list = ls())
library(EMC2)
library(cmdstanr)

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# --- Number of repetitions
n_reps <- 10


# I admit that this storing way is ugly and code-inefficent...
fit_models <- function(rep, n_trials = 50) {
  timing_results <- data.frame(
    repetition = rep,
    emc2_time_ddm5 = NA,
    stan_time_ddm5 = NA,
    emc2_time_ddm6 = NA,
    stan_time_ddm6 = NA,
    emc2_time_ddm7 = NA
  )
  ESS_results <- data.frame(
    repetition = rep,
    emc2_ESS_ddm5_mu = NA,
    emc2_ESS_ddm5_sigma = NA,
    emc2_ESS_ddm5_alpha = NA,
    stan_ESS_ddm5_mu = NA,
    stan_ESS_ddm5_sigma = NA,
    stan_ESS_ddm5_alpha = NA,
    emc2_ESS_ddm6_mu = NA,
    emc2_ESS_ddm6_sigma = NA,
    emc2_ESS_ddm6_alpha = NA,
    stan_ESS_ddm6_mu = NA,
    stan_ESS_ddm6_sigma = NA,
    stan_ESS_ddm6_alpha = NA,
    emc2_ESS_ddm7_mu = NA,
    emc2_ESS_ddm7_sigma = NA,
    emc2_ESS_ddm7_alpha = NA
  )


  # DDM5 --------------------------------------------------------------------


  # 1) Define hierarchical design for 30 subjects, conditions = {HARD, EASY}
  design_DDM5 <- design(
    factors = list(subjects = 1:30, C = c("HARD", "EASY")),
    Rlevels = c(0, 1),
    formula = list(v~C, a~1, t0~1, s~1, Z~1, sv~1),
    constants = c(s = log(1)),  # fix log(s)=0 => s=1
    model = DDM
  )

  # 2) Generate parameter means & random effects, then make data
  pmean5 <- sampled_pars(design_DDM5, doMap = FALSE)
  covs <- diag(c(.2, .075, .075, .02, .02, .02))
  # threshold=1.5, t0=exp(log(0.435))=0.435, etc.
  pmean5[] <- c(1.5, 0.1, log(1), log(.435), qnorm(.5), log(.3))
  DDM5dat <- make_data(
    make_random_effects(design_DDM5, pmean5, covariances = covs),
    design_DDM5,
    n_trials = 50
  )

  library(cmdstanr)

  stan_time5 <- system.time({
    # Compile the model with threading support enabled
    mod <- cmdstan_model(
      "hierarchical_ddm5.stan",
      cpp_options = list(stan_threads = TRUE),
      force_recompile = TRUE
    )

    # Prepare your data list (adjust according to your data)
    data_list <- list(
      N = nrow(DDM5dat),
      Nsubj = length(unique(DDM5dat$subjects)),
      Ncnds = length(unique(DDM5dat$C)),
      subj = as.integer(DDM5dat$subjects),
      cnd = as.integer(DDM5dat$C),
      resp = as.integer(DDM5dat$R) - 1,  # adjust responses to 0/1 if needed
      rt = DDM5dat$rt,
      parallel = 1
    )
    # Fit the model: 4 chains, 300 warm-up and 1000 sampling iterations,
    # and 10 threads per chain (total 40 cores)
    fit <- mod$sample(
      data = data_list,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 300,
      iter_sampling = 1000,
      threads_per_chain = 10,
      seed = 123
    )
  })[3]
  ess_all <- fit$summary()
  mu_idx <- grepl("mu", ess_all$variable)
  sigma_idx <- grepl("sigma", ess_all$variable)
  stan_ESS5_mu <- fit$summary()[mu_idx,'ess_bulk'] |> min()
  stan_ESS5_sigma <- fit$summary()[sigma_idx,'ess_bulk'] |> min()
  stan_ESS5_alpha <- fit$summary()[!(mu_idx & sigma_idx),'ess_bulk'] |> min()

  emc2_time5 <- system.time({
    DDM5_emc2 <- make_emc(DDM5dat, design_DDM5, n_chains = 4, compress = FALSE, type = "diagonal")
    DDM5_emc2 <- fit(DDM5_emc2, cores_per_chain = 10)
  })[3]
  emc2_ESS5_mu <- min(ess_summary(DDM5_emc2, selection = "mu"))
  emc2_ESS5_sigma <- min(ess_summary(DDM5_emc2, selection = "sigma2"))
  emc2_ESS5_alpha <- min(ess_summary(DDM5_emc2, selection = "alpha"))


  # DDM6 --------------------------------------------------------------------


  # 1) Define hierarchical design for 30 subjects, conditions = {HARD, EASY}
  design_DDM6 <- design(
    factors = list(subjects = 1:30, C = c("HARD", "EASY")),
    Rlevels = c(0, 1),
    formula = list(v~C, a~1, t0~1, s~1, Z~1, sv~1, SZ ~ 1),
    constants = c(s = log(1)),  # fix log(s)=0 => s=1
    model = DDM
  )

  # 2) Generate parameter means & random effects, then make data
  pmean6 <- sampled_pars(design_DDM6, doMap = FALSE)
  covs <- diag(c(.2, .075, .075, .02, .02, .02, .02))
  # threshold=1.5, t0=exp(log(0.435))=0.435, etc.
  pmean6[] <- c(1.5, 0.1, log(1), log(.435), qnorm(.5), log(.3), qnorm(.05))
  DDM6dat <- make_data(
    make_random_effects(design_DDM6, pmean6, covariances = covs),
    design_DDM6,
    n_trials = 50
  )

  library(cmdstanr)

  stan_time6 <- system.time({
    # Compile the model with threading support enabled
    mod <- cmdstan_model(
      "hierarchical_ddm6.stan",
      cpp_options = list(stan_threads = TRUE),
      force_recompile = TRUE
    )

    # Prepare your data list (adjust according to your data)
    data_list <- list(
      N = nrow(DDM6dat),
      Nsubj = length(unique(DDM6dat$subjects)),
      Ncnds = length(unique(DDM6dat$C)),
      subj = as.integer(DDM6dat$subjects),
      cnd = as.integer(DDM6dat$C),
      resp = as.integer(DDM6dat$R) - 1,  # adjust responses to 0/1 if needed
      rt = DDM6dat$rt,
      parallel = 1
    )

    init_function <- function() {
      list(
        # Group-level hyperparameters:
        mu_a_tr      = 0,
        sigma_a_tr   = 0.3,
        mu_zr_m_tr   = 0,
        sigma_zr_m_tr= 0.3,
        mu_v_m_tr    = rep(0, data_list$Ncnds),
        sigma_v_m_tr = rep(0.3, data_list$Ncnds),
        mu_v_s_tr    = 0,
        sigma_v_s_tr = 0.3,
        # Choose mu_t0_m_tr so that exp(mu_t0_m_tr) is safely below observed RTs.
        mu_t0_m_tr   = -1.5,
        sigma_t0_m_tr= 0.3,
        mu_zr_sd_tr  = 0,
        sigma_zr_sd_tr= 0.3,

        # Individual-level (non-centered) parameters:
        a_tr_raw     = rep(0, data_list$Nsubj),
        zr_m_tr_raw  = rep(0, data_list$Nsubj),
        # For each subject and condition:
        v_m_tr_raw   = matrix(0, nrow = data_list$Nsubj, ncol = data_list$Ncnds),
        v_s_tr_raw   = rep(0, data_list$Nsubj),
        t0_m_tr_raw  = rep(0, data_list$Nsubj),
        zr_sd_tr_raw = rep(0, data_list$Nsubj)
      )
    }



    # Fit the model: 4 chains, 300 warm-up and 1000 sampling iterations,
    # and 10 threads per chain (total 40 cores)
    fit <- mod$sample(
      data = data_list,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 300,
      iter_sampling = 1000,
      init = init_function,
      threads_per_chain = 10,
      seed = 123
    )
  })[3]
  ess_all <- fit$summary()
  mu_idx <- grepl("mu", ess_all$variable)
  sigma_idx <- grepl("sigma", ess_all$variable)
  stan_ESS6_mu <- fit$summary()[mu_idx,'ess_bulk'] |> min()
  stan_ESS6_sigma <- fit$summary()[sigma_idx,'ess_bulk'] |> min()
  stan_ESS6_alpha <- fit$summary()[!(mu_idx & sigma_idx),'ess_bulk'] |> min()

  emc2_time6 <- system.time({
    DDM6_emc2 <- make_emc(DDM6dat, design_DDM6, n_chains = 4, compress = FALSE, type = "diagonal")
    DDM6_emc2 <- fit(DDM6_emc2, cores_per_chain = 10)
  })[3]
  emc2_ESS6_mu <- min(ess_summary(DDM6_emc2, selection = "mu"))
  emc2_ESS6_sigma <- min(ess_summary(DDM6_emc2, selection = "sigma2"))
  emc2_ESS6_alpha <- min(ess_summary(DDM6_emc2, selection = "alpha"))
  # # DDM7 --------------------------------------------------------------------
  design_DDM7 <- design(
    factors = list(subjects = 1:30, C = c("HARD", "EASY")),
    Rlevels = c(0, 1),
    formula = list(v~C, a~1, t0~1, s~1, Z~1, sv~1, SZ ~ 1, st0 ~ 1),
    constants = c(s = log(1)),  # fix log(s)=0 => s=1
    model = DDM
  )

  # 2) Generate parameter means & random effects, then make data
  pmean7 <- sampled_pars(design_DDM7, doMap = FALSE)
  covs <- diag(c(.2, .075, .075, .02, .02, .02, .02, .03))
  pmean7[] <- c(1.5, 0.1, log(1), log(.435), qnorm(.5), log(.3), qnorm(.05), log(.3))
  DDM7dat <- make_data(
    make_random_effects(design_DDM7, pmean7, covariances = covs),
    design_DDM7,
    n_trials = 50
  )
  #
  emc2_time7 <- system.time({
    DDM7_emc2 <- make_emc(DDM7dat, design_DDM7, n_chains = 4, compress = FALSE, type = "diagonal")
    DDM7_emc2 <- fit(DDM7_emc2, cores_per_chain = 10)
  })[3]
  emc2_ESS7_mu <- min(ess_summary(DDM7_emc2, selection = "mu"))
  emc2_ESS7_sigma <- min(ess_summary(DDM7_emc2, selection = "sigma2"))
  emc2_ESS7_alpha <- min(ess_summary(DDM7_emc2, selection = "alpha"))
  #

  # First timing storage
  timing_results$emc2_time_ddm5 <- emc2_time5
  timing_results$stan_time_ddm5 <- stan_time5
  timing_results$emc2_time_ddm6 <- emc2_time6
  timing_results$stan_time_ddm6 <- stan_time6
  timing_results$emc2_time_ddm7 <- emc2_time7
  # Now ESS, again could have been cleaner.
  ESS_results$emc2_ESS_ddm5_mu <- emc2_ESS5_mu
  ESS_results$emc2_ESS_ddm5_sigma <- emc2_ESS5_sigma
  ESS_results$emc2_ESS_ddm5_alpha <- emc2_ESS5_alpha
  ESS_results$stan_ESS_ddm5_mu <- stan_ESS5_mu
  ESS_results$stan_ESS_ddm5_sigma <- stan_ESS5_sigma
  ESS_results$stan_ESS_ddm5_alpha <- stan_ESS5_alpha
  ESS_results$emc2_ESS_ddm6_mu <- emc2_ESS6_mu
  ESS_results$emc2_ESS_ddm6_sigma <- emc2_ESS6_sigma
  ESS_results$emc2_ESS_ddm6_alpha <- emc2_ESS6_alpha
  ESS_results$stan_ESS_ddm6_mu <- stan_ESS6_mu
  ESS_results$stan_ESS_ddm6_sigma <- stan_ESS6_sigma
  ESS_results$stan_ESS_ddm6_alpha <- stan_ESS6_alpha
  ESS_results$emc2_ESS_ddm7_mu <- emc2_ESS7_mu
  ESS_results$emc2_ESS_ddm7_sigma <- emc2_ESS7_sigma
  ESS_results$emc2_ESS_ddm7_alpha <- emc2_ESS7_alpha
  return(list(timing = timing_results, ESS = ESS_results))
}
out <- list()
for(i in 1:n_reps){
  out[[i]] <- fit_models(i)
  gc()
}

timing_hier <- do.call(rbind, lapply(out, \(x) x$timing))
ess_hier <- do.call(rbind, lapply(out, \(x) x$ESS))

save(timing_hier, file = "results/timing_hier.RData")
save(ess_hier, file = "results/ess_hier")
