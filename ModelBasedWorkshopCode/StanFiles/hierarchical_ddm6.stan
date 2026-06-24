functions {
  // This function sums the log likelihood contributions for a slice of trials.
  // Each trial i uses the subject-specific parameters:
  //   a[s]   = threshold,
  //   t0[s]  = non-decision time,
  //   zr[s]  = starting point,
  //   v_m[s,cond] = drift,
  //   v_s[s] = drift variability,
  //   s_z[s] = starting point variability.
  real partial_sum_ddm(array[] real rt_slice, int start, int end,
                       vector a, vector t0, vector zr, matrix v_m, vector v_s, vector s_z,
                       array[] int resp, array[] int cnd, array[] int subj) {
    real ans = 0;
    for (i in start:end) {
      int s = subj[i];       // subject index for trial i
      int cond = cnd[i];     // condition index for trial i
      if (resp[i] == 1)
        ans += wiener_lpdf(rt_slice[i+1-start] | a[s], t0[s], zr[s],
                           v_m[s, cond], v_s[s], s_z[s], 0);
      else
        ans += wiener_lpdf(rt_slice[i+1-start] | a[s], t0[s], 1 - zr[s],
                           -v_m[s, cond], v_s[s], s_z[s], 0);
    }
    return ans;
  }
}

data {
  int<lower=0> N;                        // total number of trials
  int<lower=1> Nsubj;                    // number of subjects
  int<lower=1> Ncnds;                    // number of conditions
  array[N] int<lower=1, upper=Nsubj> subj; // subject index for each trial
  array[N] int<lower=1, upper=Ncnds> cnd;   // condition index for each trial
  array[N] int<lower=0, upper=1> resp;      // response on each trial (0 or 1)
  array[N] real<lower=0> rt;                // response times (seconds)
  int<lower=0, upper=1> parallel;           // flag for parallel processing (1 = yes)
}

parameters {
  // === Group-level hyperparameters ===
  // For threshold (a)
  real mu_a_tr;
  real<lower=0> sigma_a_tr;
  
  // For starting point (zr)
  real mu_zr_m_tr;
  real<lower=0> sigma_zr_m_tr;
  
  // For drift mean (v) across conditions
  vector[Ncnds] mu_v_m_tr;
  vector<lower=0>[Ncnds] sigma_v_m_tr;
  
  // For drift variability (v_s)
  real mu_v_s_tr;
  real<lower=0> sigma_v_s_tr;
  
  // For non-decision time (t0) on the log scale
  real mu_t0_m_tr;
  real<lower=0> sigma_t0_m_tr;
  
  // For starting point variability on the logit scale
  real mu_zr_sd_tr;
  real<lower=0> sigma_zr_sd_tr;
  
  // === Individual-level (non-centered) parameters ===
  vector[Nsubj] a_tr_raw;
  vector[Nsubj] zr_m_tr_raw;
  matrix[Nsubj, Ncnds] v_m_tr_raw;
  vector[Nsubj] v_s_tr_raw;
  vector[Nsubj] t0_m_tr_raw;
  vector[Nsubj] zr_sd_tr_raw;
}

transformed parameters {
  // (a) Threshold: a = exp(a_tr)
  vector[Nsubj] a_tr = mu_a_tr + sigma_a_tr * a_tr_raw;
  vector[Nsubj] a = exp(a_tr);
  
  // (b) Starting point: zr = Phi_approx(zr_m_tr) ensures zr ∈ (0, 1)
  vector[Nsubj] zr_m_tr = mu_zr_m_tr + sigma_zr_m_tr * zr_m_tr_raw;
  vector[Nsubj] zr = Phi_approx(zr_m_tr);
  
  // (c) Drift mean: v_m for each subject and condition
  matrix[Nsubj, Ncnds] v_m;
  {
    matrix[Nsubj, Ncnds] v_m_tr;
    for (i in 1:Nsubj)
      for (j in 1:Ncnds)
        v_m_tr[i, j] = mu_v_m_tr[j] + sigma_v_m_tr[j] * v_m_tr_raw[i, j];
    v_m = v_m_tr;
  }
  
  // (d) Drift variability: v_s = exp(v_s_tr)
  vector[Nsubj] v_s_tr = mu_v_s_tr + sigma_v_s_tr * v_s_tr_raw;
  vector[Nsubj] v_s = exp(v_s_tr);
  
  // (e) Non-decision time: t0 = exp(t0_m_tr) so that t0 > 0 with no upper bound.
  vector[Nsubj] t0_m_tr = mu_t0_m_tr + sigma_t0_m_tr * t0_m_tr_raw;
  vector[Nsubj] t0 = exp(t0_m_tr);
  
  // (f) Starting point variability: s_z = Phi_approx(zr_sd_tr)
  vector[Nsubj] zr_sd_tr = mu_zr_sd_tr + sigma_zr_sd_tr * zr_sd_tr_raw;
  vector[Nsubj] s_z = Phi_approx(zr_sd_tr);
}

model {
  // === Group-level priors ===
  mu_a_tr ~ normal(0, 1);
  sigma_a_tr ~ student_t(2, 0, 0.3);
  
  mu_zr_m_tr ~ normal(0, 1);
  sigma_zr_m_tr ~ student_t(2, 0, 0.3);
  
  mu_v_m_tr ~ normal(0, 1);
  sigma_v_m_tr ~ student_t(2, 0, 0.3);
  
  mu_v_s_tr ~ normal(0, 1);
  sigma_v_s_tr ~ student_t(2, 0, 0.3);
  
  mu_t0_m_tr ~ normal(0, 1);
  sigma_t0_m_tr ~ student_t(2, 0, 0.3);
  
  mu_zr_sd_tr ~ normal(0, 1);
  sigma_zr_sd_tr ~ student_t(2, 0, 0.3);
  
  // === Individual-level (non-centered) priors ===
  a_tr_raw ~ std_normal();
  zr_m_tr_raw ~ std_normal();
  to_vector(v_m_tr_raw) ~ std_normal();
  v_s_tr_raw ~ std_normal();
  t0_m_tr_raw ~ std_normal();
  zr_sd_tr_raw ~ std_normal();
  
  // === Likelihood ===
  if (parallel == 1) {
    target += reduce_sum(partial_sum_ddm, rt, 1,
                         a, t0, zr, v_m, v_s, s_z, resp, cnd, subj);
  } else {
    for (i in 1:N) {
      int s = subj[i];
      int cond = cnd[i];
      if (resp[i] == 1)
        target += wiener_lpdf(rt[i] | a[s], t0[s], zr[s],
                              v_m[s, cond], v_s[s], s_z[s], 0);
      else
        target += wiener_lpdf(rt[i] | a[s], t0[s], 1 - zr[s],
                              -v_m[s, cond], v_s[s], s_z[s], 0);
    }
  }
}
