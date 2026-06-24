functions {
  real partial_sum_fullddm(array[] real rt_slice, int start, int end,
    real a, real t0_m, real t0_s, real zr_m, real zr_s, array[] real v_m, real v_s,
    array[] int resp, array[] int cnd) {
      real ans = 0;
      for (i in start:end) {
        if (resp[i] == 1) {
          // upper threshold
          ans += wiener_lpdf(rt_slice[i+1-start] | a, t0_m, zr_m, v_m[cnd[i]] , v_s,  zr_s, t0_s);
        } else {
          // lower threshold (mirror drift and starting point!)
          ans += wiener_lpdf(rt_slice[i+1-start] | a, t0_m, 1 - zr_m, -v_m[cnd[i]], v_s,  zr_s, t0_s);
        }
      }
      return ans;
    }
}

data {
  int<lower=0> N;                     // No trials
  int<lower=1> Ncnds;                 // No conditions
  array[N] real rt;                       // response times (seconds)
  array[N] int<lower=1, upper=Ncnds> cnd;       // stimulus type/condition
  array[N] int<lower=0, upper=1> resp;      // responses (0,1)
  int parallel;
}

parameters {
  real a_tr;                    // threshold separation (a>0)
  real zr_m_tr;        // mean starting point (0<zr<1)
  real zr_s_tr;  // sd of starting point (zr_s>0)

  array[Ncnds] real v_m_tr;                        // mean drift for Ncnds stimulus types
  real v_s_tr;                  // sd of drift
  //real t0_s_tr;     // sd of non-decision-time
  real t0_m_tr;                 // mean non-decision-time
}

transformed parameters{
  real a = exp(a_tr);
  real zr_m = Phi_approx(zr_m_tr);
  real zr_s = Phi_approx(zr_s_tr);
  array[Ncnds] real v_m = v_m_tr; // no transformation here, so this could be removed
  real v_s = exp(v_s_tr);
  //real t0_s = exp(t0_s_tr);
  real t0_m = exp(t0_m_tr);
}

model {
  a_tr ~ normal(0,1);
  zr_m_tr ~ normal(0,1);
  zr_s_tr ~ normal(0,1);
  v_m_tr ~ normal(0,1);
  v_s_tr ~ normal(0,1);
  // t0_s_tr ~ normal(0,1);
  t0_m_tr ~ normal(0,1);
  
  if (parallel) {
    target += reduce_sum(partial_sum_fullddm, rt, 1,
      a, t0_m, 0, zr_m, zr_s, v_m, v_s, resp, cnd); // st0 set to zero
  } else {
    for (i in 1:N) {
      if (resp[i] == 1) {
        // upper threshold
        target += wiener_lpdf(rt[i] | a, t0_m, zr_m, v_m[cnd[i]], v_s,  zr_s, 0); // st0 set to zero
      } else {
        // lower threshold (mirror drift and starting point!)
        target += wiener_lpdf(rt[i] | a, t0_m, 1-zr_m, -v_m[cnd[i]], v_s,  zr_s, 0); // st0 set to zero
      }
    }
  }
}
