  __host__ __device__ __ext_inline__
  double XCAT(XCAT(XCAT(XCAT(continua_core_,CO2),H2O),N2),O2)
  (ctl_t const *ctl, pos_t const *los, int const ig_co2, int const ig_h2o, int const id) {
      double const p = los->p;
      double const t = los->t;
      double const ds = los->ds;
      double beta_ds = los->k[ctl->window[id]]*ds;													// extinction
      // make sure that ig_co2 and ig_h2o are both >= 0
      if ((CO2)) beta_ds += continua_ctmco2(ctl->nu[id], p, t, los->u[ig_co2]);						// co2 continuum
      if ((H2O)) beta_ds += continua_ctmh2o(ctl->nu[id], p, t, los->q[ig_h2o], los->u[ig_h2o]);		// h2o continuum
      if ((N2))  beta_ds += continua_ctmn2(ctl->nu[id], p, t)*ds;									// n2 continuum
      if ((O2))  beta_ds += continua_ctmo2(ctl->nu[id], p, t)*ds;									// o2 continuum
      return     beta_ds;
  } // continua_core_bbbb where each b is either 0 or 1
