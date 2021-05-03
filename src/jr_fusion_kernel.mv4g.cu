    void __global__ // GPU-kernel
    XCAT(XCAT(XCAT(XCAT(fusion_kernel_GPU_,CO2),H2O),N2),O2)
        (tbl_t const *tbl, ctl_t const *ctl,
        obs_t *obs, pos_t const (*restrict los)[NLOS],
        int const np[], int const ig_co2, int const ig_h2o) {
			double tau_path[NG];
			for(int ir = blockIdx.x; ir < obs->nr; ir += gridDim.x) { // grid stride loop over blocks = rays
				for(int id = threadIdx.x; id < ND; id += blockDim.x) { // grid stride loop over threads = detectors
					obs->rad[ir][id] = 0.0;
					obs->tau[ir][id] = 1.0;
				} // id
				for(int ig = 0; ig < NG; ++ig) {
					tau_path[ig] = 1.0;
				} // ig
				for(int ip = 0; ip < np[ir]; ++ip) {
					for(int id = threadIdx.x; id < ctl->nd; id += blockDim.x) { // grid stride loop over threads = detectors
						double const beta_ds = XCAT(XCAT(XCAT(XCAT(continua_core_,CO2),H2O),N2),O2) // function name
                                               (ctl, &(los[ir][ip]), ig_co2, ig_h2o, id);           // function args
						double const tau_gas = apply_ega_core(tbl, &(los[ir][ip]), tau_path, ctl->ng, id);
						double const planck = src_planck_core(tbl, los[ir][ip].t, id);
						new_obs_core(obs, ir, id, beta_ds, planck, tau_gas);
					} // id --> parallel over detectors=threads
				} // ip --> non-parallelisable
			} // ir --> parallel over rays==blocks
    } // fusion_kernel_GPU
