#include "jr_common.h" // inline definitions of common functions for CPU and GPU code

	// ################ CPU driver routines - keep consistent with GPUdrivers.cu ##############

    __host__
	void radiance_to_brightness_CPU(ctl_t const *ctl, obs_t *obs) {
#pragma omp for
		for(int ir = 0; ir < obs->nr; ir++) { // loop over rays
			for(int id = 0; id < ctl->nd; id++) { // loop over detectors
				// convert in-place
				obs->rad[ir][id] = brightness_core(obs->rad[ir][id], ctl->nu[id]);
			} // id
		} // ir
	} // radiance_to_brightness

	__host__
	void surface_terms_CPU(tbl_t const *tbl, obs_t *obs, double const tsurf[], int const nd) {
#pragma omp for
		for(int ir = 0; ir < obs->nr; ir++) { // loop over rays
			for(int id = 0; id < nd; id++) { // loop over detectors
				add_surface_core(obs, tbl, tsurf[ir], ir, id);
			} // id
		} // ir
	} // surface_terms_CPU

	__host__
	double multi_continua_CPU(char const fourbit, ctl_t const *ctl, 
        pos_t const *los, int const ig_co2, int const ig_h2o, int const id) {
		switch (fourbit) { // multi-versioning for continua
			case  0: return continua_core_0000(ctl, los, ig_co2, ig_h2o, id);
			case  1: return continua_core_0001(ctl, los, ig_co2, ig_h2o, id);
			case  2: return continua_core_0010(ctl, los, ig_co2, ig_h2o, id);
			case  3: return continua_core_0011(ctl, los, ig_co2, ig_h2o, id);
			case  4: return continua_core_0100(ctl, los, ig_co2, ig_h2o, id);
			case  5: return continua_core_0101(ctl, los, ig_co2, ig_h2o, id);
			case  6: return continua_core_0110(ctl, los, ig_co2, ig_h2o, id);
			case  7: return continua_core_0111(ctl, los, ig_co2, ig_h2o, id);
			case  8: return continua_core_1000(ctl, los, ig_co2, ig_h2o, id);
			case  9: return continua_core_1001(ctl, los, ig_co2, ig_h2o, id);
			case 10: return continua_core_1010(ctl, los, ig_co2, ig_h2o, id);
			case 11: return continua_core_1011(ctl, los, ig_co2, ig_h2o, id);
			case 12: return continua_core_1100(ctl, los, ig_co2, ig_h2o, id);
			case 13: return continua_core_1101(ctl, los, ig_co2, ig_h2o, id);
			case 14: return continua_core_1110(ctl, los, ig_co2, ig_h2o, id);
			case 15: return continua_core_1111(ctl, los, ig_co2, ig_h2o, id);
		} // fourbit
		assert(0); return 0; // unreached return statement, but suppresses a compiler warning
	} // multi_continua_CPU

	__host__
	void apply_kernels_CPU(tbl_t const *tbl, ctl_t const *ctl, obs_t *obs, 
        pos_t const (*restrict los)[NLOS], int const np[], 
        int const ig_co2, int const ig_h2o, char const fourbit) {

#pragma omp for
		for(int ir = 0; ir < obs->nr; ir++) { // loop over independent rays
			double tau_path[ND][NG]; // private for each ray
			for(int id = 0; id < ND; id++) { // loop over detectors
				obs->rad[ir][id] = 0.0;
				obs->tau[ir][id] = 1.0;
				for(int ig = 0; ig < NG; ig++) { // loop over gases
					tau_path[id][ig] = 1.0;
				} // ig
			} //  id

			for(int ip = 0; ip < np[ir]; ++ip) { // loop over line-of-sight points
				for(int id = 0; id < ctl->nd; id++) { // loop over detector channels

					// compute extinction coefficient
					double const beta_ds = multi_continua_CPU(fourbit, ctl, &(los[ir][ip]), ig_co2, ig_h2o, id);
					// compute transmission with the EGA method
					double const tau_gas = apply_ega_core(tbl, &(los[ir][ip]), tau_path[id], ctl->ng, id);
					// compute the source term
					double const planck = src_planck_core(tbl, los[ir][ip].t, id);
					// perform integration
					new_obs_core(obs, ir, id, beta_ds, planck, tau_gas);

				} // id --> could be vectorized over detector channels
#ifdef GPUDEBUG
				assert(los[ir][ip].ip == ip); // sanity check
				assert(los[ir][ip].ir == ir); // sanity check
#endif
			} // ip --> non-parallelisable due to loup carried dependency
		} // ir --> OpenMP parallel over rays

	} // apply_kernels_CPU

	__host__
	void raytrace_rays_CPU(ctl_t const *ctl, atm_t const *atm, obs_t *obs, 
                           pos_t los[NR][NLOS], double tsurf[], int np[]) {
#pragma omp for
		for(int ir = 0; ir < obs->nr; ir++) { // loop over rays
			np[ir] = traceray(ctl, atm, obs, ir, los[ir], &tsurf[ir]);
		} // ir
	} // raytrace_rays_CPU

	__host__
	void hydrostatic1d_CPU(ctl_t const *ctl, atm_t *atm, int const nr, int const ig_h2o) {
			if(ctl->hydz < 0) return; // Check reference height
			for(int ir = 0; ir < nr; ir++) { // Apply hydrostatic equation to individual profiles
				hydrostatic_1d_h2o(ctl, atm, 0, atm->np, ig_h2o);
			} // ir
	} // hydrostatic1d_CPU

	// ################ end of CPU driver routines ##############

	// The full forward model on the CPU ////////////////////////////////////////////
	__host__
	void formod_CPU(ctl_t const *ctl, atm_t *atm, obs_t *obs) {
        if (ctl->checkmode) {
            printf("# %s: checkmode = %d, no actual computation is performed!\n", __func__, ctl->checkmode);
            return; // do nothing here
        } // checkmode

		assert(obs);

		char mask[NR][ND];
		save_mask(mask, obs, ctl);

		tbl_t const *tbl = get_tbl(ctl);
		double *t_surf = (double*)malloc((obs->nr)*sizeof(double));
		int *np = (int*)malloc((obs->nr)*sizeof(int));
		pos_t (*los)[NLOS] = (pos_t (*)[NLOS])malloc((obs->nr)*(NLOS)*sizeof(pos_t));

		// gas absorption continua configuration
		static int ig_co2 = -999, ig_h2o = -999;
		if((ctl->ctm_h2o) && (-999 == ig_h2o)) ig_h2o = find_emitter(ctl, "H2O");
		if((ctl->ctm_co2) && (-999 == ig_co2)) ig_co2 = find_emitter(ctl, "CO2");
		// binary switches for the four gases
		char const fourbit = (char) (
                  ( (1 == ctl->ctm_co2) && (ig_co2 >= 0) )*8   // CO2
				+ ( (1 == ctl->ctm_h2o) && (ig_h2o >= 0) )*4   // H2O
				+   (1 == ctl->ctm_n2)                    *2   // N2
				+   (1 == ctl->ctm_o2)                    *1); // O2

        hydrostatic1d_CPU(ctl, atm, obs->nr, ig_h2o); // in this call atm might get modified
        raytrace_rays_CPU(ctl, atm, obs, los, t_surf, np);
#pragma omp parallel
		{
			apply_kernels_CPU(tbl, ctl, obs, los, np, ig_co2, ig_h2o, fourbit);
			surface_terms_CPU(tbl, obs, t_surf, ctl->nd);
		} // parallel

		free(los);
		free(np);
		free(t_surf);

		if(ctl->write_bbt) radiance_to_brightness_CPU(ctl, obs);

		apply_mask(mask, obs, ctl);
	} // formod_CPU

    __host__ void formod_GPU(ctl_t const *ctl, atm_t *atm, obs_t *obs)
#ifdef hasGPU
    ; // declaration only, will be provided by GPUdrivers.o at link time 
#else
    { // definition here
        static int warnGPU = 1;
        if (ctl->useGPU > 0) { // USEGPU > 0 means use-GPU-always
            fprintf(stdout, "CUDA not found during compilation, cannot run on GPUs. ABORT\n\n");
            fprintf(stdout, "USEGPU = 1 (use-GPU-always) found in controls\n"
                            "USEGPU = -1 (use-GPU-if-possible) could help\n\n");
            fflush(stdout);
            fprintf(stderr, "CUDA not found during compilation, cannot run on GPUs. ABORT\n\n");
            exit(EXIT_FAILURE);
        } else {               // USEGPU < 0 means use-GPU-if-possible
            assert(ctl->useGPU < 0); 
            // automatic decision: fallback solution is CPU
            if (warnGPU) { // this warning appears only once per process
                printf("CUDA not found during compilation, continue on CPUs instead!\n");
                warnGPU = 0; // switch this warning off
            } // warnGPU
            formod_CPU(ctl, atm, obs);
        } //
    } // formod_GPU
#endif

	__host__
	void formod(ctl_t const *ctl, atm_t *atm, obs_t *obs) {
        if (ctl->checkmode) {
            static int nr_last_time = -999;
            if (obs->nr != nr_last_time) {
                printf("# %s: %d max %d rays , %d of max %d gases, %d of max %d channels\n", 
                        __func__, obs->nr, NR, ctl->ng, NG, ctl->nd, ND);
                // the number of rays can be zero if we skipped to read obs.tab, checkmode=1
                nr_last_time = obs->nr;
            } // only report if nr changed
        } // checkmode
        if (ctl->useGPU) {
            formod_GPU(ctl, atm, obs);
        } else { // USEGPU = 0 means use-GPU-never
            formod_CPU(ctl, atm, obs);
        } // useGPU
    } // formod
