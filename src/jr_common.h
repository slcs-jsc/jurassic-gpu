#pragma once

#include <assert.h> // assert


#ifdef __NVCC__
  
  extern "C" {
  #include "jurassic.h" // ...
  }

  #define UNROLL _Pragma("unroll")
  #define __ext_inline__ inline
  #define restrict __restrict
#else

  #include "jurassic.h" // ...
  
  #define UNROLL
  // Workaround if no NVIDIA compiler is available: __host__ and __device__ are simply ignored
  #define __host__
  #define __device__
  #define __global__
  #define __ext_inline__ extern inline
#endif

#define ptr const restrict


#ifdef  FAST_INVERSE_OF_U
	__host__ __device__ __ext_inline__ 
	int fast_logarithmic_index(double const x) {
        // see checks in jurassic.c for the derivation of this
        double const x2 = x*x, x4 = x2*x2, x6 = x4*x2;
        unsigned long long const bits = *((unsigned long long*)(&x6)); // reinterpret cast
        int const approximate_log2_of_x6_times8 = (bits >> (52 - 3)) - (1023 << 3);
        int const iu_fast = (int)(1.003472 * 0.125 * approximate_log2_of_x6_times8);
    } // fast_logarithmic_index
#endif


	// Clamp function ////////////////////////////////////////////////////////////
	__host__ __device__ __ext_inline__ 
	double c01(double const x) 
    {   return (x > 1.)?1.:((x < 0.)?0.:x); }

	// Linear interpolation //////////////////////////////////////////////////////
	__host__ __device__ __ext_inline__
	double lip(double const x0, double const y0, double const x1, double const y1, double const x) 
    {   return y0 + (x - x0)*(y1 - y0)/(x1 - x0); }

	// Exponential interpolation //////////////////////////////////////////////////
	__host__ __device__ __ext_inline__
	double eip(double const x0, double const y0, double const x1, double const y1, double const x) {   
        if((y0 > 0) && (y1 > 0)) return y0*exp(log(y1/y0)/(x1 - x0)*(x - x0));
	    else                     return lip(x0, y0, x1, y1, x);
	} // eip

	// Setup tables if necessary and cache them ///////////////////////////////////
	__host__ __ext_inline__
	tbl_t* get_tbl(ctl_t const *ctl) {   
        static tbl_t *tbl = NULL;
		if(!tbl) {
#pragma omp barrier
#pragma omp master
			{
#ifdef  USE_UNIFIED_MEMORY_FOR_TABLES
            printf("# call cudaMallocManaged for tables of size %.3f MByte\n", 1e-6*sizeof(tbl_t));
			int const status = cudaMallocManaged(&tbl, sizeof(tbl_t));
#else
			tbl = (tbl_t*)malloc(sizeof(tbl_t));
#endif            
			init_tbl(ctl, tbl); // CPU reads the table content from files
			}
#pragma omp barrier
		}
		return tbl;
	} // get_tbl

	// Index finding ////////////////////////////////////////////////////////////
#pragma GCC diagnostic ignored "-Wunused-parameter"
	__host__ __device__ __ext_inline__
	int locate_st(double const *ptr xx, int const n, double const x)
    {   return (int)(4*x) - 400; } // only for source temperatures

	// Table lookups ////////////////////////////////////////////////////////////
	__host__ __device__ __ext_inline__
	int locate(double const *ptr xx, int const n, double const x) {
        int ilo = 0, ihi = n - 1, i = (n - 1) >> 1;
        if(xx[i] < xx[i + 1]) {
            while (ihi > ilo + 1) { // divergent execution on GPU happens here
                i = (ihi + ilo) >> 1;
                if(xx[i] > x) ihi = i;
                else					ilo = i;
            } // while
        } else {
            while (ihi > ilo + 1) {
                i = (ihi + ilo) >> 1;
                if(xx[i] <= x) ihi = i;
                else					 ilo = i;
            } // while
        } // if
        return ilo;
    } // locate

	__host__ __device__ __ext_inline__
    int locate_id(double const (*ptr xx)[ND], int const n, double const x, int const id) {
        int ilo = 0, ihi = n - 1;
        while (ihi > ilo + 1) { // divergent execution on GPU happens here
            int i = (ihi + ilo) >> 1;
            if (xx[i][id] > x) { ihi = i; } else { ilo = i; }
        } // while
        return ilo;
    } // locate_id

	__host__ __device__ __ext_inline__
    int locate_tbl_id(real_tblND_t const (*ptr xx)[ND], int const n, double const x, int const id, int const i0) {
        int ilo = i0, ihi = n - 1;
        while (ihi > ilo + 1) { // divergent execution on GPU happens here
            int i = (ihi + ilo) >> 1;
            if (xx[i][id] > x) ihi = i;
            else ilo = i;
        } // while
        return ilo;
    } // locate_tbl_id
    
	__host__ __device__ __ext_inline__
    void locate_atm(atm_t const *atm, double const time, size_t *atmIdx, int *atmNp) {
        assert(atm);
        // Find lower bound of time stamp
//         printf("# %s (line %d) atm->np=%d\n", __func__, __LINE__, atm->np);
        int lo = 0, hi = atm->np - 1;
        int i;
        while(hi > lo + 1) {
            i = (lo + hi)/2;
//             printf("# %s %d < %d %d\n", __func__, hi, lo, i);
            if(atm->time[i] < time) lo = i; else hi = i;
        } // while
        int const lower = (0 == lo) ? lo : hi;
        *atmIdx = (unsigned)lower;

        // Find upper bound
        lo = lower;
        hi = atm->np - 1;
        while(hi > lo + 1) {
            i = (lo + hi)/2;
//             printf("# %s %d > %d %d\n", __func__, hi, lo, i);
            if(atm->time[i] > time) hi = i; else lo = i;
        } // while
        int const upper = (hi == atm->np - 1) ? atm->np : hi;
        *atmNp = upper - lower;
        
//         printf("# %s [%d, %d)\n\n", __func__, lower, upper);
    } // locate_atm

	__host__ __device__ __ext_inline__
    double get_eps(tbl_t const *tbl, int const ig, int const id, int const ip, int const it, double const u) {
        int const nu = tbl->nu[ig][ip][it][id]; // number of u grid entries
#ifdef  FAST_INVERSE_OF_U
        double const x = u * tbl->u0inv[ig][ip][it][id];
        int const ifx = (x > 1) ? fast_logarithmic_index(x) : 0;
        int const ifxc = (ifx < nu) ? ifx : (nu - 1);
        
        int const idx_ref = locate_tbl_id(tbl->u[ig][ip][it], nu, u, id, 0); // DEBUG
        if (abs(idx_ref - ifxc) > 1) printf("# FAST_INVERSE_OF_U locate= %i fast= %i\n", idx_ref, ifxc);
        
        int const guess_0 = (ifxc > 0)?(ifxc - 1):0;
        int const guess_n = ((ifxc + 2) > nu)?nu:(ifxc + 2);
        int const idx = locate_tbl_id(tbl->u[ig][ip][it], guess_n, u, id, guess_0);
        assert(idx == idx_ref); // DEBUG
#else
        int const idx = locate_tbl_id(tbl->u[ig][ip][it], nu, u, id, 0);
#endif
        return lip(tbl->u[ig][ip][it][idx		][id], tbl->eps[ig][ip][it][idx		][id],
                tbl->u[ig][ip][it][idx + 1][id], tbl->eps[ig][ip][it][idx + 1][id],
                u); // <- e_i + (e_i+1 - e_i)(u - u_i)/(u_i+1 - u_i)
    } // get_eps

	__host__ __device__ __ext_inline__
    double get_u(tbl_t const *tbl, int const ig, int const id, int const ip, int const it, double const eps) {
        int const idx = locate_tbl_id(tbl->eps[ig][ip][it], tbl->nu[ig][ip][it][id], eps, id, 0);
        return lip(tbl->eps[ig][ip][it][idx		][id], tbl->u[ig][ip][it][idx		][id],
                tbl->eps[ig][ip][it][idx + 1][id], tbl->u[ig][ip][it][idx + 1][id],
                eps);
    } // get_u

	// Convert radiance to brightness ////////////////////////////////////////////
	__host__ __device__ __ext_inline__
    double brightness_core(double const rad, double const nu)
    {   return C2*nu/log1p((C1*nu*nu*nu)/rad); }

	// Save observation mask prior to formod ////////////////////////////////////
	__host__ __ext_inline__
    void save_mask(char mask[NR][ND], obs_t const *obs, ctl_t const *ctl) {
		for(int ir = 0; ir < obs->nr; ir++) {
			for(int id = 0; id < ctl->nd; id++) {
				mask[ir][id] = !gsl_finite(obs->rad[ir][id]);
			} // id
		} // ir
	} // save_mask

	// Apply observation mask after formod ///////////////////////////////////////
	__host__ __ext_inline__
	void apply_mask(char mask[NR][ND], obs_t *obs, ctl_t const *ctl) {
		for(int ir = 0; ir < obs->nr; ir++) {
			for(int id = 0; id < ctl->nd; id++) {
				if (mask[ir][id]) obs->rad[ir][id] = GSL_NAN;
			} // id
		} // ir
	} // apply_mask

	// Gravity as a function of altitude and latitude
	__host__ __device__ __ext_inline__
    double gravity(double const z, double const lat) {
        double const deg2rad = M_PI/180., x = sin(lat*deg2rad), y = sin(2*lat*deg2rad);
        return 9.780318*(1. + 0.0053024*x*x - 5.8e-6*y*y) - 3.086e-3*z;
    } // gravity

	// Black body radiation //////////////////////////////////////////////////////
	__host__ __device__ __ext_inline__
    double src_planck_core(tbl_t const *tbl, double const t, int const id) {
        int const it = locate_st(tbl->st, TBLNS, t);
        return lip(tbl->st[it], tbl->sr[it][id], tbl->st[it + 1], tbl->sr[it + 1][id], t);
    } // src_planck_core

	// Surface emission //////////////////////////////////////////////////////////
	__host__ __device__ __ext_inline__ 
    void add_surface_core(obs_t *obs, tbl_t const *tbl, double const tsurf, int const ir, int const id) {
        if(tsurf > 0.) {
            int const it = locate_st(tbl->st, TBLNS, tsurf);
            double const src = lip(tbl->st[it], tbl->sr[it][id], tbl->st[it + 1], tbl->sr[it + 1][id], tsurf);
            obs->rad[ir][id] += src*obs->tau[ir][id];
        } // if
    } // add_surface_core

	// EGA model //////////////////////////////////////////////////////////////////
	__host__ __device__ __ext_inline__
    double ega_eps(tbl_t const *tbl, double const tau, double const t, double const u, double const p, int const ig, int const id) {
        if(tau < 1e-9)  return 0.; // opaque
        if(tbl->np[ig][id] < 2) return 1.; // no table
        int const ipr = locate_id(tbl->p[ig], tbl->np[ig][id], p, id);
        if(tbl->nt[ig][ipr		][id]			 < 2 || tbl->nt[ig][ipr + 1][id]					< 2) return 1.;
        int const it0 = locate_id(tbl->t[ig][ipr		 ], tbl->nt[ig][ipr		 ][id], t, id);
        if(tbl->nu[ig][ipr		][it0][id] < 2 || tbl->nu[ig][ipr		 ][it0 + 1][id] < 2) return 1.;
        int const it1 = locate_id(tbl->t[ig][ipr + 1], tbl->nt[ig][ipr + 1][id], t, id);
        if(tbl->nu[ig][ipr + 1][it1][id] < 2 || tbl->nu[ig][ipr + 1][it1 + 1][id] < 2) return 1.;

        double const eps = 1 - tau;
        double const u00 = get_u(tbl, ig, id, ipr,		 it0,			eps);
        double const u01 = get_u(tbl, ig, id, ipr,		 it0 + 1, eps);
        double const u10 = get_u(tbl, ig, id, ipr + 1, it1,			eps);
        double const u11 = get_u(tbl, ig, id, ipr + 1, it1 + 1, eps);

        double const eps00 = c01(get_eps(tbl, ig, id, ipr,		 it0,			u00 + u));
        double const eps01 = c01(get_eps(tbl, ig, id, ipr,		 it0 + 1, u01 + u));
        double const eps10 = c01(get_eps(tbl, ig, id, ipr + 1, it1,			u10 + u));
        double const eps11 = c01(get_eps(tbl, ig, id, ipr + 1, it1 + 1, u11 + u));

        double const eps_p0 = c01(lip(tbl->t[ig][ipr		 ][it0		][id], eps00,
                    tbl->t[ig][ipr		 ][it0 + 1][id], eps01, t));
        double const eps_p1 = c01(lip(tbl->t[ig][ipr + 1][it1		][id], eps10,
                    tbl->t[ig][ipr + 1][it1 + 1][id], eps11, t));

        double const eps_t	= c01(lip(tbl->p[ig][ipr		 ][id], eps_p0,
                    tbl->p[ig][ipr + 1][id], eps_p1, p));

        return (1. - eps_t)/tau; // if divisions are expensive and ng > 1, it could be collected for all gases first...
    } // ega_eps

	__host__ __device__ __ext_inline__
    double apply_ega_core(tbl_t const *tbl, pos_t const *los, double (*ptr tau_path), int const ng, int const id) {
        double tau_gas = 1.0;
        for(int ig = 0; ig < NG; ig++) { // to enable unrolling of this loop, static indexing into tau_path, so tau_path can stay in regfile
            double eps = 1.0;
            if (ig < ng) eps = ega_eps(tbl, tau_path[ig], los->t, los->u[ig], los->p, ig, id);
            tau_path[ig] *= eps;
            tau_gas      *= eps;
        } // ig
        return tau_gas;
    } // apply_ega_core

	__host__ __device__ __ext_inline__
    void apply_ega_kernel(tbl_t const *tbl, pos_t const *los,
            double (*ptr tau_path)[NG], // tau_path[id][ig] gets modified as well
            double *ptr tau_gas, // [ND] result
            int const ng, int const nd) {
        for(int id = 0; id < nd; id++) {
            tau_gas[id] = apply_ega_core(tbl, los, tau_path[id], ng, id);
        } // id
    } // apply_ega_kernel

	// Update observation struct /////////////////////////////////////////////////
	__host__ __device__ __ext_inline__
    void new_obs_core(obs_t *obs, int const ir, int const id, double const beta_ds, double const src, double const tau_gas) {
        if (tau_gas > 1e-50) {
            double const eps = 1. - tau_gas*exp(-beta_ds);
            obs->rad[ir][id] += src*eps*obs->tau[ir][id];
            obs->tau[ir][id] *= (1. - eps); // update bar tau
        } // if
    } // new_obs_core

	// Continuum //////////////////////////////////////////////////////////////

    #define LOC(nm) nm

    __host__ __device__ __ext_inline__ 
    double load_ro(double const *ptr arr, int const idx) {
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 350
        return __ldg(arr + idx);
#else
        return arr[idx];
#endif
    } // load_ro

    __host__ __device__ __ext_inline__
    double continua_ctmco2(double const nu, double const p, double const t, double const u) {
        #include "ctmco2.tbl"
        if (nu < 0 || nu >= 4000) return 0;
        double const xw = nu*0.5 + 1;
        int    const iw = (int) xw;
        double const dw = xw - iw;
        double const ew = 1 - dw;
        double const cw296 = ew*load_ro(LOC(co2296), iw - 1) + dw*load_ro(LOC(co2296), iw);
        double const cw260 = ew*load_ro(LOC(co2260), iw - 1) + dw*load_ro(LOC(co2260), iw);
        double const cw230 = ew*load_ro(LOC(co2230), iw - 1) + dw*load_ro(LOC(co2230), iw);
        double const dt230 = t - 230;
        double const dt260 = t - 260;
        double const dt296 = t - 296;
        double const ctw = dt260*5.050505e-4*dt296*cw230 - dt230*9.259259e-4*dt296*cw260 + dt230*4.208754e-4*dt260*cw296;
        return u*p*ctw/(GSL_CONST_NUM_AVOGADRO*1000*P0);
    } // ctmco2

    __host__ __device__ __ext_inline__
    double continua_ctmh2o(double const nu, double const p, double const t, double const q, double const u) {
        #include "ctmh2o.tbl"
        if (nu < 0 || nu >= 20000) return 0;
        double const xw = nu/10 + 1;
        int    const iw = (int) xw;
        double const dw = xw - iw;
        double const ew = 1 - dw;
        double const cw296 = ew*load_ro(LOC(h2o296), iw - 1) + dw*load_ro(LOC(h2o296), iw);
        double const cw260 = ew*load_ro(LOC(h2o260), iw - 1) + dw*load_ro(LOC(h2o260), iw);
        double const cwfrn = ew*load_ro(LOC(h2ofrn), iw - 1) + dw*load_ro(LOC(h2ofrn), iw);
        double sfac = 1.;
        if ((nu > 820.) && (nu < 960.)) { // equidistant grid of 10 cm^-1
            char const xfcrev_char[16] = {3, 9, 15, 23, 29, 33, 37, 39, 40, 46, 36, 27, 10, 2, 0, 0};
            float const xx = nu*0.1 - 82; // xx = (nu - 820)/10.;
            int   const ix = (int)xx;
            float const dx = xx - ix;
            sfac += .001*((1 - dx)*xfcrev_char[ix] + dx*xfcrev_char[ix + 1]);
        }
        double const ctwslf = sfac*cw296*pow(cw260/cw296, (296. - t)/(296. - 260.));
        double const vf1 = nu - 370.;
        double const vf2 = vf1*vf1;
        double const vf6 = vf2*vf2*vf2;
        double const fscal = 36100./(vf2 + vf6*1e-8 + 36100.)*-.25 + 1.;
        double const ctwfrn = cwfrn*fscal;
        double const a1 = nu*u*tanh(.7193876/t*nu);
        double const a2 = 296./t;
        double const a3 = p/P0*(q*ctwslf + (1 - q)*ctwfrn)*1e-20;
        return a1*a2*a3;
    } // ctmh2o

    __host__ __device__ __ext_inline__
    double continua_ctmn2(double const nu, double const p, double const t) {
        #include "ctmn2.tbl"
        if (nu < 2120 || nu > 2605) return 0;
        double const xnu = nu*0.2 - 424; // 2120/5 = 424, to be exact, use xnu = (nu - 2120.)/5.
        int    const idx = (int) xnu;
        double const a1 = xnu - idx, a0 = 1 - a1;
        double const b		= a0*load_ro(LOC(ba),		 idx) + a1*load_ro(LOC(ba),		 idx + 1);
        double const beta = a0*load_ro(LOC(betaa), idx) + a1*load_ro(LOC(betaa), idx + 1);
        double const q_n2 = 0.79, t0 = 273, tr = 296;
        // Compute absorption coefficient
        return 0.1*(p/P0)*(p/P0)*(t0/t)*(t0/t)*exp(beta*(1/tr - 1/t))*q_n2*b*(q_n2 + (1 - q_n2)*(1.294 - 0.4545*t/tr));
    } // ctmn2

    __host__ __device__ __ext_inline__
    double continua_ctmo2(double const nu, double const p, double const t) {
        #include "ctmo2.tbl"
        if (nu < 1360 || nu > 1805) return 0;
        double const xnu = nu*0.2 - 272;
        int    const idx = (int) xnu;
        double const a1 = xnu - idx, a0 = 1 - a1;
        double const b		= a0*load_ro(LOC(ba),		 idx) + a1*load_ro(LOC(ba),		 idx + 1);
        double const beta = a0*load_ro(LOC(betaa), idx) + a1*load_ro(LOC(betaa), idx + 1);
        double const q_o2 = 0.21, t0 = 273, tr = 296;
        // Compute absorption coefficient
        return 0.1*(p/P0)*(p/P0)*(t0/t)*(t0/t)*exp(beta*(1/tr - 1/t))*q_o2*b;
    } // ctmo2

// template<int CO2, int H2O, int N2, int O2> for multi-versioning
#define KERNEL "jr_continua_core.mv4g.h"
      #include "jr_multiversion4gases.h" // continua_core_0000, 0001, ..., _1111
#undef  KERNEL

    __host__ __device__ __ext_inline__
    double continua_core(ctl_t const *ctl, pos_t const *los, int const ig_co2, int const ig_h2o, int const id) {
        double const p = los->p;
        double const t = los->t;
        double const ds = los->ds;
        double beta_ds = los->k[ctl->window[id]]*ds;													// extinction
        // make sure that ig_co2 and ig_h2o are both >= 0
                  beta_ds += continua_ctmco2(ctl->nu[id], p, t, los->u[ig_co2]);						// co2 continuum
                  beta_ds += continua_ctmh2o(ctl->nu[id], p, t, los->q[ig_h2o], los->u[ig_h2o]);		// h2o continuum
                  beta_ds += continua_ctmn2(ctl->nu[id], p, t)*ds;									// n2 continuum
                  beta_ds += continua_ctmo2(ctl->nu[id], p, t)*ds;									// o2 continuum
        return    beta_ds;
    } // continua_core

	__host__ __device__ __ext_inline__ 
	void altitude_range_nn(atm_t const *atm, size_t const atmIdx, int const atmNp, double *zmin, double *zmax) {
		*zmax = *zmin = atm->z[atmIdx];
		for(size_t ipp = atmIdx;
				(ipp < atmIdx+atmNp) && (atm->lon[ipp] == atm->lon[atmIdx]) && (atm->lat[ipp] == atm->lat[atmIdx]);
				++ipp) {
			*zmax = fmax(*zmax, atm->z[ipp]);
			*zmin = fmin(*zmin, atm->z[ipp]);
		} // ipp
	} // altitude_range_nn

	__host__ __device__ __ext_inline__ 
	void write_pos_point(pos_t *los,
			double const lon, double const lat, double const z,
			double const p, double const t, double const q[], double const k[], double const ds) {
		los->lon = lon;
		los->lat = lat;
		los->z = z;
		los->p = p;
		los->t = t;
		for(int ig = 0; ig < NG; ig++) los->q[ig] = q[ig];
		for(int iw = 0; iw < NW; iw++) los->k[iw] = k[iw];
		los->ds = ds;
	} // write_pos_point

	// Change segment lengths according to trapezoid rule
	__host__ __device__ __ext_inline__ 
	void trapezoid_rule_pos(int np, pos_t los[]) {
		for(int ip = np - 1; ip >= 1; ip--) {
			los[ip].ds = 0.5*(los[ip - 1].ds + los[ip].ds);
		} // ip
		los[0].ds *= 0.5; // first point
	} // trapezoid_rule

	// Compute column density
	__host__ __device__ __ext_inline__ 
	void column_density(int const ng, pos_t los[], int const np) {
		for(int ip = 0; ip < np; ip++) {
			for(int ig = 0; ig < ng; ig++) {
				los[ip].u[ig] = 10.*los[ip].q[ig]*los[ip].p/(GSL_CONST_MKSA_BOLTZMANN*los[ip].t)*los[ip].ds;
			} // ig
		} // ip
	} // column_density

#ifdef CURTIS_GODSON
	__host__ __device__ __ext_inline__ 
	void curtis_godson(ctl_t const *ctl, pos_t los[], int const np) {
		for(int ig = 0; ig < ctl->ng; ig++) { // Compute Curtis-Godson pressure and temperature
			los[0].cgp[ig] = los[0].u[ig]*los[0].p;
			los[0].cgt[ig] = los[0].u[ig]*los[0].t;
			los[0].cgu[ig] = los[0].u[ig];
			for(int ip = 1; ip < np; ip++) {
				los[ip].cgp[ig] = los[ip - 1].cgp[ig] + los[ip].u[ig]*los[ip].p;
				los[ip].cgt[ig] = los[ip - 1].cgt[ig] + los[ip].u[ig]*los[ip].t;
				los[ip].cgu[ig] = los[ip - 1].cgu[ig] + los[ip].u[ig];
			} // ip
			for(int ip = 0; ip < np; ip++) {
				los[ip].cgp[ig] /= los[ip].cgu[ig];
				los[ip].cgt[ig] /= los[ip].cgu[ig];
			} // ip
		} // ig
	} // curtis_godson
#endif // CURTIS_GODSON

	__host__ __device__ __ext_inline__ 
	double refractivity(double const p, double const t)
    {   return 7.753e-05*p/t; }

#define rad2grd (180/M_PI)
#define grd2rad (M_PI/180)

	__host__ __device__ __ext_inline__ 
	void cart2geo(double const x[], double *alt, double *lon, double *lat) {
		double const radius = NORM(x);
		*lat = asin(x[2]/radius)*rad2grd;
		*lon = atan2(x[1], x[0])*rad2grd;
		*alt = radius - RE; // subtract radius of the earth
	} // cart2geo

	__host__ __device__ __ext_inline__ 
	double cart2alt(double const x[]) 
    {   return NORM(x) - RE; } // compute the altitude only

	__host__ __device__ __ext_inline__
	void geo2cart(double const alt, double const lon, double const lat, double x[]) {
		double const radius = alt + RE, clat = cos(lat*grd2rad);
		x[0] = radius*clat*cos(lon*grd2rad);
		x[1] = radius*clat*sin(lon*grd2rad);
		x[2] = radius*sin(lat*grd2rad);
	} // geo2cart

	__host__ __device__ __ext_inline__
	void tangent_point(pos_t const los[], const int np, const int ip, double *tpz, double *tplon, double *tplat) {
		// ip (=) gsl_stats_min_index(los->z, 1, (size_t) los->np), found while tracing!
		if(ip <= 0 || ip >= np-1) {		// Nadir or zenith
			*tpz	 = los[np-1].z;
			*tplon = los[np-1].lon;
			*tplat = los[np-1].lat;
		} else {																			// Limb
			// Determine interpolating polynomial y=a*x^2+b*x+c
			double const
				yy0 = los[ip - 1].z,
                yy1 = los[ip].z,
                yy2 = los[ip + 1].z,
                ds0 = los[ip].ds,
                ds1 = los[ip + 1].ds,
                dyy10 = yy1 - yy0,
                dyy21 = yy2 - yy1,
                x1	=      sqrt(ds0*ds0 - dyy10*dyy10),
                x2	= x1 + sqrt(ds1*ds1 - dyy21*dyy21),
                dx12	= x1 - x2,
                a		= (dyy10*x2 + (yy0 - yy2)*x1)/(x1*x2*dx12),
                b		= dyy10/x1 - a*x1,
                c		= yy0,
                x		= -b/(2*a);													// Get tangent point location
            
			*tpz = (a*x + b)*x + c;
			double v[3], v0[3], v2[3], dummy;
			geo2cart(los[ip - 1].z, los[ip - 1].lon, los[ip - 1].lat, v0);
			geo2cart(los[ip + 1].z, los[ip + 1].lon, los[ip + 1].lat, v2);
//             printf("# %s v0= %g %g %g\n", __func__, v0[0], v0[1], v0[2]);
//             printf("# %s v2= %g %g %g\n", __func__, v2[0], v2[1], v2[2]);
//             printf("# %s x2= %g\n", __func__, x2);
			UNROLL
              for(int i = 0; i < 3; i++) v[i] = lip(0.0, v0[i], x2, v2[i], x);
//          printf("# %s v= %g %g %g\n", __func__, v[0], v[1], v[2]);
			cart2geo(v, &dummy, tplon, tplat);
		}
	} // tangent_point_pos

	__host__ __device__ __ext_inline__ 
	void last_point(pos_t const *los, double *tpz, double *tplon, double *tplat) {
		// Nadir sounder uses the last point as tangent point
		*tpz	 = los->z;	 // altitude
		*tplon = los->lon; // longitude
		*tplat = los->lat; // latitude
	} // last_point

	__host__ __device__ __ext_inline__ 
	void intpol_atm_1d_pt(ctl_t const *ctl, atm_t const *atm,
			int const idx0, int const n, double const z0, double *p, double *t) {
		int const ip = idx0 + locate(&atm->z[idx0], n, z0);						        // Get array index
		*p = eip(atm->z[ip], atm->p[ip], atm->z[ip + 1], atm->p[ip + 1], z0);			 // Interpolate
		*t = lip(atm->z[ip], atm->t[ip], atm->z[ip + 1], atm->t[ip + 1], z0);
	} // intpol_atm_1d_pt

	__host__ __device__ __ext_inline__ 
	void intpol_atm_1d_qk(ctl_t const *ctl, atm_t const *atm,
			int const idx0, int const n, double const z0, double q[], double k[]) {
		int const ip = idx0 + locate(&atm->z[idx0], n, z0);																	// Get array index
		for(int ig = 0; ig < ctl->ng; ig++) {
			q[ig] = lip(atm->z[ip], atm->q[ig][ip], atm->z[ip + 1], atm->q[ig][ip + 1], z0);	// Interpolate
		} // ig
		for(int iw = 0; iw < ctl->nw; iw++) {
			k[iw] = lip(atm->z[ip], atm->k[iw][ip], atm->z[ip + 1], atm->k[iw][ip + 1], z0);
		} // iw
	} // intpol_atm_1d_qk

	__host__ __device__ __ext_inline__ 
	void intpol_atm_geo_pt(ctl_t const *ctl, atm_t const *atm, int const atmIdx, int const atmNp,
			double const z0, double const lon0, double const lat0,
			double *p, double *t) {
		assert(ctl->ip == 1);
		intpol_atm_1d_pt(ctl, atm, atmIdx, atmNp, z0, p, t); // 1D interpolation (vertical profile)
	} // intpol_atm_geo_pt

	__host__ __device__ __ext_inline__ 
	void intpol_atm_geo_qk(ctl_t const *ctl, atm_t const *atm, int const atmIdx, int const atmNp,
			double const z0, double const lon0, double const lat0,
			double q[], double k[]) { // results
		assert(ctl->ip == 1);
		intpol_atm_1d_qk(ctl, atm, atmIdx, atmNp, z0, q, k); // 1D interpolation (vertical profile)
	} // intpol_atm_geo_qk

	__host__ __device__ __ext_inline__ 
	int traceray(ctl_t const *ctl, atm_t const *atm, obs_t *obs, int const ir, pos_t los[], double *tsurf) {
		double ex0[3], ex1[3], q[NG], k[NW], lat, lon, p, t, x[3], xobs[3], xvp[3], z = 1e99, z_low=z, zmax, zmin, zrefrac = 60;
		// Initialize
		*tsurf = -999;
		for(int ig = 0; ig < NG; ig++) q[ig] = 0;
		for(int iw = 0; iw < NW; iw++) k[iw] = 0;
		obs->tpz[ir]   = obs->vpz[ir];
		obs->tplon[ir] = obs->vplon[ir];
		obs->tplat[ir] = obs->vplat[ir];
		size_t atmIdx=0; int atmNp=0;
		locate_atm(atm, obs->time[ir], &atmIdx, &atmNp);
		altitude_range_nn(atm, atmIdx, atmNp, &zmin, &zmax);
		if(obs->obsz[ir] < zmin)        return 0;																		// Check observer altitude
		if(obs->vpz[ir] > zmax - 0.001) return 0;																		// Check view point altitude
		geo2cart(obs->obsz[ir], obs->obslon[ir], obs->obslat[ir], xobs);						// Cart. coordinates of observer
		geo2cart(obs->vpz[ir],	obs->vplon[ir],  obs->vplat[ir],	xvp);							// and view point
		UNROLL
			for(int i = 0; i < 3; i++) ex0[i] = xvp[i] - xobs[i];											// Determine initial tangent vector
		double const norm = NORM(ex0);
		UNROLL
			for(int i = 0; i < 3; i++) {
				ex0[i] /= norm;
				x[i] = xobs[i];																													// Observer within atmosphere
			} // i
		if(obs->obsz[ir] > zmax) {																									// Above atmosphere, search entry point
			double dmax = norm, dmin = 0.;
			while(fabs(dmin - dmax) > 0.001) {
				double const d = 0.5*(dmax + dmin);
				UNROLL
					for(int i = 0; i < 3; i++) x[i] = xobs[i] + d*ex0[i];
				z = cart2alt(x); // no need to compute lat and lon here
				if((z <= zmax) && (z > zmax - 0.001)) break;
				if(z < zmax - 0.0005) dmax = d;
				else									dmin = d;
			} // while
		}

		int np = 0, z_low_idx=-1;
		for(int stop = 0; np < NLOS; ++np) {																				// Ray-tracing
			double ds = ctl->rayds, dz = ctl->raydz;																	// Set step length
			if(dz > 0.) {
				double const norm_x = 1.0/NORM(x);
				double dot = 0.;
				UNROLL
					for(int i = 0; i < 3; i++) {
						dot += ex0[i]*x[i]*norm_x;
					}
				double const cosa = fabs(dot);
				if(cosa != 0.) ds = fmin(ds, dz/cosa);
			}
			cart2geo(x, &z, &lon, &lat);																							// Determine geolocation
			if((z < zmin) || (z > zmax)) {																						// LOS escaped
				double xh[3];
				stop = (z < zmin) ? 2 : 1;
				geo2cart(los[np - 1].z, los[np - 1].lon, los[np - 1].lat, xh);
				double const zfrac = (z < zmin) ? zmin : zmax;
				double const frac = (zfrac - los[np - 1].z)/(z - los[np - 1].z);
				UNROLL
					for(int i = 0; i < 3; i++) x[i] = xh[i] + frac*(x[i] - xh[i]);
				cart2geo(x, &z, &lon, &lat);
				los[np - 1].ds = ds*frac;
				ds = 0.;
			}

			intpol_atm_geo_pt(ctl, atm, (int) atmIdx, atmNp, z, lon, lat, &p, &t);					// Interpolate atmospheric data
			intpol_atm_geo_qk(ctl, atm, (int) atmIdx, atmNp, z, lon, lat, q, k);						// Interpolate atmospheric data
//     printf("ray #%i point#%i %g %g %g\n", ir, np, lon, lat, z); // to debug the raytracer
			write_pos_point(los + np, lon, lat, z, p, t, q, k, ds);
			if(z < z_low) {
				z_low = z;
				z_low_idx = np; // store the index of the point where the altitude is lowest, used to compute the tangent point later
			}
#ifdef GPUDEBUG
			los[np].ir = ir; los[np].ip = np; // for DEBUGging
#endif

			if(stop) { *tsurf = (stop == 2 ? t : -999); break; }											// Hit ground or space?

			double n = 1., ng[] = {0., 0., 0.};
			if(ctl->refrac && z <= zrefrac) {																					// Compute gradient of refractivity
				n += refractivity(p, t);
				double xh[3];
				UNROLL
					for(int i = 0; i < 3; i++) xh[i] = x[i] + 0.5*ds*ex0[i];
				cart2geo(xh, &z, &lon, &lat);
				intpol_atm_geo_pt(ctl, atm, (int) atmIdx, atmNp, z, lon, lat, &p, &t);
				double const n2 = refractivity(p, t);
				for(int i = 0; i < 3; i++) {
					double const h = 0.02;
					xh[i] += h;
					cart2geo(xh, &z, &lon, &lat);
					intpol_atm_geo_pt(ctl, atm, (int) atmIdx, atmNp, z, lon, lat, &p, &t);
					ng[i] = (refractivity(p, t) - n2)/h;
					xh[i] -= h;
				} // i
			}
			UNROLL
				for(int i = 0; i < 3; i++) ex1[i] = ex0[i]*n + ds*ng[i];								// Construct new tangent vector

			double const norm_ex1 = NORM(ex1);
			for(int i = 0; i < 3; i++) {
				ex1[i] /= norm_ex1;
				x[i] += 0.5*ds*(ex0[i] + ex1[i]);																				// Determine next point of LOS
				ex0[i] = ex1[i];																												// Copy tangent vector
			} // i
		} // np
		++np;
#ifndef __NVCC__
        if (NLOS <= np) ERRMSG("Too many LOS points!");
#endif

		// Get tangent point (before changing segment lengths!)
		tangent_point(los, np, z_low_idx, &obs->tpz[ir], &obs->tplon[ir], &obs->tplat[ir]);
		trapezoid_rule_pos(np, los);
		column_density(ctl->ng, los, np);
#ifdef CURTIS_GODSON
		if(ctl->formod == 1) curtis_godson(ctl, los, np);
		// this could be done during the while loop using the aux-variables:
		//					double cgpxu[NG];	/*! Curtis-Godson pressure times column density */
		//					double cgtxu[NG];	/*! Curtis-Godson temperature times column density */
#else
		assert(1 != ctl->formod);
#endif

		return np;
	} // traceray

    // Find air parcel next to reference height
    __host__ __device__ __ext_inline__
    int find_reference_parcel(ctl_t const *ctl, atm_t const *atm, int const ip0, int const ip1) {
        double dzmin = 1e99;
        int ipref = 0;
        for(int ip = ip0; ip < ip1; ip++) {
            double const dz = fabs(atm->z[ip] - ctl->hydz);
            if(dz < dzmin) {
                dzmin = dz;
                ipref = ip;
            }
        } // ip
        return ipref;
    } // find_reference_parcel

    __host__ __device__ __ext_inline__
    void hydrostatic_1d_h2o(ctl_t const *ctl, atm_t *atm, int const ip0, int const ip1, int const ig_h2o) {
        int const npts = 20;
        int const ipref = find_reference_parcel(ctl, atm, ip0, ip1);
        double const lat = atm->lat[ipref];
        double const mmair = 28.96456e-3, mmh2o = 18.0153e-3;
        double e = 0.;

        // Upper part of profile
        for(int ip = ipref + 1; ip < ip1; ip++) {
            double mean = 0.;
            for(int i = 0; i < npts; i++) {
                double const z = lip(0.0, atm->z[ip - 1], npts - 1.0, atm->z[ip], (double) i);
                double const grav = gravity(z, lat);
                if(ig_h2o >= 0) e = lip(0.0, atm->q[ig_h2o][ip - 1], npts - 1.0, atm->q[ig_h2o][ip], (double) i);
                double const temp = lip(0.0, atm->t[ip - 1], npts - 1.0, atm->t[ip], (double) i);
                mean += (e*mmh2o + (1 - e)*mmair)*grav/(GSL_CONST_MKSA_MOLAR_GAS*temp*npts);
            }
            atm->p[ip] = atm->p[ip - 1]*exp(-1000*mean*(atm->z[ip] - atm->z[ip - 1])); // Compute p(z,T)
        } // ip

        // Lower part of profile
        for(int ip = ipref - 1; ip >= ip0; ip--) {
            double mean = 0.;
            for(int i = 0; i < npts; i++) {
                double const z = lip(0.0, atm->z[ip + 1], npts - 1.0, atm->z[ip], (double) i);
                double const grav = gravity(z, lat);
                if(ig_h2o >= 0) e = lip(0.0, atm->q[ig_h2o][ip + 1], npts - 1.0, atm->q[ig_h2o][ip], (double) i);
                double const temp = lip(0.0, atm->t[ip + 1], npts - 1.0, atm->t[ip], (double) i);
                mean += (e*mmh2o + (1 - e)*mmair)*grav/(GSL_CONST_MKSA_MOLAR_GAS*temp*npts);
            } // i
            atm->p[ip] = atm->p[ip + 1]*exp(-1000*mean*(atm->z[ip] - atm->z[ip + 1])); // Compute p(z,T)
        } // ip
    } // hydrostatic_1d
