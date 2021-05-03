/*
	 This file is part of JURASSIC.

	 JURASSIC is free software: you can redistribute it and/or modify
	 it under the terms of the GNU General Public License as published by
	 the Free Software Foundation, either version 3 of the License, or
	 (at your option) any later version.

	 JURASSIC is distributed in the hope that it will be useful,
	 but WITHOUT ANY WARRANTY; without even the implied warranty of
	 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	 GNU General Public License for more details.

	 You should have received a copy of the GNU General Public License
	 along with JURASSIC. If not, see <http://www.gnu.org/licenses/>.

	 Copright (C) 2003-2015 Forschungszentrum Juelich GmbH
	 */

/*!
	\file
	JURASSIC library definitions.
	*/
#include "jurassic.h"
#include "jr_binary_tables_io.h" // jr_write_binary_tables, jr_read_binary_tables

    double gravity(double const z, double const lat);
//     {
//         double const deg2rad = M_PI/180., x = sin(lat*deg2rad), y = sin(2*lat*deg2rad);
//         return 9.780318*(1. + 0.0053024*x*x - 5.8e-6*y*y) - 3.086e-3*z;
//     } // gravity

	int locate(double const *xx, int const n, double const x); 
//     {
//         int ilo = 0, ihi = n - 1, i = (n - 1) >> 1;
//         if(xx[i] < xx[i + 1]) {
//             while (ihi > ilo + 1) { // divergent execution on GPU happens here
//                 i = (ihi + ilo) >> 1;
//                 if(xx[i] > x) ihi = i;
//                 else					ilo = i;
//             } // while
//         } else {
//             while (ihi > ilo + 1) {
//                 i = (ihi + ilo) >> 1;
//                 if(xx[i] <= x) ihi = i;
//                 else					 ilo = i;
//             } // while
//         } // if
//         return ilo;
//     } // locate

// #define grd2rad (M_PI/180)
	void geo2cart(double const alt, double const lon, double const lat, double x[]);
//     {
// 		double const radius = alt + RE, clat = cos(lat*grd2rad);
// 		x[0] = radius*clat*cos(lon*grd2rad);
// 		x[1] = radius*clat*sin(lon*grd2rad);
// 		x[2] = radius*sin(lat*grd2rad);
// 	} // geo2cart
    

#define rad2grd (180/M_PI)
#define grd2rad (M_PI/180)


// /////////////////////////////////////////////////////////////////////////
FILE* mkFile(char const* dn, char const* fn, char const* acc) {
	FILE* out;
	char file[LEN];
	if(dn) sprintf(file, "%s/%s", dn, fn);
	else	 sprintf(file, "%s", fn);
	if(!(out = fopen(file, acc))) ERRMSG(file);
	return out;
}

void shell(char const* cmd, char const* err) { if(system(cmd)) ERRMSG(err); }

//***************************************************************************
void climatology(ctl_t const *ctl, atm_t *atm) {
    #include "climatology.tbl"
	static int ig_co2 = -999;
	double co2;
    double const *q[NG] = { NULL };
	// Find emitter index of CO2
	if (ig_co2 == -999) ig_co2 = find_emitter(ctl, "CO2");
	// Identify variable
	for(int ig = 0; ig < ctl->ng; ig++) {
		q[ig] = NULL;
        char const *gasname = ctl->emitter[ig];
        if (0) { ; }
		else if(0 == strcasecmp(gasname, "C2H2")  ) q[ig] = c2h2;
		else if(0 == strcasecmp(gasname, "C2H6")  ) q[ig] = c2h6;
		else if(0 == strcasecmp(gasname, "CCl4")  ) q[ig] = ccl4;
		else if(0 == strcasecmp(gasname, "CH4")   ) q[ig] = ch4;
		else if(0 == strcasecmp(gasname, "ClO")   ) q[ig] = clo;
		else if(0 == strcasecmp(gasname, "ClONO2")) q[ig] = clono2;
		else if(0 == strcasecmp(gasname, "CO")    ) q[ig] = co;
		else if(0 == strcasecmp(gasname, "COF2")  ) q[ig] = cof2;
		else if(0 == strcasecmp(gasname, "F11")   ) q[ig] = f11;
		else if(0 == strcasecmp(gasname, "F12")   ) q[ig] = f12;
		else if(0 == strcasecmp(gasname, "F14")   ) q[ig] = f14;
		else if(0 == strcasecmp(gasname, "F22")   ) q[ig] = f22;
		else if(0 == strcasecmp(gasname, "H2O")   ) q[ig] = h2o;
		else if(0 == strcasecmp(gasname, "H2O2")  ) q[ig] = h2o2;
		else if(0 == strcasecmp(gasname, "HCN")   ) q[ig] = hcn;
		else if(0 == strcasecmp(gasname, "HNO3")  ) q[ig] = hno3;
		else if(0 == strcasecmp(gasname, "HNO4")  ) q[ig] = hno4;
		else if(0 == strcasecmp(gasname, "HOCl")  ) q[ig] = hocl;
		else if(0 == strcasecmp(gasname, "N2O")   ) q[ig] = n2o;
		else if(0 == strcasecmp(gasname, "N2O5")  ) q[ig] = n2o5;
		else if(0 == strcasecmp(gasname, "NH3")   ) q[ig] = nh3;
		else if(0 == strcasecmp(gasname, "NO")    ) q[ig] = no;
		else if(0 == strcasecmp(gasname, "NO2")   ) q[ig] = no2;
		else if(0 == strcasecmp(gasname, "O3")    ) q[ig] = o3;
		else if(0 == strcasecmp(gasname, "OCS")   ) q[ig] = ocs;
		else if(0 == strcasecmp(gasname, "SF6")   ) q[ig] = sf6;
		else if(0 == strcasecmp(gasname, "SO2")   ) q[ig] = so2;
        else { 
              printf("# Warning! no climatology table for found emitter %s\n", gasname);
        }
        
        if (ctl->checkmode) printf("# found emitter %s at gas index %d\n", gasname, ig);
	}
	
	if (ctl->checkmode) return; // early return
	
	for(int ip = 0; ip < atm->np; ip++) {						// Loop over atmospheric data points
		int iz = locate(z, 121, atm->z[ip]);					// Get altitude index
		atm->p[ip] = EXP(z[iz], pre[iz], z[iz + 1], pre[iz + 1], atm->z[ip]);	// Interpolate pressure
		atm->t[ip] = LIN(z[iz], tem[iz], z[iz + 1], tem[iz + 1], atm->z[ip]);	// Interpolate temperature
		for(int ig = 0; ig < ctl->ng; ig++) {					// Interpolate trace gases
			atm->q[ig][ip] = q[ig] ? LIN(z[iz], q[ig][iz], z[iz + 1], q[ig][iz + 1], atm->z[ip]) : 0;
		}
		if(ig_co2 >= 0) {								// Set CO2
			co2 = 371.789948e-6 + 2.026214e-6 *(atm->time[ip] - 63158400.)/31557600.;
			atm->q[ig_co2][ip] = co2;
		}
		for(int iw = 0; iw < ctl->nw; iw++) atm->k[iw][ip] = 0;			// Set extinction to zero
	}
}


//***************************************************************************
void copy_atm(ctl_t const *const ctl, atm_t * atm_dest, atm_t const *const atm_src, int const init) {
	size_t s = (size_t) atm_src->np*sizeof(double);		// Data size
	// Copy data
	atm_dest->np = atm_src->np;
	memcpy(atm_dest->time, atm_src->time,  s);
	memcpy(atm_dest->z,		 atm_src->z,		 s);
	memcpy(atm_dest->lon,  atm_src->lon,	 s);
	memcpy(atm_dest->lat,  atm_src->lat,	 s);
	memcpy(atm_dest->p,		 atm_src->p,		 s);
	memcpy(atm_dest->t,		 atm_src->t,		 s);
	for(int ig = 0; ig < ctl->ng; ig++) memcpy(atm_dest->q[ig], atm_src->q[ig], s);
	for(int iw = 0; iw < ctl->nw; iw++) memcpy(atm_dest->k[iw], atm_src->k[iw], s);
	atm_dest->init = atm_src->init;
	if(init) {	 // Initialize
		for(int ip = 0; ip < atm_dest->np; ip++) {
			atm_dest->p[ip] = 0;
			atm_dest->t[ip] = 0;
			for(int ig = 0; ig < ctl->ng; ig++) atm_dest->q[ig][ip] = 0;
			for(int iw = 0; iw < ctl->nw; iw++) atm_dest->k[iw][ip] = 0;
		}
	}
}

//***************************************************************************
void copy_obs(ctl_t const *const ctl, obs_t *obs_dest, obs_t const *const obs_src, int const init) {
	size_t s = (size_t) obs_src->nr *sizeof(double);
	// Copy data
	obs_dest->nr = obs_src->nr;
	memcpy(obs_dest->time,	obs_src->time,					 s);
	memcpy(obs_dest->obsz,	obs_src->obsz,					 s);
	memcpy(obs_dest->obslon,	obs_src->obslon,				 s);
	memcpy(obs_dest->obslat,	obs_src->obslat,				 s);
	memcpy(obs_dest->vpz,		obs_src->vpz,						 s);
	memcpy(obs_dest->vplon,	obs_src->vplon,					 s);
	memcpy(obs_dest->vplat,	obs_src->vplat,					 s);
	memcpy(obs_dest->tpz,		obs_src->tpz,						 s);
	memcpy(obs_dest->tplon,	obs_src->tplon,					 s);
	memcpy(obs_dest->tplat,	obs_src->tplat,					 s);
	size_t l = (size_t) ctl->nd *sizeof(double);
	for(int ir = 0; ir < obs_src->nr; ir++) memcpy(obs_dest->rad[ir], obs_src->rad[ir], l);
	for(int ir = 0; ir < obs_src->nr; ir++) memcpy(obs_dest->tau[ir], obs_src->tau[ir], l);
	if(init) {	 // Initialize
		for(int ir = 0; ir < obs_dest->nr; ir++) {
			for(int id = 0; id < ctl->nd; id++) {
				if(gsl_finite(obs_dest->rad[ir][id])) {
					obs_dest->rad[ir][id] = 0;
					obs_dest->tau[ir][id] = 0;
				}
			}
		}
	}
}

//***************************************************************************
int find_emitter(ctl_t const *ctl, char const *emitter) {
	for(int ig = 0; ig < ctl->ng; ig++) {
		if(0 == strcasecmp(ctl->emitter[ig], emitter)) {
            if (ctl->checkmode) printf("# find_emitter %s at gas index %d\n", emitter, ig);
            return ig;
        }
	}
    if (ctl->checkmode) printf("# Warning! did not find_emitter for %s\n", emitter);
	return -1;
}

//***************************************************************************
inline double brightness(double rad, double nu) { return C2*nu/gsl_log1p(C1*gsl_pow_3(nu)/rad); }


//***************************************************************************
void formod_fov(ctl_t const *const ctl, obs_t *obs) {
	static double dz[NSHAPE], w[NSHAPE];
	static int init = 0, n;
	obs_t obs2;
	double rad[NR][ND], tau[NR][ND], z[NR], zfov;
	if(ctl->fov[0] == '-') return;							// Do not take into account FOV
	if(!init) {										// Initialize FOV data
		init = 1;
		n = read_shape(ctl->fov, dz, w, ctl->checkmode);
	}
	copy_obs(ctl, &obs2, obs, 0);								// Copy observation data
	for(int ir = 0; ir < obs->nr; ir++) {							// Loop over ray paths
		int nz = 0;
		for(int ir2 = GSL_MAX(ir - NFOV, 0); ir2 < GSL_MIN(ir + 1 + NFOV, obs->nr); ir2++)	// Get radiance and transmittance profiles
			if(obs->time[ir2] == obs->time[ir]) {
				z[nz] = obs2.vpz[ir2];
				for(int id = 0; id < ctl->nd; id++) {
					rad[nz][id] = obs2.rad[ir2][id];
					tau[nz][id] = obs2.tau[ir2][id];
				}
				nz++;
			}
		if(nz < 2) ERRMSG("Cannot apply FOV convolution!");
		// Convolute profiles with FOV
		double wsum = 0;
		for(int id = 0; id < ctl->nd; id++) {
			obs->rad[ir][id] = 0;
			obs->tau[ir][id] = 0;
		}
		for(int i = 0; i < n; i++) {
			zfov = obs->vpz[ir] + dz[i];
			const int idx = locate(z, nz, zfov);
			for(int id = 0; id < ctl->nd; id++) {
				obs->rad[ir][id] += w[i]*LIN(z[idx], rad[idx][id], z[idx + 1], rad[idx + 1][id], zfov);
				obs->tau[ir][id] += w[i]*LIN(z[idx], tau[idx][id], z[idx + 1], tau[idx + 1][id], zfov);
			}
			wsum += w[i];
		}
		// normalize (could maybe be done in advance by normalizing the weights w[i] beforehand)
		for(int id = 0; id < ctl->nd; id++) {
			obs->rad[ir][id] /= wsum;
			obs->tau[ir][id] /= wsum;
		}
	}
}



//***************************************************************************
void hydrostatic(ctl_t const *const ctl, atm_t *atm) {
	if(ctl->hydz < 0) return; // Check reference height
    if (ctl->checkmode) { printf("# apply hydrostatic equation to individual profiles\n"); return; }
	double lat0 = -999, lon0 = -999;
	int ip0 = -999;
	for(int ip = 0; ip < atm->np; ip++)		// Apply hydrostatic equation to individual profiles
		if((atm->lon[ip] != lon0) || (atm->lat[ip] != lat0)) {
			if(ip > 0) hydrostatic_1d(ctl, atm, ip0, ip);
			lon0 = atm->lon[ip];
			lat0 = atm->lat[ip];
			ip0 = ip;
		}
	hydrostatic_1d(ctl, atm, ip0, atm->np);
}

//***************************************************************************
void hydrostatic_1d(ctl_t const *const ctl, atm_t *atm, int const ip0, int const ip1) {
	static int ig_h2o = -999;
	double dzmin = 1e99, e = 0, mmair = 28.96456e-3, mmh2o = 18.0153e-3;
	int ipref = 0, ipts = 20;
	if(ig_h2o == -999) ig_h2o = find_emitter(ctl, "H2O");					// Determine emitter index of H2O
	for(int ip = ip0; ip < ip1; ip++)							// Find air parcel next to reference height
		if(fabs(atm->z[ip] - ctl->hydz) < dzmin) {
			dzmin = fabs(atm->z[ip] - ctl->hydz);
			ipref = ip;
		}
	for(int ip = ipref + 1; ip < ip1; ip++) {						// Upper part of profile
		double mean = 0;
		for(int i = 0; i < ipts; i++) {
			const double z = LIN(0.0, atm->z[ip - 1], ipts - 1.0, atm->z[ip], (double) i);
			if(ig_h2o >= 0) e = LIN(0.0, atm->q[ig_h2o][ip - 1], ipts - 1.0, atm->q[ig_h2o][ip], (double) i);
			mean += (e*mmh2o + (1 - e)*mmair)*gravity(z, atm->lat[ipref])/GSL_CONST_MKSA_MOLAR_GAS/LIN(0.0, atm->t[ip - 1], ipts - 1.0, atm->t[ip], (double) i)/ipts;
		}
		atm->p[ip] = exp(log(atm->p[ip - 1]) - mean*1000*(atm->z[ip] - atm->z[ip - 1]));		// Compute p(z,T)
	}
	for(int ip = ipref - 1; ip >= ip0; ip--) {						// Lower part of profile
		double mean = 0;
		for(int i = 0; i < ipts; i++) {
			const double z = LIN(0.0, atm->z[ip + 1], ipts - 1.0, atm->z[ip], (double) i);
			if(ig_h2o >= 0) e = LIN(0.0, atm->q[ig_h2o][ip + 1], ipts - 1.0, atm->q[ig_h2o][ip], (double) i);
			mean += (e*mmh2o + (1 - e)*mmair)*gravity(z, atm->lat[ipref])/GSL_CONST_MKSA_MOLAR_GAS/LIN(0.0, atm->t[ip + 1], ipts - 1.0, atm->t[ip], (double) i)/ipts;
		}
		atm->p[ip] = exp(log(atm->p[ip + 1]) - mean*1000*(atm->z[ip] - atm->z[ip + 1]));		// Compute p(z,T)
	}
}


//***************************************************************************
void init_tbl(ctl_t const *ctl, tbl_t *tbl) {
    if (ctl->read_binary) {
        int const binary_reading_status = jr_read_binary_tables(tbl, ctl);
        if (0 == binary_reading_status) {
            printf("matching binary tables file found\n");
            return;
        } else if (ctl->read_binary > 0) {
            ERRMSG("Failed to read binary file while READ_BINARY > 0");
        }
    }
  
    TIMER("INIT_TBL", 1);
    int const logg = 0;
    size_t warn_nu_ignored = 0;
    int max_nu_in_file = 0, warn_missing_file = 0;
	for(int ig = 0; ig < ctl->ng; ig++) { // Loop over trace gases and channels
      int warn_missing_file_gas = 0;
      if (0 == ctl->checkmode) {
#pragma omp parallel for
		for(int id = 0; id < ctl->nd; id++) {
          
			char filename[LEN];
#ifdef      DIRECTORY_WITH_GAS_NAME
			sprintf(filename, "%s/%s_%s/boxcar_%.4f_%s.tab", ctl->tblbase, 
                    ctl->tblbase, ctl->emitter[ig], ctl->nu[id], ctl->emitter[ig]);
#else            
			sprintf(filename, "%s_%.4f_%s.tab", ctl->tblbase, ctl->nu[id], ctl->emitter[ig]);
#endif
			// Try to open file
			FILE *in;
			if(!(in = fopen(filename, "r"))) {
				if(logg) printf("Missing emissivity table: %s\n", filename);
                ++warn_missing_file; ++warn_missing_file_gas; // summarize
				continue;
			}
#define IDX_P (tbl->np[ig][id])
#define IDX_T (tbl->nt[ig][IDX_P][id])
#define IDX_U (tbl->nu[ig][IDX_P][IDX_T][id])
			// Initialize
			IDX_P = -1;
			if(logg) printf("Read emissivity table: %s\n", filename);
            char line[LEN];
			double eps_old = -999, press_old = -999, temp_old = -999, u_old = -999;
            int iu_max = 0;
			while (fgets(line, LEN, in)) {							// Read data
                double eps = 0, press = 0, temp = 0, u = 0;
				if(sscanf(line, "%lg %lg %lg %lg", &press, &temp, &u, &eps) != 4) continue;	// Parse line
				if(press != press_old) {							// Determine pressure index
					press_old = press;
					if((++IDX_P) >= TBLNP) ERRMSG("Too many pressure levels! max. " xstr(TBLNP) " defined in jurassic.h");
					IDX_T = -1;
				}
				if(temp != temp_old) {								// Determine temperature index
					temp_old = temp;
					if((++IDX_T) >= TBLNT) ERRMSG("Too many temperatures! max. " xstr(TBLNT) " defined in jurassic.h");
					IDX_U = -1;
                    iu_max = 0;
				}
				if((eps > eps_old && u > u_old) || IDX_U < 0) {     // Determine column density index
					eps_old = eps;
					u_old = u;
                    ++iu_max;
					if((++IDX_U) >= TBLNU) { 
                        ++warn_nu_ignored; // warning will be shown later
                        max_nu_in_file = GSL_MAX(max_nu_in_file, iu_max);
                        IDX_U--; 
                        continue;
                    } 
				}
				// Store data
				tbl->p[ig][IDX_P][id]                  = press;
				tbl->t[ig][IDX_P][IDX_T][id]           = temp;
				tbl->u[ig][IDX_P][IDX_T][IDX_U][id]    = (real_tblND_t)u;
				tbl->eps[ig][IDX_P][IDX_T][IDX_U][id]  = (real_tblND_t)eps;
//				// Warnings: preparation for the case when p and t become id-independent --> works only without omp loop over id-loop
//				if (press != tbl->p[ig][IDX_P][0]) printf("In emissivity table: %s, pressure differs for channel %d\n", filename, id);
//				if (temp != tbl->t[ig][IDX_P][IDX_T][0]) printf("In emissivity table: %s, temperature differs for channel %d, pressure %f\n", filename, id, press);
			}
			// Increment counters
			IDX_P++;
			for(int ip = 0; ip < IDX_P; ip++) {
				tbl->nt[ig][ip][id]++;
				for(int it = 0; it < tbl->nt[ig][ip][id]; it++) tbl->nu[ig][ip][it][id]++;
			}
			fclose(in);
#undef IDX_U
#undef IDX_T
#undef IDX_P

		} // id
      } else { // checkmode
			char filenames[LEN];
#ifdef      DIRECTORY_WITH_GAS_NAME
			sprintf(filenames, "%s/%s_%s/boxcar_%s_%s.tab", ctl->tblbase, 
                    ctl->tblbase, ctl->emitter[ig], "<nu.4>", ctl->emitter[ig]);
#else            
			sprintf(filenames, "%s_%s_%s.tab", ctl->tblbase, "<nu.4>", ctl->emitter[ig]);
#endif
            printf("# try to initialize tables for gas %d %s from filenames %s\n", 
                   ig, ctl->emitter[ig], filenames);
            printf("# jurassic.h: max dim for table[%d g][%d p][%d T][%d u][%d nu]\n",
                                                     NG, TBLNP, TBLNT, TBLNU, ND);
      } // checkmode
      if (warn_missing_file_gas > 0) printf("%d missing emissivity tables for gas %i %s\n", 
                                              warn_missing_file_gas, ig, ctl->emitter[ig]);
	} // ig
    TIMER("INIT_TBL", 3); // stop timer

	if (warn_nu_ignored > 0) {
        fprintf(stderr, "Warning! %ld table entries ignored, increase TBLNU from %d to %d\n", warn_nu_ignored, TBLNU, max_nu_in_file);
        fprintf(stdout, "Warning! %ld table entries ignored, increase TBLNU from %d to %d\n", warn_nu_ignored, TBLNU, max_nu_in_file);
    }

	if (warn_missing_file > 0) {
        fprintf(stdout, "Warning! %d files were not found!\n", warn_missing_file);
        fprintf(stderr, "Warning! %d files were not found!\n", warn_missing_file);
    }
    
	if (0 == ctl->checkmode) { // scope: show how far we can reduce the limits TBLNP, TBLNT and TBLNU in jurassic.h
        int np_max = 0, np_gd[2]   = {0,0};
        int nt_max = 0, nt_gdp[3]  = {0,0,0};
        int nu_max = 0, nu_gdpt[4] = {0,0,0,0};
        size_t mem_used = 0;
        for(int ig = 0; ig < ctl->ng; ig++) {
            for(int id = 0; id < ctl->nd; id++) {
                int const np = tbl->np[ig][id];
                if (np > np_max) { 
                    np_max = np;
                    np_gd[0] = ig; np_gd[1] = id;
                } // max
                for(int ip = 0; ip < np; ip++) {
                    int const nt = tbl->nt[ig][ip][id];
                    if (nt > nt_max) { 
                        nt_max = nt;
                        nt_gdp[0] = ig; nt_gdp[1] = id; nt_gdp[2] = ip;
                    } // max
                    for(int it = 0; it < nt; it++) {
                        int const nu = tbl->nu[ig][ip][it][id];
                        if (nu > nu_max) { 
                            nu_max = nu;
                            nu_gdpt[0] = ig; nu_gdpt[1] = id; nu_gdpt[2] = ip; nu_gdpt[3] = it;
                        } // max
                        mem_used += (unsigned)nu;
                    } // it
                } // ip
            } // ig
        } // id

        printf("\n# jurassic.h could be configured minimally with\n");
        printf("# NG = %d  \t now %d\n", ctl->ng, NG);
        printf("# ND = %d  \t now %d\n", ctl->nd, ND);
        {   
            int const ig = np_gd[0], id = np_gd[1];
            printf("# TBLNP = %d  \t now %d \t(gas[%d]=%s  nu[%d]=%.4f)\n", 
                   np_max, TBLNP,  ig, ctl->emitter[ig],  id, ctl->nu[id]);
        }
        {   
            int const ig = nt_gdp[0], id = nt_gdp[1], ip = nt_gdp[2];
            printf("# TBLNT = %d  \t now %d \t(gas[%d]=%s  nu[%d]=%.4f  pressure[%d]=%.2e)\n",
                    nt_max, TBLNT,  ig, ctl->emitter[ig],  id, ctl->nu[id],  ip, tbl->p[ig][ip][id]);
        }
        {
            int const ig = nt_gdp[0], id = nt_gdp[1], ip = nt_gdp[2], it = nu_gdpt[3];
            printf("# TBLNU = %d  \t now %d \t(gas[%d]=%s  nu[%d]=%.4f  pressure[%d]=%.2e  temperature[%d]=%g)\n",
                    nu_max, TBLNU,  ig, ctl->emitter[ig],  id, ctl->nu[id],  ip, tbl->p[ig][ip][id],  it, tbl->t[ig][ip][it][id]);
        }
        double const f = 1e-9*sizeof(real_tblND_t)*2; // two arrays: u/eps[NG][TBLNP][TBLNT][TBLNU][ND]
        double const total_now = NG*TBLNP*f*TBLNT*TBLNU*ND;
        printf("# now main table arrays (tbl->u + tbl->eps) consume %.6f GByte\n", total_now);
        double const total_min = (ctl->ng)*np_max*f*nt_max*nu_max*(ctl->nd);
        printf("# main table arrays could consume %.6f GByte\n", total_min);
        double const sparse_mem = f*mem_used;
        printf("# with sparse storage only %.6f GByte (%.1f %%)\n\n", sparse_mem, 100*sparse_mem/total_min);
    } // scope
        

#ifdef  FAST_INVERSE_OF_U
    if (ctl->checkmode) {
        printf("# prepare the fast inverse tbl->u0inv[ig][ip][it][id] due to -D FAST_INVERSE_OF_U\n");
    } else { // scope: 

        int echo = 5; // output level, 0:silent, 1:minimum summary, 2:more detials, ..., 9:all summary levels, >9:every value
        assert(8 == sizeof(unsigned long long));
        unsigned long long ndev = 0, ndenom = 0; // number of integer deviations
        
        int i5dgptu[5] = {-1, -1, -1, -1, -1};
        double maxdev_all = 0, maxdev_val = 0;

        double max_inc = 0, min_inc = 99;
        
        double maxdev_d = 0;
        for(int id = 0; id < ctl->nd; id++) {
            double maxdev_g = 0;
            for(int ig = 0; ig < ctl->ng; ig++) {
    //          printf("In emissivity table: gas %i channel %i #pressures %i\n", ig, id, tbl->np[ig][id]);
                double maxdev_p = 0;
                for(int ip = 0; ip < tbl->np[ig][id]; ip++) {
    //              printf("In emissivity table: gas %i channel %i pressure %i #temperatures %i\n", ig, id, ip, tbl->nt[ig][ip][id]);
                    double maxdev_t = 0;
                    for(int it = 0; it < tbl->nt[ig][ip][id]; it++) {
                        int const nu = tbl->nu[ig][ip][it][id];
    // 					printf("In emissivity table: gas %i pressure %i temperature %i channel %i #u %i\n", ig, ip, it, id, nu);
                      
//                      double const u0 = tbl->u[ig][ip][it][0][id];
//                      double const uM = sqrt(u0 * tbl->u[ig][ip][it][nu - 1][id]); // geometric mean between lowest and highest 
//                         int const iuM = nu/2;
                        
                        double const increment_factor = pow( tbl->u[ig][ip][it][nu - 1][id]
                                                           / tbl->u[ig][ip][it][  0   ][id], 1./(nu - 1.) );
//                      printf("# FAST_INVERSE_OF_U det=%i gas=%i pre=%i tmp=%i increment %.9f\n", id,ig,ip,it, increment_factor);
                        max_inc = GSL_MAX_DBL(max_inc, increment_factor);
                        min_inc = GSL_MIN_DBL(min_inc, increment_factor);
                        
                        int const iuM = nu - 1; // highest!!
                        double const uM = tbl->u[ig][ip][it][iuM][id]; // approximate mid-point (on a log scale)

                        assert(uM > 0);
                        double const u0inv = pow(2.0, iuM/6.)/uM;
                        double const u0 = 1./u0inv; // this u0 does not match u[0]
                        tbl->u0inv[ig][ip][it][id] = u0inv;
                        // now check that the fast inverse of the u-grid works
                        // we assume that the u-grid is defined by u(i) = u0*2^{i/6.020832}
                        double maxdev_u = 0;
                        assert(uM * u0inv < 1.5e51); // since 1.5777e51^6.020832 = 1.79769e308 is the largest representable number in fp64
                        // this is also the reason why fp32 does not work here:
                        //    the largest representable number in fp32 is 3.4e38 
                        //    so covering a log range with 204 increments of 2^(1/6.) gives 2^34
                        //    and 2^(34*6) is 2.57e61, i.e. out of range for fp32
                        for(int iu = 0; iu < nu; ++iu) {
                            double const u_ref = tbl->u[ig][ip][it][iu][id]; // load the value specified in the file
                            double const x = u_ref * u0inv, x2 = x*x, x4 = x2*x2, x6 = x4*x2;

                            unsigned long long const bits = *((unsigned long long*)(&x6)); // reinterpret cast
//                          // now get rid of all the 52 mantissa bits of the fp64-representation of x^6 and the exponent offset 1023
//                          int const iu_fast = (bits >> 52) - 1023; // works but 0.3% index deviations for limb.ctl case
   
                            // now get rid of all but 3 mantissa bits of the fp64-representation of x^6 
                            // and take care of the exponent offset 1023
                            int const approximate_log2_of_x6_times8 = (bits >> (52 - 3)) - (1023 << 3);
                            int const iu_fast = (int)(1.003472 * 0.125 * approximate_log2_of_x6_times8); // conversion to int floors
                            // the correction factor  1.003472 comes from the increment of 1.122 = 2^0.16609641 
                            // compared to 2^(1/6.)
                            // the factor 0.125 makes up for the 3 bits leftover from the mantissa

                            ndev += (iu_fast != iu); ++ndenom;
                            if (iu_fast != iu) printf("# FAST_INVERSE_OF_U indx deviation at det=%i gas=%i pre=%i tmp=%i u %i %i\n",
                                                        id,ig,ip,it,iu,iu_fast);

                            double const u_new = tbl->u[ig][ip][it][0][id] * pow(2.0, iu/(6. * 1.003472));
                            double const reldev = (u_new/u_ref - 1); // relative deviation
                            double const absreldev = fabs(reldev); // absolute of relative deviation
                            if (echo > 9) printf("# FAST_INVERSE_OF_U channel %i gas %i pressure %i temperature %i u[%i] = %g and u[%i] = %g relative dev. %.3e\n", 
                                                                              id,    ig,         ip,            it,  iu,   u_ref,   iu_fast, u_new,         reldev);
                            maxdev_u = GSL_MAX_DBL(maxdev_u, absreldev);
                            if (absreldev > maxdev_all) {
                                i5dgptu[0] = id;
                                i5dgptu[1] = ig;
                                i5dgptu[2] = ip;
                                i5dgptu[3] = it;
                                i5dgptu[4] = iu;
                                maxdev_all = absreldev;
                                maxdev_val = reldev;
                            } //
                            
                        } // iu
                        
                        if (echo > 8) printf("# FAST_INVERSE_OF_U channel %i gas %i pressure %i temperature %i max relative dev. %.3e\n",
                                                                          id,    ig,         ip,            it,                maxdev_u);
                        if (echo > 9) --echo;
                        maxdev_t = GSL_MAX_DBL(maxdev_t, maxdev_u);
                    } // it
                    if (echo > 6) printf("# FAST_INVERSE_OF_U channel %i gas %i pressure %i max relative dev. %.3e\n",
                                                                      id,    ig,         ip,                maxdev_t);
                    if (echo > 8) --echo;
                    maxdev_p = GSL_MAX_DBL(maxdev_p, maxdev_t);
                } // ip
                if (echo > 4) printf("# FAST_INVERSE_OF_U channel %i gas %i max relative dev. %.3e\n",
                                                                  id,    ig,                maxdev_d);
                if (echo > 6) --echo;
                maxdev_g = GSL_MAX_DBL(maxdev_g, maxdev_p);
            } // ig
            if (echo > 2) printf("# FAST_INVERSE_OF_U channel %i max relative dev. %.3e\n",
                                                              id,                maxdev_g);
            if (echo > 4) --echo;
            maxdev_d = GSL_MAX_DBL(maxdev_d, maxdev_g);
        } // id
        if (echo > 0) printf("# FAST_INVERSE_OF_U total largest relative dev. %.3e\n", maxdev_d);
        if (echo > 0) printf("# FAST_INVERSE_OF_U increment factors vary between %.9f and %.9f\n", min_inc, max_inc);
        // found 1.122015007 and 1.122021491, geometric mean is 1.12201825 == 0.99960462 * 2^(1/6.) == 
        if (echo > 1) {
            printf("# FAST_INVERSE_OF_U total largest absolute of relative deviation is %.3e"
                   "  found at channel %i gas %i pressure %i temperature %i u[%i]\n", maxdev_val,
                       i5dgptu[0], i5dgptu[1], i5dgptu[2], i5dgptu[3], i5dgptu[4]);
        } // echo
        assert(maxdev_all == maxdev_d);
        if (echo > 0) printf("# FAST_INVERSE_OF_U %.2f %% index deviations: %lld of %lld\n", ndev/(.01*ndenom), ndev, ndenom);
        assert(0 == ndev);
    } // end of scope
#endif
	

	printf("Initialize source function table...\n");
	for(int it = 0; it < TBLNS; it++) {																	// Compute source function table
		tbl->st[it] = LIN(0.0, 100, TBLNS - 1.0, 400, (double) it); // Set temperature axis: equidistant steps of 300/(TBLNS - 1) Kelvin	, e.g. 0.25 K if TBLNS == 1201
	}

	
	printf("Count initialized tables...\n");
    {
        printf("# per channel  ");
        int total = 0;
        for(int id = 0; id < ctl->nd; id++) {
            int ntab = 0;
            for(int ig = 0; ig < ctl->ng; ig++) 
                ntab += (tbl->np[ig][id] > 1);
            printf(" %d", ntab);
            total += ntab;
        } // id
        printf("  total=%d (of %d)\n", total, ctl->ng*ctl->nd);
    }
    {
        printf("# per gas      ");
        int total = 0;
        for(int ig = 0; ig < ctl->ng; ig++) {
            int ntab = 0;
            for(int id = 0; id < ctl->nd; id++)
                ntab += (tbl->np[ig][id] > 1);
            printf(" %d", ntab);
            total += ntab;
        } // id
        printf("  total=%d (of %d)\n", total, ctl->ng*ctl->nd);
    }

	
#pragma omp parallel for if(0 == ctl->checkmode)
	for(int id = 0; id < ctl->nd; id++) {
		char filename[LEN];
#ifdef  DIRECTORY_WITH_GAS_NAME
        sprintf(filename, "%s/%s_%s/boxcar_%.4f.filt", ctl->tblbase, ctl->tblbase, "CO2", ctl->nu[id]);
#else        
	sprintf(filename, "%s_%.4f.filt", ctl->tblbase, ctl->nu[id]);				// Read filter function
#endif
		double f[NSHAPE], nu[NSHAPE];																 // Arguments for read_shape
		int const n = read_shape(filename, nu, f, ctl->checkmode);
		if (n > NSHAPE) ERRMSG("Increase NSHAPE defined as " xstr(NSHAPE) " in jurassic.h");
      if (0 == ctl->checkmode) {
		for(int it = 0; it < TBLNS; it++) {																	// Compute source function table
			// Integrate Planck function
			double fsum = 0, fpsum = 0;
			for(int i = 0; i < n; i++) {
				fsum  += f[i];
				fpsum += f[i]*planck(tbl->st[it], nu[i]);
			}
			tbl->sr[it][id] = fpsum / fsum;
		}
      } // checkmode
	}

	if (ctl->write_binary) {
        jr_write_binary_tables(tbl, ctl);
    }
}

//***************************************************************************
void intpol_atm(ctl_t *ctl, atm_t *atm_dest, atm_t *atm_src) {
	for(int ip = 0; ip < atm_dest->np; ip++) {	 // Interpolate atmospheric data
		double k[NW], q[NG];
		intpol_atm_geo(ctl, atm_src, atm_dest->z[ip], atm_dest->lon[ip], atm_dest->lat[ip], &atm_dest->p[ip], &atm_dest->t[ip], q, k);
		for(int ig = 0; ig < ctl->ng; ig++) atm_dest->q[ig][ip] = q[ig];
		for(int iw = 0; iw < ctl->nw; iw++) atm_dest->k[iw][ip] = k[iw];
	}
}

//***************************************************************************
void intpol_atm_geo(ctl_t const *const ctl, atm_t *atm, double const z0, double const lon0, double const lat0,
		double *p, double *t, double *q, double *k) {
	if		 (ctl->ip == 1) intpol_atm_1d(ctl, atm, 0, atm->np, z0, p, t, q, k); // 1D interpolation (vertical profile)
	else if(ctl->ip == 2) intpol_atm_2d(ctl, atm, z0, lon0, lat0, p, t, q, k); // 2D interpolation (satellite track)
	else if(ctl->ip == 3) intpol_atm_3d(ctl, atm, z0, lon0, lat0, p, t, q, k); // 3D interpolation (Lagrangian grid)
	else									ERRMSG("Unknown interpolation method, check IP!");	 // Wrong parameter
}

//***************************************************************************
void intpol_atm_1d(ctl_t const *const ctl, atm_t const *const atm, int const idx0, int const n, double const z0,
		double *p, double *t, double *q, double *k) {
	const int ip = idx0 + locate(&atm->z[idx0], n, z0);												 // Get array index
	*p = EXP(atm->z[ip], atm->p[ip], atm->z[ip + 1], atm->p[ip + 1], z0);			 // Interpolate
	*t = LIN(atm->z[ip], atm->t[ip], atm->z[ip + 1], atm->t[ip + 1], z0);
	for(int ig = 0; ig < ctl->ng; ig++) q[ig] = LIN(atm->z[ip], atm->q[ig][ip], atm->z[ip + 1], atm->q[ig][ip + 1], z0);
	for(int iw = 0; iw < ctl->nw; iw++) k[iw] = LIN(atm->z[ip], atm->k[iw][ip], atm->z[ip + 1], atm->k[iw][ip + 1], z0);
}

//***************************************************************************
void intpol_atm_2d(ctl_t const *const ctl, atm_t * atm,
		double const z0, double const lon0, double const lat0,
		double *p, double *t, double *q, double *k) {
	static double x1[NP][3];
	static int idx[NP], nx, nz[NP];
	double dhmin0 = 1e99, dhmin1 = 1e99, dlat = 10, k0[NW], k1[NW], lat1 = -999, lon1 = -999, p0, p1, q0[NG], q1[NG], r, t0, t1, x0[3];
	int ix0 = 0, ix1 = 0;
	if(!atm->init) {																													 // Initialize
		atm->init = 1;
		// Determine grid dimensions
		nx = 0;
		for(int ip = 0; ip < atm->np; ip++) {
			if((atm->lon[ip] != lon1) || (atm->lat[ip] != lat1)) {
				if((++nx) > NP) ERRMSG("Too many profiles!");
				nz[nx - 1] = 0;
				lon1 = atm->lon[ip];
				lat1 = atm->lat[ip];
				geo2cart(0, lon1, lat1, x1[nx - 1]);
				idx[nx - 1] = ip;
			}
			++nz[nx - 1];
		}
		for(int ix = 0; ix < nx; ix++) {
			if(nz[ix] <= 1)																													 ERRMSG("Cannot identify profiles. Check ordering of data points!");
			if((ix > 0) && (fabs(atm->lat[idx[ix - 1]] - atm->lat[idx[ix]]) > dlat)) ERRMSG("Distance of profiles is too large!");
		}
	}
	geo2cart(0, lon0, lat0, x0);				// Get Cartesian coordinates
	for(int ix = 0; ix < nx; ix++) {			// Find next neighbours
		if(fabs(lat0 - atm->lat[idx[ix]]) <= dlat) {	// Get squared horizontal distance
			const double dh = DIST2(x0, x1[ix]);
			if(dh <= dhmin0) {				// Find neighbours
				dhmin1 = dhmin0;
				ix1 = ix0;
				dhmin0 = dh;
				ix0 = ix;
			} else if(dh <= dhmin1) {
				dhmin1 = dh;
				ix1 = ix;
			}
		}
	}
	// Interpolate vertically
	intpol_atm_1d(ctl, atm, idx[ix0], nz[ix0], z0, &p0, &t0, q0, k0);
	intpol_atm_1d(ctl, atm, idx[ix1], nz[ix1], z0, &p1, &t1, q1, k1);
	// Interpolate horizontally
	const double x2 = DIST2(x1[ix0], x1[ix1]);
	const double x = sqrt(x2);
	const double r0 = (dhmin0 - dhmin1 + x2)/(2*x);
	const double r1 = x - r0;
	if(r0 <= 0) r = 0;
	else				r = (r1 <= 0) ? 1 : r0/(r0 + r1);
	*p = (1 - r)*p0 + r*p1;
	*t = (1 - r)*t0 + r*t1;
	for(int ig = 0; ig < ctl->ng; ig++) q[ig] = (1 - r)*q0[ig] + r*q1[ig];
	for(int iw = 0; iw < ctl->nw; iw++) k[iw] = (1 - r)*k0[iw] + r*k1[iw];
}

//***************************************************************************
void intpol_atm_3d(ctl_t const *const ctl, atm_t *atm,
		double const z0, double const lon0, double const lat0,
		double *p, double *t, double *q, double *k) {
	static double rm2, x1[NP][3];
	if(!atm->init) {						// Initialize
		atm->init = 1;
		// Get Cartesian coordinates
		for(int ip = 0; ip < atm->np; ip++) geo2cart(0, atm->lon[ip], atm->lat[ip], x1[ip]);
		rm2 = gsl_pow_2(ctl->cx);					// Get squared influence radius
	}
	// Initialize for interpolation
	double wsum = 0, x0[3];
	*p = *t = 0.;
	for(int ig = 0; ig < ctl->ng; ig++) q[ig] = 0;
	for(int iw = 0; iw < ctl->nw; iw++) k[iw] = 0;
	for(int ip = 0; ip < atm->np; ip++) {				// Loop over grid points
		const double dz = fabs(atm->z[ip] - z0);			// Get vertical distance
		if(dz >= ctl->cz) continue;
		if(fabs(atm->lat[ip] - lat0)*111.13 >= ctl->cx) continue;		// Check latitude distance
		geo2cart(0, lon0, lat0, x0);				// Get horizontal distance
		const double dx2 = DIST2(x0, x1[ip]);
		if(dx2 >= rm2) continue;
		const double w = (1 - dz/ctl->cz)*(rm2 - dx2)/(rm2 + dx2);	// Set distance-based weighting factor
		// Average data
		wsum += w;
		*p += w*atm->p[ip];
		*t += w*atm->t[ip];
		for(int ig = 0; ig < ctl->ng; ig++) q[ig] += w*atm->q[ig][ip];
		for(int iw = 0; iw < ctl->nw; iw++) k[iw] += w*atm->k[iw][ip];
	}
	// Check sum of weights
	if(wsum >= 1e-6) {						// Normalize
		*p /= wsum;
		*t /= wsum;
		for(int ig = 0; ig < ctl->ng; ig++) q[ig] /= wsum;
		for(int iw = 0; iw < ctl->nw; iw++) k[iw] /= wsum;
	} else {							// Set to nan
		*p = *t = GSL_NAN;
		for(int ig = 0; ig < ctl->ng; ig++) q[ig] = GSL_NAN;
		for(int iw = 0; iw < ctl->nw; iw++) k[iw] = GSL_NAN;
	}
}


//***************************************************************************
// Refractivity of air at 4 to 15 micron
// inline double refractivity(double p, double t) { return 7.753e-05*p/t; } // moved to jr_common.h

//***************************************************************************
void kernel(ctl_t const *const ctl, atm_t *atm, obs_t *obs, gsl_matrix *k) {
	int iqa[N];
	formod(ctl, atm, obs);													 // Compute radiance for undisturbed atmospheric data
	// Compose vectors
	size_t m = k->size1, n = k->size2;
	//printf("k->size2=%ld\n", (long int) n);
	gsl_vector *x0, *yy0;
	x0	= gsl_vector_alloc(n);
	yy0 = gsl_vector_alloc(m);
	atm2x(ctl, atm, x0, iqa, NULL);
	obs2y(ctl, obs, yy0, NULL, NULL);
	gsl_matrix_set_zero(k);													// Initialize kernel matrix
	//	#pragma omp parallel default(shared)
	{
		gsl_vector *x1 = gsl_vector_alloc(n), *yy1 = gsl_vector_alloc(m);
		atm_t* atm1 = malloc(sizeof(atm_t));
		obs_t* obs1 = malloc(sizeof(obs_t));
		//#pragma omp for
		for(int j = 0; j < (int) n; j++) {							// Loop over state vector elements
			// Set perturbation size
			double h;
			if(iqa[j] == IDXP)																	 h = GSL_MAX(fabs(0.01*gsl_vector_get(x0, (size_t) j)), 1e-7);
			else if(iqa[j] == IDXT)															 h = 1;
			else if(iqa[j] >= IDXQ(0) && iqa[j] < IDXQ(ctl->ng)) h = GSL_MAX(fabs(0.01*gsl_vector_get(x0, (size_t) j)), 1e-15);
			else if(iqa[j] >= IDXK(0) && iqa[j] < IDXK(ctl->nw)) h = 1e-4;
			else																								 ERRMSG("Cannot set perturbation size!");
			// Disturb state vector element
			gsl_vector_memcpy(x1, x0);
			gsl_vector_set(x1, (size_t) j, gsl_vector_get(x1, (size_t) j) + h);
			copy_atm(ctl, atm1, atm, 0);
			copy_obs(ctl, obs1, obs, 0);
			x2atm(ctl, x1, atm1);
			formod(ctl, atm1, obs1);											// Compute radiance for disturbed atmospheric data
			obs2y(ctl, obs1, yy1, NULL, NULL);						// Compose measurement vector for disturbed radiance data
			for(size_t i = 0; i < m; i++) {								// Compute derivatives
				gsl_matrix_set(k, i, (size_t)j, (gsl_vector_get(yy1, i) - gsl_vector_get(yy0, i))/h);
			}
		}
		gsl_vector_free(x1);
		gsl_vector_free(yy1);
		free(obs1);
		free(atm1);
	}
	gsl_vector_free(x0);
	gsl_vector_free(yy0);
}

//***************************************************************************
inline double planck(double t, double nu) { return C1*gsl_pow_3(nu)/gsl_expm1(C2 *nu/t); }

//***************************************************************************
void altitude_range(atm_t const *atm, double *zmin, double *zmax) {
	*zmax = *zmin = atm->z[0];
	for(int ipp = 0;
			(ipp < atm->np) && (atm->lon[ipp] == atm->lon[0]) && (atm->lat[ipp] == atm->lat[0]);
			++ipp) {
		*zmax = fmax(*zmax, atm->z[ipp]);
		*zmin = fmin(*zmin, atm->z[ipp]);
	}
}

// Change segment lengths according to trapezoid rule
void trapezoid_rule(int const np, double ds[]) {
	for(int ip = np - 1; ip >= 1; ip--) {
        ds[ip] = 0.5*(ds[ip - 1] + ds[ip]);
    }
	ds[0] *= 0.5;
}

//***************************************************************************
void read_atm(char const *dirname, char const *filename, ctl_t *ctl, atm_t *atm) {
	FILE *in;
	char line[LEN], *tok, *saveptr;
	int ig, iw;
	// Init
	atm->init = 0;
	atm->np = 0;
	printf("Read atmospheric data: %s/%s\n", dirname, filename);
	// Open file
	in = mkFile(dirname, filename, "r");
    if (ctl->checkmode) { 
        printf("# read_atm can read max %d points\n", NP); 
        printf("# read_atm found file %s/%s but skip\n", dirname, filename); 
        fclose(in); return; 
    } // checkmode
	// Read line
	while (fgets(line, LEN, in)) {
		// Read data
		TOK(line, tok, "%lg", atm->time[atm->np], &saveptr);
		TOK(NULL, tok, "%lg", atm->z[atm->np], &saveptr);
		TOK(NULL, tok, "%lg", atm->lon[atm->np], &saveptr);
		TOK(NULL, tok, "%lg", atm->lat[atm->np], &saveptr);
		TOK(NULL, tok, "%lg", atm->p[atm->np], &saveptr);
		TOK(NULL, tok, "%lg", atm->t[atm->np], &saveptr);
		for(ig = 0; ig < ctl->ng; ig++) TOK(NULL, tok, "%lg", atm->q[ig][atm->np], &saveptr);
		for(iw = 0; iw < ctl->nw; iw++) TOK(NULL, tok, "%lg", atm->k[iw][atm->np], &saveptr);
		// Increment data point counter
		if((++atm->np) > NP) ERRMSG("Too many data points!");
	}
	// Close file
	fclose(in);
	// Check number of points
	if(atm->np < 1) ERRMSG("Could not read any data!");
	printf("Read atmospheric data found %d height levels, max %d\n", atm->np, NP);
}


//***************************************************************************
void read_ctl(int argc, char *argv[], ctl_t *ctl) {
	// Write info
	printf("\nJuelich Rapid Spectral Simulation Code (JURASSIC)\n"
			"(executable: %s | compiled: %s, %s)\n\n", argv[0], __DATE__, __TIME__);
#ifdef SHOW_GIT_KEY
	printf("# JURASSIC git commit " xstr(SHOW_GIT_KEY) "\n\n");
#endif
	// Emitters
	ctl->ng = (int) scan_ctl(argc, argv, "NG", -1, "0", NULL);
	if(ctl->ng < 0 || ctl->ng > NG) ERRMSG("Set 0 <= NG <= " xstr(NG) " (max. defined in jurassic.h)");
	for(int ig = 0; ig < ctl->ng; ig++) {
        scan_ctl(argc, argv, "EMITTER", ig, "", ctl->emitter[ig]);
    } // ig
	// Radiance channels
	ctl->nd = (int) scan_ctl(argc, argv, "ND", -1, "0", NULL);
	if(ctl->nd < 0 || ctl->nd > ND) ERRMSG("Set 0 <= ND <= " xstr(ND) " (max. defined in jurassic.h)");
	for(int id = 0; id < ctl->nd; id++) {
        ctl->nu[id] = scan_ctl(argc, argv, "NU", id, "", NULL);
    } // id
	// Spectral windows
	ctl->nw = (int) scan_ctl(argc, argv, "NW", -1, "1", NULL);
	if(ctl->nw < 0 || ctl->nw > NW) ERRMSG("Set 0 <= NW <= " xstr(NW) " (max. defined in jurassic.h)");
	for(int id = 0; id < ctl->nd; id++) {
        ctl->window[id] = (int) scan_ctl(argc, argv, "WINDOW", id, "0", NULL);
    } // id
	// Emissivity look-up tables
	scan_ctl(argc, argv, "TBLBASE", -1, "-", ctl->tblbase);
	// Hydrostatic equilibrium
	ctl->hydz = scan_ctl(argc, argv, "HYDZ", -1, "-999", NULL);
	// Continua
	ctl->ctm_co2 = (int) scan_ctl(argc, argv, "CTM_CO2", -1, "1", NULL);
	ctl->ctm_h2o = (int) scan_ctl(argc, argv, "CTM_H2O", -1, "1", NULL);
	ctl->ctm_n2 = (int) scan_ctl(argc, argv, "CTM_N2", -1, "1", NULL);
	ctl->ctm_o2 = (int) scan_ctl(argc, argv, "CTM_O2", -1, "1", NULL);
	if (1) {
		// automatic control of gases: CTM_...
		int in_co2 = 0, in_h2o = 0, in_n2 = 0, in_o2 = 0; // counters how many frequencies are in range
		for(int id = 0; id < ctl->nd; id++) {
			double const nu = ctl->nu[id]; // abbreviate
			in_co2 += (nu <  4000); // xw = nu/2	+ 1; if(xw >= 1.0 && xw < 2001.0) non-zero
			in_h2o += (nu < 20000); // xw = nu/10 + 1; if(xw >= 1.0 && xw < 2001.0) non-zero
			in_n2  += (nu >= 2120 && nu <= 2605); // if(nu < 2120 || nu > 2605) return 0;
			in_o2  += (nu >= 1360 && nu <= 1805); // if(nu < 1360 || nu > 1805) return 0;
		}
		if(0 == in_co2 && ctl->ctm_co2) { ctl->ctm_co2 = 0; printf("No frequency in CO2 range, automatically set CTM_CO2 = 0\n"); }
		if(0 == in_h2o && ctl->ctm_h2o) { ctl->ctm_h2o = 0; printf("No frequency in H2O range, automatically set CTM_H20 = 0\n"); }
		if(0 == in_n2  && ctl->ctm_n2)	{ ctl->ctm_n2  = 0; printf("No frequency in N2 range, automatically set CTM_N2 = 0\n"); }
		if(0 == in_o2  && ctl->ctm_o2)	{ ctl->ctm_o2  = 0; printf("No frequency in O2 range, automatically set CTM_O2 = 0\n"); }
	}
	// Interpolation of atmospheric data
	ctl->ip = (int) scan_ctl(argc, argv, "IP", -1, "1", NULL);
	ctl->cz = scan_ctl(argc, argv, "CZ", -1, "0", NULL);
	ctl->cx = scan_ctl(argc, argv, "CX", -1, "0", NULL);
	// Ray-tracing
	ctl->refrac = (int) scan_ctl(argc, argv, "REFRAC", -1, "1", NULL);
	ctl->rayds = scan_ctl(argc, argv, "RAYDS", -1, "10", NULL);
	ctl->raydz = scan_ctl(argc, argv, "RAYDZ", -1, "0.5", NULL);
	// Field of view
	scan_ctl(argc, argv, "FOV", -1, "-", ctl->fov);
	// Retrieval interface
	ctl->retp_zmin = scan_ctl(argc, argv, "RETP_ZMIN", -1, "-999", NULL);
	ctl->retp_zmax = scan_ctl(argc, argv, "RETP_ZMAX", -1, "-999", NULL);
	ctl->rett_zmin = scan_ctl(argc, argv, "RETT_ZMIN", -1, "-999", NULL);
	ctl->rett_zmax = scan_ctl(argc, argv, "RETT_ZMAX", -1, "-999", NULL);
	for(int ig = 0; ig < ctl->ng; ig++) {
		ctl->retq_zmin[ig] = scan_ctl(argc, argv, "RETQ_ZMIN", ig, "-999", NULL);
		ctl->retq_zmax[ig] = scan_ctl(argc, argv, "RETQ_ZMAX", ig, "-999", NULL);
	}
	for(int iw = 0; iw < ctl->nw; iw++) {
		ctl->retk_zmin[iw] = scan_ctl(argc, argv, "RETK_ZMIN", iw, "-999", NULL);
		ctl->retk_zmax[iw] = scan_ctl(argc, argv, "RETK_ZMAX", iw, "-999", NULL);
	}
	// Output flags
	ctl->write_bbt = (int) scan_ctl(argc, argv, "WRITE_BBT", -1, "0", NULL);
	ctl->write_matrix =
		(int) scan_ctl(argc, argv, "WRITE_MATRIX", -1, "0", NULL);
	// External forward models
	ctl->formod = (int) scan_ctl(argc, argv, "FORMOD", -1, "2", NULL);
	scan_ctl(argc, argv, "RFMBIN", -1, "-", ctl->rfmbin);
	scan_ctl(argc, argv, "RFMHIT", -1, "-", ctl->rfmhit);
	for(int ig = 0; ig < ctl->ng; ig++) {
        scan_ctl(argc, argv, "RFMXSC", ig, "-", ctl->rfmxsc[ig]);
    } // ig

	ctl->useGPU = (int) scan_ctl(argc, argv, "USEGPU", -1, "0", NULL);
#ifndef hasGPU  
	if (ctl->useGPU > 0) {
		ERRMSG("Requested USEGPU = 1 (always) but compiled without -D hasGPU");
	} else if (ctl->useGPU < 0) {
		fprintf(stderr, "\n\nRequested USEGPU = %d but compiled without CUDA, default to CPU\n\n", ctl->useGPU);
		ctl->useGPU = 0;
	}
#endif
  	
	ctl->checkmode = (int) scan_ctl(argc, argv, "CHECKMODE", -1, "0", NULL);
    	printf("CHECKMODE = %d (%s)\n", ctl->checkmode, 
           (0 == ctl->checkmode)?"run":((ctl->checkmode > 0)?"skip":"obs"));
        
    ctl->read_binary  = (int) scan_ctl(argc, argv, "READ_BINARY", -1, "-1", NULL);
    ctl->write_binary = (int) scan_ctl(argc, argv, "WRITE_BINARY", -1, "1", NULL);

    ctl->gpu_nbytes_shared_memory = (int) scan_ctl(argc, argv, "GPU_SHARED_MEMORY", -1, "0", NULL);
}

//***************************************************************************
void read_matrix(char const *dirname, char const *filename, gsl_matrix *matrix) {
	char dum[LEN], line[LEN];
	double value;
	int i, j;
	printf("Read matrix: %s/%s\n", dirname, filename);
	FILE* in = mkFile(dirname, filename, "r");
	// Read data
	gsl_matrix_set_zero(matrix);
	while (fgets(line, LEN, in))
		if(sscanf(line, "%d %s %s %s %s %s %d %s %s %s %s %s %lg",
					&i, dum, dum, dum, dum, dum,
					&j, dum, dum, dum, dum, dum, &value) == 13) gsl_matrix_set(matrix, (size_t) i, (size_t) j, value);
	fclose(in);
}

//***************************************************************************
void read_obs(char const *dirname, char const *filename, ctl_t *ctl, obs_t *obs) {
	char line[LEN], *tok, *saveptr;
	obs->nr = 0;						// Init
	printf("Read observation data: %s/%s\n", dirname, filename);
	FILE* in = mkFile(dirname, filename, "r");
    if (ctl->checkmode > 0) { 
        printf("# read_obs can read max %d rays\n", NR); 
        printf("# read_obs found file %s/%s but skip\n", dirname, filename); 
        fclose(in); return;
    } // checkmode
	while (fgets(line, LEN, in)) {			// Read line
		TOK(line, tok, "%lg", obs->time[obs->nr], &saveptr);		// Read data
		TOK(NULL, tok, "%lg", obs->obsz[obs->nr], &saveptr);
		TOK(NULL, tok, "%lg", obs->obslon[obs->nr], &saveptr);
		TOK(NULL, tok, "%lg", obs->obslat[obs->nr], &saveptr);
		TOK(NULL, tok, "%lg", obs->vpz[obs->nr], &saveptr);
		TOK(NULL, tok, "%lg", obs->vplon[obs->nr], &saveptr);
		TOK(NULL, tok, "%lg", obs->vplat[obs->nr], &saveptr);
		TOK(NULL, tok, "%lg", obs->tpz[obs->nr], &saveptr);
		TOK(NULL, tok, "%lg", obs->tplon[obs->nr], &saveptr);
		TOK(NULL, tok, "%lg", obs->tplat[obs->nr], &saveptr);
		for(int id = 0; id < ctl->nd; id++) TOK(NULL, tok, "%lg", obs->rad[obs->nr][id], &saveptr);
		for(int id = 0; id < ctl->nd; id++) TOK(NULL, tok, "%lg", obs->tau[obs->nr][id], &saveptr);
		if((++obs->nr) > NR) ERRMSG("Too many rays!");	// Increment counter
	}
	fclose(in);
	if(obs->nr < 1) ERRMSG("Could not read any data!");		// Check number of points
}

//***************************************************************************
double read_obs_rfm(char const *basename, double z, double *nu, double *f, int n) {
	FILE *in;
	char filename[LEN];
	double filt, fsum = 0, nu2[NSHAPE], *nurfm, *rad, radsum = 0;
	int npts;
	ALLOC(nurfm, double, RFMNPTS);
	ALLOC(rad,	 double, RFMNPTS);
	// Search RFM spectrum
	sprintf(filename, "%s_%05d.asc", basename, (int) (z *1000));
	if(!(in = fopen(filename, "r"))) {
		sprintf(filename, "%s_%05d.asc", basename, (int) (z *1000) + 1);
		in = mkFile(NULL, filename, "r");
	}
	fclose(in);
	// Read RFM spectrum
	read_rfm_spec(filename, nurfm, rad, &npts);
	// Set wavenumbers
	nu2[0] = nu[0];
	for(int i = 1; i < n - 1; i++) nu2[i] = LIN(0.0, nu2[0], n - 1.0, nu2[n - 1], i);
	nu2[n - 1] = nu[n - 1];
	// Convolute
	for(int ipts = 0; ipts < npts; ipts++) {
		if(nurfm[ipts] >= nu2[0] && nurfm[ipts] <= nu2[n - 1]) {
			const int idx = locate(nu2, n, nurfm[ipts]);
			filt = LIN(nu2[idx], f[idx], nu2[idx + 1], f[idx + 1], nurfm[ipts]);
			fsum += filt;
			radsum += filt *rad[ipts];
		}
	}
	free(nurfm);
	free(rad);
	// Return radiance
	return radsum/fsum;
}

//***************************************************************************
void read_rfm_spec(char const *filename, double *nu, double *rad, int *npts) {
	char line[RFMLINE], *tok;
	double dnu, nu0, nu1;
	int ipts = 0;
	printf("Read RFM data: %s\n", filename);
	FILE* in = mkFile(NULL, filename, "r");
	// Read header...
	for(int i = 0; i < 4; i++) {
		if(fgets(line, RFMLINE, in) == NULL) ERRMSG("Error while reading file header!");
	}
	sscanf(line, "%d %lg %lg %lg", npts, &nu0, &dnu, &nu1);
	if(*npts > RFMNPTS) ERRMSG("Too many spectral grid points!");
	// Read radiance data
	while (fgets(line, RFMLINE, in) && ipts < *npts - 1) {
		if((tok = strtok(line, " \t\n")))
			if(sscanf(tok, "%lg", &rad[ipts]) == 1) ipts++;
		while ((tok = strtok(NULL, " \t\n")))
			if(sscanf(tok, "%lg", &rad[ipts]) == 1) ipts++;
	}
	if(ipts != *npts) ERRMSG("Error while reading RFM data!");
	// Compute wavenumbers
	for(ipts = 0; ipts < *npts; ipts++) nu[ipts] = LIN(0.0, nu0, (double) (*npts - 1), nu1, (double) ipts);
	fclose(in);
}

//***************************************************************************

int read_shape(char const *filename, double *x, double *y, int const checkmode) {
	char line[LEN];
	printf("Read shape function: %s\n", filename);
	FILE *in = mkFile(NULL, filename, "r");
    if (checkmode) { fclose(in); printf("# read_shape found %s\n", filename); return 0; }
	// Read data
	int n = 0;
	while (fgets(line, LEN, in)) {
		if(sscanf(line, "%lg %lg", &x[n], &y[n]) == 2) {
			if((++n) > NSHAPE) ERRMSG("Too many data points!");
		}
	}
	// Check number of points
	if(n < 1) ERRMSG("Could not read any data!");
	fclose(in);
	return n;
}

//***************************************************************************
double scan_ctl(int argc, char *argv[], char const *varname, int arridx, char const *defvalue, char *value) {
	FILE *in = NULL;
	char dummy[LEN], fullname1[LEN], fullname2[LEN], line[LEN], msg[LEN], rvarname[LEN], rval[LEN];
	int contain = 0;
	// Open file
	if((argv[1][0] != '-')) in = mkFile(NULL, argv[1], "r");
	// Set full variable name
	if(arridx >= 0) {
		sprintf(fullname1, "%s[%d]", varname, arridx);
		sprintf(fullname2, "%s[*]", varname);
	} else {
		sprintf(fullname1, "%s", varname);
		sprintf(fullname2, "%s", varname);
	}
	if(in) {	 // Read data (TH: This check must be redundant, since mkFile checks for NULL)
		while (fgets(line, LEN, in)) {
			if(sscanf(line, "%s %s %s", rvarname, dummy, rval) == 3) {
				if(0 == strcasecmp(rvarname, fullname1) ||
						0 == strcasecmp(rvarname, fullname2)) {
					contain = 1;
					break;
				}
			}
		}
	}
	for(int i = 1; i < argc - 1; i++) {
		if(0 == strcasecmp(argv[i], fullname1) ||
				0 == strcasecmp(argv[i], fullname2)) {
			sprintf(rval, "%s", argv[i + 1]);
			contain = 1;
			break;
		}
	}
	if(in) fclose(in);
	// Check for missing variables
	if(!contain) {
		if(strlen(defvalue) > 0) {
			sprintf(rval, "%s", defvalue);
		} else {
			sprintf(msg, "Missing variable %s!\n", fullname1);
			ERRMSG(msg);
		}
	}
	// Write info
	if(arridx < 0) printf("%s = %s\n", fullname1, rval);
	// Return values
	if(value) sprintf(value, "%s", rval);
	return atof(rval);
}

//***************************************************************************
void time2jsec(int year, int mon, int day, int hour, int min, int sec, double remain, double *jsec) {
	struct tm t0 = { .tm_year=100,				 .tm_mon=0,				.tm_mday=1,		.tm_hour=0,		 .tm_min=0,		.tm_sec=0		};
	struct tm t1 = { .tm_year=year - 1900, .tm_mon=mon - 1, .tm_mday=day, .tm_hour=hour, .tm_min=min, .tm_sec=sec };
	*jsec = (double) timegm(&t1) - (double) timegm(&t0) + remain;
}

void jsec2time(double jsec, int *year, int *mon, int *day, int *hour, int *min, int *sec, double *remain) {
	struct tm t0 = {.tm_year = 100, .tm_mon = 0, .tm_mday = 1, .tm_hour = 0, .tm_min = 0, .tm_sec = 0};
	time_t jsec0 = (time_t) jsec + timegm(&t0);
	struct tm *t1 = gmtime(&jsec0);
	*year		= t1->tm_year + 1900;
	*mon		= t1->tm_mon + 1;
	*day		= t1->tm_mday;
	*hour		= t1->tm_hour;
	*min		= t1->tm_min;
	*sec		= t1->tm_sec;
	*remain	= jsec - floor(jsec);
}

//***************************************************************************
double timer(char const *name, char const *file, char const *func, int line, int mode) {
    #define MaxNumTimers 10
	static int    nt = 0; // number of timers active
	static double w0[MaxNumTimers]; // start time in seconds
	static int    l0[MaxNumTimers]; // source line where the timer was started
	struct timeval tim;
    double dt_w = 0;
	if(mode == 1) {										// Start new timer
		gettimeofday(&tim, NULL);
		w0[nt] = (double) tim.tv_sec + (double) tim.tv_usec*1e-6;
		l0[nt] = line;
		if((++nt) >= MaxNumTimers) ERRMSG("Too many timers! max. is " xstr(MaxNumTimers));
	} else {
		if (nt - 1 < 0) ERRMSG("Coding error!");							// Check timer index
		gettimeofday(&tim, NULL);									// Get time diff
		dt_w = (double) tim.tv_sec + (double) tim.tv_usec*1e-6 - w0[nt - 1];
        if (mode != -3) // -3:silent stop
		printf("Timer '%s' (%s, %s, l%d-%d): %.3f sec\n", name, file, func, l0[nt - 1], line, dt_w); // Write elapsed time
	}
	if(abs(mode) == 3) nt--;										// Stop timer
#undef MaxNumTimers
    return dt_w;
}

//***************************************************************************
void write_atm(char const *dirname, char const *filename, ctl_t *ctl, atm_t *atm) {
    if (ctl->checkmode) {
        printf("# skip writing target file name for atmospheric data: %s/%s\n", dirname, filename);
        return; // return before creating the file
    } // checkmode
    printf("Write atmospheric data: %s/%s", dirname, filename);
	FILE* out = mkFile(dirname, filename, "w");
	// Write header
	fprintf(out,
			"# $1 = time (seconds since 2000-01-01T00:00Z)\n"
			"# $2 = altitude [km]\n"
			"# $3 = longitude [deg]\n"
			"# $4 = latitude [deg]\n"
			"# $5 = pressure [hPa]\n"
            "# $6 = temperature [K]\n");
	int n = 6;
	for(int ig = 0; ig < ctl->ng; ig++) fprintf(out, "# $%d = %s volume mixing ratio\n", ++n, ctl->emitter[ig]);
	for(int iw = 0; iw < ctl->nw; iw++) fprintf(out, "# $%d = window %d: extinction [1/km]\n", ++n, iw);
	// Write data
	for(int ip = 0; ip < atm->np; ip++) {
		if(ip == 0 || atm->time[ip] != atm->time[ip - 1]) fprintf(out, "\n");
		fprintf(out, "%.2f %g %g %g %g %g", atm->time[ip], atm->z[ip],
				atm->lon[ip], atm->lat[ip], atm->p[ip], atm->t[ip]);
		for(int ig = 0; ig < ctl->ng; ig++) fprintf(out, " %g", atm->q[ig][ip]);
		for(int iw = 0; iw < ctl->nw; iw++) fprintf(out, " %g", atm->k[iw][ip]);
		fprintf(out, "\n");
	}
	fclose(out);
}

//***************************************************************************
void write_atm_rfm(char const *filename, ctl_t const *const ctl, atm_t const *const atm) {
	printf("Write RFM data: %s\n", filename);
	FILE *out = mkFile(NULL, filename, "w");
	// Write data
	fprintf(out, "%d\n", atm->np);
	fprintf(out, "*HGT [km]\n");
	for(int ip = 0; ip < atm->np; ip++) fprintf(out, "%g\n", atm->z[ip]);
	fprintf(out, "*PRE [mb]\n");
	for(int ip = 0; ip < atm->np; ip++) fprintf(out, "%g\n", atm->p[ip]);
	fprintf(out, "*TEM [K]\n");
	for(int ip = 0; ip < atm->np; ip++) fprintf(out, "%g\n", atm->t[ip]);
	for(int ig = 0; ig < ctl->ng; ig++) {
		fprintf(out, "*%s [ppmv]\n", ctl->emitter[ig]);
		for(int ip = 0; ip < atm->np; ip++) fprintf(out, "%g\n", atm->q[ig][ip] *1e6);
	}
	fprintf(out, "*END\n");
	fclose(out);
}

//***************************************************************************
void idx2name(ctl_t *ctl, int idx, char *quantity) {
	if(idx == IDXP) sprintf(quantity, "PRESSURE");
	if(idx == IDXT) sprintf(quantity, "TEMPERATURE");
	for(int ig = 0; ig < ctl->ng; ig++)
		if(idx == IDXQ(ig)) sprintf(quantity, "%s", ctl->emitter[ig]);
	for(int iw = 0; iw < ctl->nw; iw++)
		if(idx == IDXK(iw)) sprintf(quantity, "EXTINCT_WINDOW%d", iw);
}

void write_matrix(char const *dirname, char const *filename,
		ctl_t *ctl, gsl_matrix *matrix,
		atm_t *atm, obs_t *obs,
		char const *rowspace, char const *colspace, char const *sort) {
	char quantity[LEN];
	int *cida, *ciqa, *cipa, *cira, *rida, *riqa, *ripa, *rira;
	size_t i, j, nc, nr;
	// Check output flag
	if(!ctl->write_matrix) return;
	ALLOC(cida, int, M);
	ALLOC(ciqa, int, N);
	ALLOC(cipa, int, N);
	ALLOC(cira, int, M);
	ALLOC(rida, int, M);
	ALLOC(riqa, int, N);
	ALLOC(ripa, int, N);
	ALLOC(rira, int, M);
	// Open output file
	printf("Write matrix: %s/%s", dirname, filename);
	FILE *out = mkFile(dirname, filename, "w");
	// Write header (row space)
	if(rowspace[0] == 'y') {
		fprintf(out,
				"#	$1 = Row: index (measurement space)\n"
				"#	$2 = Row: channel wavenumber [cm^-1]\n"
				"#	$3 = Row: time (seconds since 2000-01-01T00:00Z)\n"
				"#	$4 = Row: view point altitude [km]\n"
				"#	$5 = Row: view point longitude [deg]\n"
				"#	$6 = Row: view point latitude [deg]\n");
		nr = obs2y(ctl, obs, NULL, rida, rira);			// Get number of rows
	} else {
		fprintf(out,
				"#	$1 = Row: index (state space)\n"
				"#	$2 = Row: name of quantity\n"
				"#	$3 = Row: time (seconds since 2000-01-01T00:00Z)\n"
				"#	$4 = Row: altitude [km]\n"
				"#	$5 = Row: longitude [deg]\n" "# $6 = Row: latitude [deg]\n");
		nr = atm2x(ctl, atm, NULL, riqa, ripa);			// Get number of rows
	}
	// Write header (column space)
	if(colspace[0] == 'y') {
		fprintf(out,
				"#	$7 = Col: index (measurement space)\n"
				"#	$8 = Col: channel wavenumber [cm^-1]\n"
				"#	$9 = Col: time (seconds since 2000-01-01T00:00Z)\n"
				"# $10 = Col: view point altitude [km]\n"
				"# $11 = Col: view point longitude [deg]\n"
				"# $12 = Col: view point latitude [deg]\n");
		nc = obs2y(ctl, obs, NULL, cida, cira);			// Get number of columns
	} else {
		fprintf(out,
				"#	$7 = Col: index (state space)\n"
				"#	$8 = Col: name of quantity\n"
				"#	$9 = Col: time (seconds since 2000-01-01T00:00Z)\n"
				"# $10 = Col: altitude [km]\n"
				"# $11 = Col: longitude [deg]\n" "# $12 = Col: latitude [deg]\n");
		nc = atm2x(ctl, atm, NULL, ciqa, cipa);			// Get number of columns
	}
	fprintf(out, "# $13 = Matrix element\n\n");		// Write header entry
	// Write matrix data
	i = j = 0;
	while (i < nr && j < nc) {
		if(gsl_matrix_get(matrix, i, j) != 0) {			// Check matrix value
			// Write info about the row
			if(rowspace[0] == 'y') {
				fprintf(out, "%d %g %.2f %g %g %g",
						(int) i, ctl->nu[rida[i]],
						obs->time[rira[i]], obs->vpz[rira[i]],
						obs->vplon[rira[i]], obs->vplat[rira[i]]);
			} else {
				idx2name(ctl, riqa[i], quantity);
				fprintf(out, "%d %s %.2f %g %g %g", (int) i, quantity,
						atm->time[ripa[i]], atm->z[ripa[i]],
						atm->lon[ripa[i]], atm->lat[ripa[i]]);
			}
			// Write info about the column
			if(colspace[0] == 'y') {
				fprintf(out, " %d %g %.2f %g %g %g",
						(int) j, ctl->nu[cida[j]],
						obs->time[cira[j]], obs->vpz[cira[j]],
						obs->vplon[cira[j]], obs->vplat[cira[j]]);
			} else {
				idx2name(ctl, ciqa[j], quantity);
				fprintf(out, " %d %s %.2f %g %g %g", (int) j, quantity,
						atm->time[cipa[j]], atm->z[cipa[j]],
						atm->lon[cipa[j]], atm->lat[cipa[j]]);
			}
			fprintf(out, " %g\n", gsl_matrix_get(matrix, i, j));			 // Write matrix entry
		}
		if(sort[0] == 'r') {		 // Set matrix indices
			j++;
			if(j >= nc) {
				j = 0;
				i++;
				fprintf(out, "\n");
			}
		} else {
			i++;
			if(i >= nr) {
				i = 0;
				j++;
				fprintf(out, "\n");
			}
		}
	}
	fclose(out);
	free(cida);
	free(ciqa);
	free(cipa);
	free(cira);
	free(rida);
	free(riqa);
	free(ripa);
	free(rira);
}

//***************************************************************************
void write_obs(char const *dirname, char const *filename, ctl_t *ctl, obs_t *obs) {
    if (ctl->checkmode) {
        printf("# skip writing target file name for observation data: %s/%s\n", dirname, filename);
        return; // return before creating the file
    } // checkmode
	printf("Write observation data: %s/%s\n", dirname, filename);
	FILE* out = mkFile(dirname, filename, "w");
	// Write header
	fprintf(out,
			"# $1 = time (seconds since 2000-01-01T00:00Z)\n"
			"# $2 = observer altitude [km]\n"
			"# $3 = observer longitude [deg]\n"
			"# $4 = observer latitude [deg]\n"
			"# $5 = view point altitude [km]\n"
			"# $6 = view point longitude [deg]\n"
			"# $7 = view point latitude [deg]\n"
			"# $8 = tangent point altitude [km]\n"
			"# $9 = tangent point longitude [deg]\n"
			"# $10 = tangent point latitude [deg]\n");
	int n = 10;
    char const *rad_or_bt = (ctl->write_bbt)? "brightness temperature [K]" : "radiance [W/(m^2 sr cm^-1)]";
    for(int id = 0; id < ctl->nd; id++) {
        fprintf(out, "# $%d = channel %g: %s\n", ++n, ctl->nu[id], rad_or_bt);
    } // id
	for(int id = 0; id < ctl->nd; id++) {
        ++n;
        if ((ctl->nd < 65) || (id < 1) || (id > ctl->nd - 2)) {
            fprintf(out, "# $%d = channel %g: transmittance\n", n, ctl->nu[id]);
        } else if (1 == id) {
            fprintf(out, "# $%d through $%d transmittance\n", n, n + ctl->nd - 3);
        }
    } // id
	for(int ir = 0; ir < obs->nr; ir++) { // Write data
		if(ir == 0 || (NR > 1 && obs->time[ir] != obs->time[ir - 1])) fprintf(out, "\n");
		fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g",
				obs->time[ir],
				obs->obsz[ir], obs->obslon[ir], obs->obslat[ir],
				obs->vpz[ir],  obs->vplon[ir],	obs->vplat[ir],
				obs->tpz[ir],  obs->tplon[ir],	obs->tplat[ir]);
		for(int id = 0; id < ctl->nd; id++) fprintf(out, " %g", obs->rad[ir][id]);
		for(int id = 0; id < ctl->nd; id++) fprintf(out, " %g", obs->tau[ir][id]);
        fprintf(out, "\n");
	}
	fclose(out);
}

//***************************************************************************
void x2atm(ctl_t const *const ctl, gsl_vector *x, atm_t *atm) {
	size_t n = 0;
	x2atm_help(atm, ctl->retp_zmin, ctl->retp_zmax, atm->p, x, &n);							// Set pressure
	x2atm_help(atm, ctl->rett_zmin, ctl->rett_zmax, atm->t, x, &n);							// Set temperature
	for(int ig = 0; ig < ctl->ng; ig++) x2atm_help(atm, ctl->retq_zmin[ig], ctl->retq_zmax[ig], atm->q[ig], x, &n);	// Set volume mixing ratio
	for(int iw = 0; iw < ctl->nw; iw++) x2atm_help(atm, ctl->retk_zmin[iw], ctl->retk_zmax[iw], atm->k[iw], x, &n);	// Set extinction
}

void x2atm_help(atm_t *atm, double zmin, double zmax, double *value, gsl_vector *x, size_t *n) {
	for(int ip = 0; ip < atm->np; ip++) {		// Extract state vector elements
		if((atm->z[ip] >= zmin) && (atm->z[ip] <= zmax)) {
			value[ip] = gsl_vector_get(x, *n);
			(*n)++;
		}
	}
}

//***************************************************************************
size_t atm2x(ctl_t const *const ctl, atm_t const *const atm, gsl_vector *x, int *iqa, int *ipa) {
	size_t n = 0;
	atm2x_help(atm, ctl->retp_zmin, ctl->retp_zmax, atm->p, IDXP, x, iqa, ipa, &n);		// Add pressure
	atm2x_help(atm, ctl->rett_zmin, ctl->rett_zmax, atm->t, IDXT, x, iqa, ipa, &n);		// Add temperature
	for(int ig = 0; ig < ctl->ng; ig++) {																							// Add volume mixing ratios
		atm2x_help(atm, ctl->retq_zmin[ig], ctl->retq_zmax[ig], atm->q[ig], IDXQ(ig), x, iqa, ipa, &n);
	}
	for(int iw = 0; iw < ctl->nw; iw++) {																						 // Add extinction
		atm2x_help(atm, ctl->retk_zmin[iw], ctl->retk_zmax[iw], atm->k[iw], IDXK(iw), x, iqa, ipa, &n);
	}
	return n;
}

void atm2x_help(atm_t const *const atm, double zmin, double zmax, double const *const value, int val_iqa, gsl_vector *x, int *iqa, int *ipa, size_t *n) {
	for(int ip = 0; ip < atm->np; ip++) {																						 // Add elements to state vector
		if(atm->z[ip] >= zmin && atm->z[ip] <= zmax) {
			if(x)		gsl_vector_set(x, *n, value[ip]);
			if(iqa) iqa[*n] = val_iqa;
			if(ipa) ipa[*n] = ip;
			(*n)++;
		}
	}
}

//***************************************************************************
void y2obs(ctl_t *ctl, gsl_vector *y, obs_t *obs) {
	size_t m = 0;
	for(int ir = 0; ir < obs->nr; ir++) { // Decompose measurement vector
		for(int id = 0; id < ctl->nd; id++) {
			if(gsl_finite(obs->rad[ir][id])) {
				obs->rad[ir][id] = gsl_vector_get(y, m);
				m++;
			}
		}
	}
}

size_t obs2y(ctl_t const *const ctl, obs_t const *const obs, gsl_vector *y, int *ida, int *ira) {
	size_t m = 0;
	for(int ir = 0; ir < obs->nr; ir++) { // Determine measurement vector
		for(int id = 0; id < ctl->nd; id++) {
			if(gsl_finite(obs->rad[ir][id])) {
				if(y)		gsl_vector_set(y, m, obs->rad[ir][id]);
				if(ida) ida[m] = id;
				if(ira) ira[m] = ir;
				++m;
			}
		}
	}
	return m;
}
