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
	JURASSIC library declarations.

	\mainpage

	The JUelich RApid Spectral SImulation Code (JURASSIC) is a fast radiative
	transfer model for the mid-infrared spectral region.

	This reference manual provides information on the algorithms
	and data structures used in the code. Further information can be found at:
http://www.fz-juelich.de/ias/jsc/jurassic
*/
#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>

/* ------------------------------------------------------------
	 Macros...
	 ------------------------------------------------------------ */
/*! Allocate memory. */
#define ALLOC(ptr, type, n)				\
	if((ptr=(type*)malloc((size_t)(n)*sizeof(type)))==NULL)	\
ERRMSG("Out of memory!");

/*! Compute Cartesian distance between two vectors. */
#define DIST(a, b) sqrt(DIST2(a, b))

/*! Compute squared distance between two vectors. */
#define DIST2(a, b)							\
	((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))

/*! Compute dot product of two vectors. */
#define DOTP(a, b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

/*! Print error message and quit program. */
#define ERRMSG(msg) {							\
	printf("\nError (%s, %s, l%d): %s\n\n",				\
			__FILE__, __func__, __LINE__, msg);				\
	exit(EXIT_FAILURE);							\
}

/*! Compute exponential interpolation. */
#define EXP(x0, y0, x1, y1, x)					\
	(((y0)>0 && (y1)>0)						\
	 ? ((y0)*exp(log((y1)/(y0))/((x1)-(x0))*((x)-(x0))))					\
	 : LIN(x0, y0, x1, y1, x))

/*! Compute linear interpolation. */
#define LIN(x0, y0, x1, y1, x) ((y0) + ((x) - (x0))*((y1) - (y0))/((x1) - (x0)))

/*! Compute norm of a vector. */
#define NORM(a) sqrt(DOTP(a, a))

/*! Print macro for debugging. */
#define PRINT(format, var)						\
	printf("Print (%s, %s, l%d): %s= " format "\n",				\
			__FILE__, __func__, __LINE__, #var, var);

/*! Start or stop a timer. */
#define TIMER(name, mode) timer(name, __FILE__, __func__, __LINE__, mode)

/*! Read string tokens. */
#define TOK(line, tok, format, var, saveptr) {			\
	if(((tok)=strtok_r((line), " \t",saveptr))) {			\
		if(sscanf(tok, format, &(var))!=1) continue;	\
	} else ERRMSG("Error while reading!");		\
}

#define __deprecated__ __attribute__((deprecated))

// double eip(const double, const double, const double, const double, const double) __deprecated__;
// double lip(const double, const double, const double, const double, const double) __deprecated__;

/* ------------------------------------------------------------
	 Constants...
	 ------------------------------------------------------------ */

/*! First spectroscopic constant (c_1 = 2 h c^2) [W/(m^2 sr cm^-4)]. */
#define C1 1.19104259e-8

/*! Second spectroscopic constant (c_2 = h c / k) [K/cm^-1]. */
#define C2 1.43877506

/*! Standard gravity [m/s^2]. */
#define G0 9.80665

/*! Standard pressure [hPa]. */
#define P0 1013.25

/*! Standard temperature [K]. */
#define T0 273.15

/*! Mean radius of Earth [km]. */
#define RE 6367.421

/*! Mass of Earth [kg]. */
#define ME 5.976e24



/* ------------------------------------------------------------
	 Dimensions...
	 ------------------------------------------------------------ */

/*! Maximum number of radiance channels. */
#ifndef ND
  #define ND 32 
#endif

/*! Maximum number of emitters. */
#ifndef NG
  #define NG 27
#endif  

/*! Maximum number of atmospheric data points. */
#define NP 9600

/*! Maximum number of ray paths. */
#define NR 1088

/*! Maximum number of spectral windows. */
#define NW 1

/*! Maximum length of ASCII data lines. */
#define LEN 5000

/*! Maximum size of measurement vector. */
#define M (NR*ND)

/*! Maximum size of state vector. */
#define N (NQ*NP)

/*! Maximum number of quantities. */
#define NQ (2+NG+NW)

/*! Maximum number of LOS points. */
#define NLOS 400

/*! Maximum number of shape function grid points. */
#define NSHAPE 2048

/*! Number of ray paths used for FOV calculations. */
#define NFOV 5

/*! Maximum number of pressure levels in emissivity tables. */
#define TBLNP 40

/*! Maximum number of temperatures in emissivity tables. */
#define TBLNT 30

/*! Maximum number of column densities in emissivity tables. */
#define TBLNU 304

/*! Maximum number of source function temperature levels. */
#define TBLNS 1201

/*! Maximum number of RFM spectral grid points. */
#define RFMNPTS 10000000

/*! Maximum length of RFM data lines. */
#define RFMLINE 100000

/* ------------------------------------------------------------
	 Quantity indices...
	 ------------------------------------------------------------ */

/*! Index for pressure. */
#define IDXP 0

/*! Index for temperature. */
#define IDXT 1

/*! Indices for volume mixing ratios. */
#define IDXQ(ig) (2+ig)

/*! Indices for extinction. */
#define IDXK(iw) (2+ctl->ng+iw)

/* ------------------------------------------------------------
	 Structs...
	 ------------------------------------------------------------ */

typedef struct { /// Atmospheric data. //////////////////////////////////////////
	double time[NP];			 /// Time (seconds since 2000-01-01T00:00Z).
	double z[NP];					 /// Altitude [km].
	double lon[NP];				 /// Longitude [deg].
	double lat[NP];				 /// Latitude [deg].
	double p[NP];					 /// Pressure [hPa].
	double t[NP];					 /// Temperature [K].
	double q[NG][NP];			 /// Volume mixing ratio.
	double k[NW][NP];			 /// Extinction [1/km].
	int np;                /// Number of data points.
	int init;							 /// Init flag for interpolation (0=no, 1=yes).
} atm_t; // ////////////////////////////////////////////////////////////////////

/*! Forward model control parameters. */
typedef struct {

	/*! Number of emitters. */
	int ng;

	/*! Name of each emitter. */
	char emitter[NG][LEN];

	/*! Number of radiance channels. */
	int nd;

	/*! Number of spectral windows. */
	int nw;

	/*! Centroid wavenumber of each channel [cm^-1]. */
	double nu[ND];

	/*! Window index of each channel. */
	int window[ND];

	/*! Basename for table files and filter function files. */
	char tblbase[LEN];

	/*! Reference height for hydrostatic pressure profile (-999 to skip) [km]. */
	double hydz;

	/*! Compute CO2 continuum (0=no, 1=yes). */
	int ctm_co2;

	/*! Compute H2O continuum (0=no, 1=yes). */
	int ctm_h2o;

	/*! Compute N2 continuum (0=no, 1=yes). */
	int ctm_n2;

	/*! Compute O2 continuum (0=no, 1=yes). */
	int ctm_o2;

	/*! Interpolation method (1=profile, 2=satellite track, 3=Lagrangian grid). */
	int ip;

	/*! Influence length for vertical interpolation [km]. */
	double cz;

	/*! Influence length for horizontal interpolation [km]. */
	double cx;

	/*! Take into account refractivity (0=no, 1=yes). */
	int refrac;

	/*! Maximum step length for raytracing [km]. */
	double rayds;

	/*! Vertical step length for raytracing [km]. */
	double raydz;

	/*! Field-of-view data file. */
	char fov[LEN];

	/*! Minimum altitude for pressure retrieval [km]. */
	double retp_zmin;

	/*! Maximum altitude for pressure retrieval [km]. */
	double retp_zmax;

	/*! Minimum altitude for temperature retrieval [km]. */
	double rett_zmin;

	/*! Maximum altitude for temperature retrieval [km]. */
	double rett_zmax;

	/*! Minimum altitude for volume mixing ratio retrieval [km]. */
	double retq_zmin[NG];

	/*! Maximum altitude for volume mixing ratio retrieval [km]. */
	double retq_zmax[NG];

	/*! Minimum altitude for extinction retrieval [km]. */
	double retk_zmin[NW];

	/*! Maximum altitude for extinction retrieval [km]. */
	double retk_zmax[NW];

	/*! Use brightness temperature instead of radiance (0=no, 1=yes). */
	int write_bbt;

	/*! Write matrix file (0=no, 1=yes). */
	int write_matrix;

	/*! Forward model (1=CGA, 2=EGA, 3=RFM). */
	int formod;

	/*! Path to RFM binary. */
	char rfmbin[LEN];

	/*! HITRAN file for RFM. */
	char rfmhit[LEN];

	/*! Emitter cross-section files for RFM. */
	char rfmxsc[NG][LEN];

	/*! Use GPU-accelerated formod implementation (0=no, 1=yes) */
	int useGPU;
    
    /*! do not perform input, computation, nor output, just make sure files are there */
	int checkmode;

	/*! MPI rank information */
	int MPIglobrank;  // global rank
	int MPIlocalrank; // node-local Rank

	/* binary IO */
    int read_binary;
    int write_binary;

    /*! Shared memory controler for GPU kernels */
    int gpu_nbytes_shared_memory;

} ctl_t;


/*! Point on the Line-of-sight data without storing */
typedef struct {
	double z;		/*! Altitude [km]. */
	double lon;	/*! Longitude [deg]. */
	double lat;	/*! Latitude [deg]. */
	double p;		/*! Pressure [hPa]. */
	double t;		/*! Temperature [K]. */
	double q[NG];	/*! Volume mixing ratio. */
	double k[NW];	/*! Extinction [1/km]. */
	double ds;	/*! Segment length [km]. */
	double u[NG];	/*! Column density [molecules/cm^2]. */
#ifdef CURTIS_GODSON
	double cgp[NG];	/*! Curtis-Godson pressure [hPa]. */
	double cgt[NG];	/*! Curtis-Godson temperature [K]. */
	double cgu[NG];	/*! Curtis-Godson column density [molecules/cm^2]. */
#endif
#ifdef GPUDEBUG
	int ip, ir;  // debug helpers
#endif
} pos_t;

typedef struct { /// Observation geometry and radiance data. ///////////////////
	double time[NR];		/// Time (seconds since 2000-01-01T00:00Z). 
	double obsz[NR];		/// Observer altitude [km]. 
	double obslon[NR];	/// Observer longitude [deg]. 
	double obslat[NR];	/// Observer latitude [deg]. 
	double vpz[NR];			/// View point altitude [km]. 
	double vplon[NR];		/// View point longitude [deg]. 
	double vplat[NR];		/// View point latitude [deg]. 
	double tpz[NR];			/// Tangent point altitude [km]. 
	double tplon[NR];		/// Tangent point longitude [deg]. 
	double tplat[NR];		/// Tangent point latitude [deg]. 
	double tau[NR][ND]; /// Transmittance of ray path.		// transposed
	double rad[NR][ND]; /// Radiance [W/(m^2 sr cm^-1)].	// transposed
	int nr;							/// Number of ray paths.
} obs_t; // ////////////////////////////////////////////////////////////////////

typedef float real_tblND_t;

/*! Emissivity look-up tables. */
typedef struct {

	/*! Number of pressure levels. */
	int32_t np[NG][ND];

	/*! Number of temperatures. */
	int32_t nt[NG][TBLNP][ND];

	/*! Number of column densities. */
	int32_t nu[NG][TBLNP][TBLNT][ND];

	/*! Pressure [hPa]. */
	double p[NG][TBLNP][ND];

	/*! Temperature [K]. */
	double t[NG][TBLNP][TBLNT][ND];

	/*! Column density [molecules/cm^2]. */
	real_tblND_t u[NG][TBLNP][TBLNT][TBLNU][ND];
    
	/*! Emissivity. */
	real_tblND_t eps[NG][TBLNP][TBLNT][TBLNU][ND];

	/*! Source function radiance [W/(m^2 sr cm^-1)]. */
	double sr[TBLNS][ND];

	/*! Source function temperature [K]. */
	double st[TBLNS];

#ifdef  FAST_INVERSE_OF_U
	/*! u0inv[g][p][t][d] * u[g][p][t][0][d] == 1 must hold! */ // FAST_INVERSE_OF_U
	double u0inv[NG][TBLNP][TBLNT][ND];                         // FAST_INVERSE_OF_U
    /*! We assume a logarithmic increment by 2^(1/6) */         // FAST_INVERSE_OF_U
#endif

} tbl_t;

/* ------------------------------------------------------------
	 Functions...
	 ------------------------------------------------------------ */

/*! Compose state vector or parameter vector. */
size_t atm2x(
		ctl_t const *ctl,
		atm_t const *atm,
		gsl_vector * x,
		int *iqa,
		int *ipa);

/*! Add elements to state vector. */
void atm2x_help(
		atm_t const *atm,
		double zmin,
		double zmax,
		double const *value,
		int val_iqa,
		gsl_vector * x,
		int *iqa,
		int *ipa,
		size_t * n);

/*! Compute brightness temperature. */
double brightness(
		double const rad,
		double const nu) __deprecated__;

/*! Convert Cartesian coordinates to geolocation. */
// void cart2geo(
// 		double const *x,
// 		double *z,
// 		double *lon,
// 		double *lat) __deprecated__;

/*! Interpolate climatological data. */
void climatology(
		ctl_t const * ctl,
		atm_t * atm_mean);

/*! Compute carbon dioxide continuum (optical depth). */
double ctmco2(
		double nu,
		double p,
		double t,
		double u) __deprecated__;

/*! Compute water vapor continuum (optical depth). */
double ctmh2o(
		double nu,
		double p,
		double t,
		double q,
		double u) __deprecated__;

/*! Compute nitrogen continuum (absorption coefficient). */
double ctmn2(
		double nu,
		double p,
		double t) __deprecated__;

/*! Compute oxygen continuum (absorption coefficient). */
double ctmo2(
		double nu,
		double p,
		double t) __deprecated__;

/*! Copy and initialize atmospheric data. */
void copy_atm(
		ctl_t const *ctl,
		atm_t * atm_dest,
		atm_t const *atm_src,
		int const init);

/*! Copy and initialize observation data. */
void copy_obs(
		ctl_t const *ctl,
		obs_t * obs_dest,
		obs_t const *obs_src,
		int const init);

/*! Find index of an emitter. */
int find_emitter(
		ctl_t const * const ctl,
		char const * const emitter);

/*! Determine ray paths and compute radiative transfer. */
void formod(
		ctl_t const *ctl,
		atm_t * atm,
		obs_t * obs);

/*! Apply field of view convolution. */
void formod_fov(
		ctl_t const *ctl,
		obs_t * obs);

/*! Compute radiative transfer for a pencil beam. */
void formod_pencil(
		ctl_t const *ctl,
		atm_t * atm,
		obs_t * obs,
		int const ir);

/*! Apply RFM for radiative transfer calculations. */
void formod_rfm(
		ctl_t const *ctl,
		atm_t * atm,
		obs_t * obs);

/*! Compute Planck source function. */
void formod_srcfunc(
		ctl_t const *ctl,
		tbl_t const *tbl,
		double const t,
		double *src);

/*! Convert geolocation to Cartesian coordinates. */
// void geo2cart(
// 		double z,
// 		double lon,
// 		double lat,
// 		double *x) __deprecated__;

// /*! Determine gravity of Earth. */
// double gravity(
// 		double z,
// 		double lat) __deprecated__;

/*! Set hydrostatic equilibrium. */
void hydrostatic(
		ctl_t const *ctl,
		atm_t * atm) __deprecated__;

/*! Set hydrostatic equilibrium for individual profile. */
void hydrostatic_1d(
		ctl_t const *ctl,
		atm_t * atm,
		int const ip0,
		int const ip1) __deprecated__;

/*! Determine name of state vector quantity for given index. */
void idx2name(
		ctl_t * ctl,
		int idx,
		char *quantity);

/*! Initialize look-up tables. */
void init_tbl(
		ctl_t const * ctl,
		tbl_t * tbl);

/*! Interpolate complete atmospheric data set. */
void intpol_atm(
		ctl_t * ctl,
		atm_t * atm_dest,
		atm_t * atm_src);

/*! Interpolate atmospheric data for given geolocation. */
void intpol_atm_geo(
		ctl_t const *ctl,
		atm_t *atm,
		double const z0,
		double const lon0,
		double const lat0,
		double *p,
		double *t,
		double *q,
		double *k);

/*! Interpolate 1D atmospheric data (vertical profile). */
void intpol_atm_1d(
		ctl_t const *ctl,
		atm_t const *atm,
		int const idx0,
		int const n,
		double const z0,
		double *p,
		double *t,
		double *q,
		double *k);

/*! Interpolate 2D atmospheric data (satellite track). */
void intpol_atm_2d(
		ctl_t const *ctl,
		atm_t *atm,
		double const z0,
		double const lon0,
		double const lat0,
		double *p,
		double *t,
		double *q,
		double *k);

/*! Interpolate 3D atmospheric data (Lagrangian grid). */
void intpol_atm_3d(
		ctl_t const *ctl,
		atm_t *atm,
		double const z0,
		double const lon0,
		double const lat0,
		double *p,
		double *t,
		double *q,
		double *k);

/*! Interpolate emissivity from look-up tables. */
double intpol_tbl_eps(
		tbl_t const *tbl,
		int const ig,
		int const id,
		int const ip,
		int const it,
		double const u);

/*! Interpolate column density from look-up tables. */
double intpol_tbl_u(
		tbl_t const *tbl,
		int const ig,
		int const id,
		int const ip,
		int const it,
		double const eps);

/*! Convert seconds to date. */
void jsec2time(
		double const jsec,
		int *year,
		int *mon,
		int *day,
		int *hour,
		int *min,
		int *sec,
		double *remain);

/*! Compute Jacobians. */
void kernel(
		ctl_t const *ctl,
		atm_t *atm,
		obs_t *obs,
		gsl_matrix * k) __deprecated__;

/*! Compose measurement vector. */
size_t obs2y(
		ctl_t const *ctl,
		obs_t const *obs,
		gsl_vector * y,
		int *ida,
		int *ira);

/*! Compute Planck function. */
double planck(
		double const t,
		double const nu);

void altitude_range(
		atm_t const *atm,
		double *zmin,
		double *zmax) __deprecated__;

/*! Change segment lengths according to trapezoid rule */
void trapezoid_rule(
		int const np,
		double ds[]) __deprecated__;

/*! Read atmospheric data. */
void read_atm(
		const char *dirname,
		const char *filename,
		ctl_t * ctl,
		atm_t * atm);

/*! Read forward model control parameters. */
void read_ctl(
		int argc,
		char *argv[],
		ctl_t * ctl);

/*! Read matrix. */
void read_matrix(
		const char *dirname,
		const char *filename,
		gsl_matrix * matrix);

/*! Read observation data. */
void read_obs(
		const char *dirname,
		const char *filename,
		ctl_t * ctl,
		obs_t * obs);

/*! Read observation data in RFM format. */
double read_obs_rfm(
		const char *basename,
		double z,
		double *nu,
		double *f,
		int n);

/*! Read RFM spectrum. */
void read_rfm_spec(
		const char *filename,
		double *nu,
		double *rad,
		int *npts);

/*! Read shape function. */
int read_shape(
		const char *filename,
		double *x,
		double *y,
        int const checkmode);

// /*! Compute refractivity (return value is n - 1). */
// double refractivity(
// 		double p,
// 		double t) __deprecated__;

/*! Search control parameter file for variable entry. */
double scan_ctl(
		int argc,
		char *argv[],
		const char *varname,
		int arridx,
		const char *defvalue,
		char *value);

/*! Convert date to seconds. */
void time2jsec(
		int year,
		int mon,
		int day,
		int hour,
		int min,
		int sec,
		double remain,
		double *jsec);

/*! Measure wall-clock time. */
double timer(
		const char *name,
		const char *file,
		const char *func,
		int line,
		int mode);

/*! Write atmospheric data. */
void write_atm(
		const char *dirname,
		const char *filename,
		ctl_t * ctl,
		atm_t * atm);

/*! Write atmospheric data in RFM format. */
void write_atm_rfm(
		const char *filename,
		ctl_t const *ctl,
		atm_t const *atm);

/*! Write matrix. */
void write_matrix(
		const char *dirname,
		const char *filename,
		ctl_t * ctl,
		gsl_matrix * matrix,
		atm_t * atm,
		obs_t * obs,
		const char *rowspace,
		const char *colspace,
		const char *sort);

/*! Write observation data. */
void write_obs(
		const char *dirname,
		const char *filename,
		ctl_t * ctl,
		obs_t * obs);

/*! Decompose parameter vector or state vector. */
void x2atm(
		ctl_t const *ctl,
		gsl_vector *x,
		atm_t * atm);

/*! Extract elements from state vector. */
void x2atm_help(
		atm_t * atm,
		double zmin,
		double zmax,
		double *value,
		gsl_vector * x,
		size_t * n);

/*! Decompose measurement vector. */
void y2obs(ctl_t * ctl, gsl_vector * y, obs_t * obs);

/* Helpers */
FILE* mkFile(const char*, const char*, const char*);
void shell(const char*, const char*);
// double c01(double);

/* stringify the value of a macro, two expansion levels needed */
#define xstr(a) str(a)
#define str(b) #b
