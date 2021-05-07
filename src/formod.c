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
  JURASSIC forward model.
*/

#include <string.h> // memcpy
#include <omp.h>
#include "jurassic.h" // ctl_t, obs_t, atm_t, read_*, write_*, formod_*PU

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
	 int argc,
	 char *argv[]) {
  
  static ctl_t ctl; // why static?

  atm_t *atm;

  obs_t *obs;
  
  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <obs> <atm> <rad>");

  /* Allocate... */
  ALLOC(atm, atm_t, 1);
  ALLOC(obs, obs_t, 1);
  
  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  
  /* Read observation geometry... */
  read_obs(".", argv[2], &ctl, obs);

  /* Read atmospheric data... */
  read_atm(".", argv[3], &ctl, atm);
  // to verify that we read in the atms correctly
  //   write_atm(".", "atm_ref.tab", &ctl, atm);

  /* Call forward model... */
  // warmp-up, some initialization will happen during the first call
  if (0 == ctl.checkmode) TIMER("warm-up", 1);
  formod(&ctl, atm, obs);
  if (0 == ctl.checkmode) TIMER("warm-up", 3);

  /* Write radiance data... for reference */
  write_obs(".", argv[4], &ctl, obs);

#ifdef BENCHMARK_FORMOD

  int const niterations = GSL_MAX(1, ctl.useGPU*ctl.useGPU);
  if (niterations > 1) printf("# always run %d iterations for benchmarking\n", niterations);
  
  obs_t *obs_bench;
  ALLOC(obs_bench, obs_t, 1);
  memcpy(obs_bench, obs, sizeof(obs_t));

  ctl_t *ctl_bench;
  ALLOC(ctl_bench, ctl_t, 1);
  memcpy(ctl_bench, &ctl, sizeof(ctl_t));

#ifdef BENCH_FORMOD_SCALING_TESTS  
  // scaling tests: modify the number of rays and channels
  for(int nd = 1; nd <= ctl.nd; nd *= 2) {  ctl_bench->nd = nd;
    printf("# with channels\n# with %d channels measure formod time\n", nd);
    for(int nr = 1; nr <= obs->nr; nr *= 2) {  obs_bench->nr = nr;
        printf("\n\n\nscaling test: runs with %d rays and %d channels\n\n", nr, nd);
#else
    {{
#endif        

        double time_stats[3] = {0, 0, 0};
        int deviations = 0;

        // run some iterations to improve on the stats
        for(int it = 0; it < niterations; ++it) {
            if (0 == ctl.checkmode) TIMER("formod", 1);
            formod(ctl_bench, atm, obs_bench); // benchmarking this routine
            double const dt = (0 == ctl.checkmode) ? TIMER("formod", -3) : 0; // -3:silent timer
            time_stats[0] += 1;
            time_stats[1] += dt;
            time_stats[2] += dt*dt;
            
          if (0 == ctl.checkmode && 0 == it) { // only check after the 1st iteration
              { /* Compare obs_bench using obs as reference */
                  char const *rad_or_bt = (ctl.write_bbt)? "brightness temperature" : "radiance";
                  {
                      double mdev_tau_a = 0, mdev_rad_a = 0; int ndev_tau_a = 0, ndev_rad_a = 0;
                      for(int ir = 0; ir < obs_bench->nr; ++ir) {
                          double mdev_tau = 0, mdev_rad = 0; int ndev_tau = 0, ndev_rad = 0;
                          for(int id = 0; id < ctl_bench->nd; ++id) {
                              double const dev_tau = obs_bench->tau[ir][id] - obs->tau[ir][id];
                              mdev_tau = GSL_MAX_DBL(mdev_tau, fabs(dev_tau));
                              ndev_tau += (dev_tau != 0);
                              double const dev_rad = obs_bench->rad[ir][id] - obs->rad[ir][id];
                              mdev_rad = GSL_MAX_DBL(mdev_rad, fabs(dev_rad));
                              ndev_rad += (dev_rad != 0);
                          } // id
                          if (ndev_tau > 0) printf("# %d deviations in transmittance in ray #%i, largest %.1e\n", ndev_tau, ir, mdev_tau);
                          if (ndev_rad > 0) printf("# %d deviations in %s in ray #%i, largest %.1e\n", ndev_rad, rad_or_bt, ir, mdev_rad);
                          ndev_tau_a += (ndev_tau > 0);
                          ndev_rad_a += (ndev_rad > 0);
                          mdev_tau_a = GSL_MAX_DBL(mdev_tau_a, mdev_tau);
                          mdev_rad_a = GSL_MAX_DBL(mdev_rad_a, mdev_rad);
                      } // ir
                      if (ndev_tau_a > 0) printf("# %d deviations in transmittance, largest %.1e\n", ndev_tau_a, mdev_tau_a);
                      if (ndev_rad_a > 0) printf("# %d deviations in %s largest %.1e\n", ndev_rad_a, rad_or_bt, mdev_rad_a);
                      if (ndev_tau_a > 0 || ndev_rad_a > 0) ++deviations;
                  }
                  if (deviations > 0) { // now check also the transposed to see better where the errors are
                      double mdev_tau_a = 0, mdev_rad_a = 0; int ndev_tau_a = 0, ndev_rad_a = 0;
                      for(int id = 0; id < ctl_bench->nd; ++id) {
                          double mdev_tau = 0, mdev_rad = 0; int ndev_tau = 0, ndev_rad = 0;
                          for(int ir = 0; ir < obs_bench->nr; ++ir) {
                              double const dev_tau = obs_bench->tau[ir][id] - obs->tau[ir][id];
                              mdev_tau = GSL_MAX_DBL(mdev_tau, fabs(dev_tau));
                              ndev_tau += (dev_tau != 0);
                              double const dev_rad = obs_bench->rad[ir][id] - obs->rad[ir][id];
                              mdev_rad = GSL_MAX_DBL(mdev_rad, fabs(dev_rad));
                              ndev_rad += (dev_rad != 0);
                          } // ir
                          if (ndev_tau > 0) printf("# %d deviations in transmittance in channel #%i, largest %.1e\n", ndev_tau, id, mdev_tau);
                          if (ndev_rad > 0) printf("# %d deviations in %s in channel #%i, largest %.1e\n", ndev_rad, rad_or_bt, id, mdev_rad);
                          ndev_tau_a += (ndev_tau > 0);
                          ndev_rad_a += (ndev_rad > 0);
                          mdev_tau_a = GSL_MAX_DBL(mdev_tau_a, mdev_tau);
                          mdev_rad_a = GSL_MAX_DBL(mdev_rad_a, mdev_rad);
                      } // id
                      if (ndev_tau_a > 0) printf("# %d deviations in transmittance, largest %.1e\n", ndev_tau_a, mdev_tau_a);
                      if (ndev_rad_a > 0) printf("# %d deviations in %s largest %.1e\n", ndev_rad_a, rad_or_bt, mdev_rad_a);
                      if (ndev_tau_a > 0 || ndev_rad_a > 0) ++deviations;
                  }
                  printf("# Compare obs-results: %s and transmittance for %d rays times %d channels shows%s deviations\n", 
                                                    rad_or_bt, obs_bench->nr, ctl_bench->nd, deviations?"":" no");
              } // compare

          } // checkmode
          
          if (deviations) it = niterations + 9; // stop the loop early, we do not need to benchmark wrong results.
        } // it
        
        if (deviations) {
          printf("# timing results are not shown due to deviations (%d) in obs-results!", deviations);
        } else {
          if (niterations > 1) printf("# ran %d iterations for benchmark\n", (int)time_stats[0]);
          double const weight = (time_stats[0] > 0)? 1./time_stats[0] : 0;
          double const meantime = time_stats[1]*weight;
          double const variance = time_stats[2]*weight - meantime*meantime;
          double const sigma = (variance > 0)? sqrt(variance) : 0;
          fprintf(stderr, "# with %d rays and %d channels formod took %g +/- %g seconds on the %cPU\n",
                                  obs_bench->nr, ctl_bench->nd, meantime, sigma, ctl.useGPU?'G':'C');
                   printf("# with %d rays and %d channels formod took %g +/- %g seconds on the %cPU\n",
                                  obs_bench->nr, ctl_bench->nd, meantime, sigma, ctl.useGPU?'G':'C');
        } // deviations
  
    } // scaling tests: modify the number of rays
  } // scaling tests: modify the number of channels

  /* Free... */
  free(obs_bench);
  free(ctl_bench);

#endif // BENCHMARK_FORMOD

  free(atm);
  free(obs);
  
  return EXIT_SUCCESS;
}
