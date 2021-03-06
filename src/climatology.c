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
  Prepare atmospheric data file from climatological data.
*/

#include "jurassic.h"
#include "gsl/gsl_rng.h"

int main(
  int argc,
  char *argv[]) {

  static atm_t atm;
  static ctl_t ctl;

  double dt, dz, t, t0, t1, z, z0, z1, dpress, dtemp;

  int ip, rand;
  
  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <atm>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  t0 = scan_ctl(argc, argv, "T0", -1, "0", NULL);
  t1 = scan_ctl(argc, argv, "T1", -1, "0", NULL);
  dt = scan_ctl(argc, argv, "DT", -1, "1", NULL);
  z0 = scan_ctl(argc, argv, "Z0", -1, "0", NULL);
  z1 = scan_ctl(argc, argv, "Z1", -1, "90", NULL);
  dz = scan_ctl(argc, argv, "DZ", -1, "1", NULL);
  rand = (int)scan_ctl(argc, argv, "RAND", -1, "0", NULL);
  
  /* Set atmospheric grid... */
  for (t = t0; t <= t1; t += dt)
    for (z = z0; z <= z1; z += dz) {
      atm.time[atm.np] = t;
      atm.z[atm.np] = z;
      if ((++atm.np) >= NP)
	ERRMSG("Too many atmospheric grid points!");
    }

  /* Interpolate climatological data... */
  climatology(&ctl, &atm);

  /* Random modifications of atmospheric profiles... */
  if(rand) {
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc (gsl_rng_default);
    for(ip=0; ip<atm.np; ip++) {
      if(ip==0 || atm.time[ip-1]!=atm.time[ip]) {
	dpress=0.05-0.1*gsl_rng_uniform_pos(r);
	dtemp=30.-60.*gsl_rng_uniform_pos(r);
      }
      atm.p[ip]*=(1.0+dpress);
      atm.t[ip]+=dtemp;
    }
    gsl_rng_free(r);
  }
  
  /* Write data to disk... */
  write_atm(".", argv[2], &ctl, &atm);

  return EXIT_SUCCESS;
}
