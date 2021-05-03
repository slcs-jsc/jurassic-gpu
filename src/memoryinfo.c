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
  Show detailed memory consumption of structs with the current jurassic.h
*/

#include "jurassic.h"

int main(int argc, char *argv[]) {
  if (argc > 1) printf("# %s: command line arguments %s ... have no effect", __FILE__, argv[0]);
  
  printf("\njurassic.h is configured as  ND=%d  NG=%d  NP=%d  NR=%d  NW=%d\n", ND, NG, NP, NR, NW);
  printf("   tables are configured as  TBLNP=%d  TBLNT=%d  TBLNU=%d\n", TBLNP, TBLNT, TBLNU);
  printf("   tables are in FP%d (%s)\n", (int)(8*sizeof(real_tblND_t)), (8 == sizeof(real_tblND_t)) ? "double" : "float");
  printf("   NLOS=%d\n", NLOS);

  double const kByte = 1e-3, MByte = 1e-6, GByte = 1e-9;
#define showMemory(TYPE, UNIT) printf("%s \t takes %12.3f %s\n", #TYPE, UNIT*sizeof(TYPE), #UNIT);
  showMemory(ctl_t, MByte);
  showMemory(atm_t, kByte);
  showMemory(pos_t, kByte);
  showMemory(obs_t, MByte);
  showMemory(tbl_t, GByte);
#undef  showMemory
  
  return EXIT_SUCCESS;
}
