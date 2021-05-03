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
  use the internal hash function onto string provided
*/

#include <stdio.h> // printf
#include <stdint.h> // uint64_t

uint64_t jr_simple_string_hash(char const *string); // declaration ...
// ... will be linked from jurassic.o
  
int main(int argc, char *argv[]) {
  if (argc < 2) { printf("usage: hash <string>"); return 1; }
  printf("0x%lx\n", jr_simple_string_hash(argv[1]));
  return 0;
}
