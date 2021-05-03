#pragma once

#include <stdio.h> // printf, fprintf, stderr, sprintf
#include <stdint.h> // uint64_t

#include "jr_simple_string_hash.h" // jr_simple_string_hash

#define BINARY_TABLES_HEADER_LEN 16384
#define BINARY_TABLES_FILENAME_LEN 256
#define BINARY_TABLES_VERSION 20200211

void jr_binary_tables_filename(char *filename) {
    // convention for file naming binary table files
    sprintf(filename, "bin.jurassic-fp%d-tables-g%d-p%d-T%d-u%d-d%d",
        (int)(sizeof(real_tblND_t)*8), NG, TBLNP, TBLNT, TBLNU, ND);
} // jr_binary_tables_filename

int jr_binary_tables_header(char header[], size_t const header_length, ctl_t const *ctl) {

  for(int c = 0; c < header_length; ++c) header[c] = 0; // clear

  char *h = header;
  h += sprintf(h, "JURASSIC 0  Binary Emissivity Tables\n"
       "version %d\n"
#ifdef SHOW_GIT_KEY
       "git_key 0 " xstr(SHOW_GIT_KEY) "\n"
#endif
       "%s   %d\n"
       "NG     %d gases\n"
       "TBLNP  %d pressures\n"
       "TBLNT  %d temperatures\n"
       "TBLNU  %d column_densities\n"
       "ND     %d detector_channels\n"
       "TBLNS  %d radiances",
       BINARY_TABLES_VERSION,
       (sizeof(real_tblND_t) > 4)?"double":"float", (int)sizeof(real_tblND_t),
       NG, TBLNP, TBLNT, TBLNU, ND, TBLNS);

  h += sprintf(h, "\n\nfile_size   %lld\nheader_size %lld\ntable_size  %lld",
          header_length + sizeof(tbl_t), header_length, sizeof(tbl_t));

  char const *const header_end = header + header_length - 16;
  h += sprintf(h, "\n\n\nng %d emitter gases:\n", ctl->ng);
  for(int ig = 0; ig < ctl->ng; ++ig) {
      if (h > header_end) return __LINE__;
      h += sprintf(h, "%s %i\n", ctl->emitter[ig], ig);
  } // ig
  h += sprintf(h, "\n\n\nnd %d channels [cm^-1]:\n", ctl->nd);
  for(int id = 0; id < ctl->nd; ++id) {
      if (h > header_end) return __LINE__;
      h += sprintf(h, "%.4f %i\n", ctl->nu[id], id);
  } // id
#ifdef  FAST_INVERSE_OF_U
  h += sprintf(h, "\n\nFAST_INVERSE_OF_U 1\n");
#else
  h += sprintf(h, "\n\nFAST_INVERSE_OF_U 0\n");
#endif
  h += sprintf(h, "\n\nheader_end 0\n"); // amke sure to finish with \n

  int const nbytes_written = h - header;
  printf("# binary header uses %.3f of %.3f kByte\n", 1e-3*nbytes_written, 1e-3*header_length);
  return 0;
} // jr_binary_tables_header

int jr_binary_tables_check_header(ctl_t const *ctl, char const header[], size_t *header_length) {
    int const echo = 1; // verbose output besides warnings and errors
    int nerrors = 0;
  
    int has_fast_inverse = 0;
    size_t header_size = 0, table_size = 0;
    if(header_length) *header_length = 0;
    int linenumber = 1;
    char gas_found[NG]; for(int ig = 0; ig < NG; ++ig) gas_found[ig] = 0;
    char nu_found[ND];  for(int id = 0; id < ND; ++id) nu_found[id] = 0;
    char *h = header; // set to the beginning of the header
    char varname[32];
    while (*h) { // Read from the buffer up to a null-char
      
        int64_t err = 0;
        long long lli = 0;
        int const n_elements_found = sscanf(h, "%s %lld", varname, &lli);
        if (n_elements_found < 2) {
            printf("# found only %d elements in line number %d\n", n_elements_found, linenumber);
            err = n_elements_found - 2;
        } else {
            uint64_t const hash = jr_simple_string_hash(varname);
            switch(hash) {
                     case 0x1505: // "", ignore
              break; case 0x1ae5ec19df29a9: // "JURASSIC"
              break; case 0xd0b773006c4b: err = (BINARY_TABLES_VERSION < lli); // "version"
                      if (echo) printf("# found header version %ld\n", lli);
              break; case 0xd0b2f9c4c031: // "git_key" (only for reference)
              break; case 0x310f71e19b:   err = sizeof(real_tblND_t) - 4; assert(4 == lli); // "float"  4
              break; case 0x652f93d5b20:  err = sizeof(real_tblND_t) - 8; assert(8 == lli); // "double" 8
              break; case 0x59749a:       err = NG - lli; // "NG"
              break; case 0x310e148925:   err = TBLNP - lli; // "TBLNP"
              break; case 0x310e148929:   err = TBLNT - lli; // "TBLNT"
              break; case 0x310e14892a:   err = TBLNU - lli; // "TBLNU"
              break; case 0x597497:       err = ND - lli; // "ND"
              break; case 0x310e148928:   err = TBLNS - lli; // "TBLNS"
              break; case 0x377c80eaed4807f: // "file_size"
              break; case 0xc09431708dbaaba8: header_size = lli; // "header_size"
                    if (lli != BINARY_TABLES_HEADER_LEN)
                        printf("# header sizes differ: %lld and %d\n", BINARY_TABLES_HEADER_LEN, lli);
              break; case 0x72730e3c0af490e7: table_size = lli; // "table_size"
                    err = (table_size < sizeof(tbl_t));
                    if (err) printf("# table sizes differ: %lld < %lld\n", table_size, sizeof(tbl_t));
              break; case 0x5978da: err = (lli < ctl->ng); // "ng"
              break; case 0x5978d7: err = (lli < ctl->nd); // "nd"
              break; case 0x18b2e43830995656: has_fast_inverse = lli; // "FAST_INVERSE_OF_U"
                      if (echo) printf("# found FAST_INVERSE_OF_U %ld\n", lli);
              break; case 0x727118c559a093e4: // "header_end"
              break;
              // maybe one out of 30 known gas names
              case 0xb87da49:     // "CO2"
              case 0xb87ebee:     // "H2O"
              case 0x597485:      // "N2"
              case 0x5974a6:      // "O2"
              case 0x5974a7:      // "O3"
              case 0x17c82ab14:   // "C2H2"
              case 0x17c82ab18:   // "C2H6"
              case 0x17c82f80b:   // "CCl4"
              case 0xb87d964:     // "CH4"
              case 0xb87de23:     // "ClO"
              case 0x652abf7a572: // "ClONO2"
              case 0x597337:      // "CO"
              case 0x17c83262f:   // "COF2"
              case 0xb87e32d:     // "F11"
              case 0xb87e32e:     // "F12"
              case 0xb87e330:     // "F14"
              case 0xb87e34f:     // "F22"
              case 0x17c8569e0:   // "H2O2"
              case 0xb87ee1e:     // "HCN"
              case 0x17c85e0fd:   // "HNO3"
              case 0x17c85e0fe:   // "HNO4"
              case 0x17c85e3eb:   // "HOCl"
              case 0xb880574:     // "N2O"
              case 0x17c88b429:   // "N2O5"
              case 0xb88082e:     // "NH3"
              case 0x5974a2:      // "NO"
              case 0xb880914:     // "NO2"
              case 0xb880bea:     // "OCS"
              case 0xb881d34:     // "SF6"
              case 0xb881e59:     // "SO2"
              { // treat varname as a gas name
                  int const ig = (int)lli; // assume that lli is the gas index
                  if (ig < 0) {                 err = ig;
                      printf("# try to interpret %s as gas name but found index %i out of range\n", varname, ig);
                  } else if (ig >= NG) {        err = ig;
                      printf("# try to interpret %s as gas name but found index %i out of range %d=NG\n", varname, ig, NG);
                  } else if (ig >= ctl->ng) { // this is NOT an error, we should be able to run with less gases
                      if (echo) printf("# interpret %s as gas name, found index %i larger than current setting ng=%d\n", varname, ig, ctl->ng);
                  } else {
                      char const *gas_ref = ctl->emitter[ig];
                      err = (hash != jr_simple_string_hash(gas_ref));
                      if (err) {
                          printf("# interpret %s as gas name but found gas %s at index %i\n", varname, gas_ref, ig);
                      } else {
                          gas_found[ig] = 1;
                      }
                  }
              }
              break; default: 
              { // treat as detector frequency channel
                  int const id = (int)lli; // assume that lli is the frequency index
                  if (id < 0) {                 err = id;
                      printf("# try to interpret %s as channel frequency but found index %i out of range\n", varname, id);
                  } else if (id >= ND) {        err = id;
                      printf("# try to interpret %s as channel frequency but found index %i out of range %d=ND\n", varname, id, ND);
                  } else if (id >= ctl->nd) { // this is NOT an error, we should be able to run with less channels
                     if (echo) printf("# interpret %s as channel frequency, found index %i larger than current setting nd=%d\n", varname, id, ctl->nd);
                  } else {
                      char nu_ref[16]; sprintf(nu_ref, "%.4f", ctl->nu[id]);
                      err = (hash != jr_simple_string_hash(nu_ref));
                      if (err) {
                          printf("# interpret %s as channel frequency but found %s at index %i\n", varname, nu_ref, id);
                      } else {
                          nu_found[id] = 1;
                      }
                  }
              }
            } // switch hash

        } // found less than 2 elements in line
            
        if (err) {
            printf("# found and error (code=%i) in linenumber %i related to entry \"%s %d\" checking the header\n", 
                    err, linenumber, varname, lli);
            ++nerrors;
        }

        while (*h != '\n') { ++h; } // forward to the next line start
        while (*h == '\n') { ++h; } // forward to the next char which is not a newline
        ++linenumber;

    } // while (*h)

    int n_gases_found = 0; for(int ig = 0; ig < ctl->ng; ++ig) n_gases_found += gas_found[ig];
    if (n_gases_found < ctl->ng) { ++nerrors; printf("# found only %d of %d gas species!\n", n_gases_found, ctl->ng); }
    int n_nu_found = 0; for(int id = 0; id < ctl->nd; ++id) n_nu_found += nu_found[id];
    if (n_nu_found < ctl->nd) { ++nerrors; printf("# found only %d of %d frequency channels!\n", n_nu_found, ctl->nd); }

    if (nerrors) {
        printf("# found %d error%s checking the header\n", nerrors, (1 == nerrors)?"":"s");
    } else if (echo) {
        printf("# found no errors checking the header\n", nerrors);
    }
    
    if (0 == nerrors && header_length) *header_length = header_size;
    return nerrors;
} // jr_binary_tables_check_header

int jr_binary_posix_write(tbl_t const *tbl, char const *filename, 
           char const *header, size_t const header_length) {
    printf("# try to open \"%s\" for writing binary\n", filename);
    FILE *fptr = fopen(filename, "wb");
    if (fptr == NULL) return __LINE__;
    printf("# opened \"%s\" for writing binary\n", filename);
    if (NULL == header) return __LINE__;
    printf("# writing header of length %.3f kByte binary\n", 1e-3*header_length);
    size_t const hl = fwrite(header, sizeof(char), header_length, fptr);
    printf("# header of length %.3f kByte written binary\n", 1e-3*header_length);
    if (hl != header_length) return __LINE__;
    if (NULL == tbl) return __LINE__;         

    printf("# writing tables with size %.3f MByte binary\n", 1e-6*sizeof(tbl_t));
    size_t const one = fwrite(tbl, sizeof(tbl_t), 1, fptr); // binary write of tbl
    printf("# tables with size %.3f MByte written binary\n", 1e-6*sizeof(tbl_t));

    fclose(fptr); 
    printf("# file \"%s\" written\n", filename);
    return one - 1;
} // jr_binary_posix_write

int jr_binary_posix_read(tbl_t *tbl, char const *filename, ctl_t const *ctl) {
    if (NULL == filename) return __LINE__; // error
    printf("# try to open \"%s\" for reading binary\n", filename);
    FILE *fptr = fopen(filename, "rb");
    if (NULL == fptr) return __LINE__;
    char header[BINARY_TABLES_HEADER_LEN];
    size_t const hl1 = BINARY_TABLES_HEADER_LEN;
    size_t const hl2 = fread(header, sizeof(char), hl1, fptr); // binary reading of header
    if (hl2 != hl1) return __LINE__; // error
    size_t hl3;
    int const status1 = jr_binary_tables_check_header(ctl, header, &hl3);
    if (status1) { // check failed
        if (hl3 > hl1) { // reason for failure was too short hl1
            char *h = malloc(hl3);
            if (NULL == h) return __LINE__;
            size_t const hl4 = fread(h, sizeof(char), hl3, fptr); // binary reading of header (again)
            if (hl4 != hl3) return __LINE__; // error
            size_t hl5;
            int const status2 = jr_binary_tables_check_header(ctl, h, &hl5);
            if (hl5 != hl4) return __LINE__;
            free(h);
            if (status2) return __LINE__; // still fails
        } else {
            return __LINE__; // some other failure reason
        }
    }
    if (NULL == tbl) { fclose(fptr); return __LINE__; } // error

    size_t const one = fread(tbl, sizeof(tbl_t), 1, fptr); // binary read of tbl

    fclose(fptr);
    printf("# file \"%s\" read, status = %lld\n", filename, one - 1);
    return one - 1;
} // jr_binary_posix_read

int jr_write_binary_tables(tbl_t const *tbl, ctl_t const *ctl) {
    TIMER("WRITE", 1); // start timer
    char filename[BINARY_TABLES_FILENAME_LEN];
    jr_binary_tables_filename(filename);
    char header[BINARY_TABLES_HEADER_LEN];
    if (jr_binary_tables_header(header, BINARY_TABLES_HEADER_LEN, ctl)) return __LINE__;
    int const status = jr_binary_posix_write(tbl, filename, header, BINARY_TABLES_HEADER_LEN);
    TIMER("WRITE", 3); // stop timer
    printf("# jr_write_binary_tables returns status %d\n", status);
    return status;
} // jr_write_binary_tables

int jr_read_binary_tables(tbl_t *tbl, ctl_t const *ctl) {
    TIMER("READ", 1); // start timer
    char filename[BINARY_TABLES_FILENAME_LEN];
    jr_binary_tables_filename(filename);
    int const status = jr_binary_posix_read(tbl, filename, ctl);
    TIMER("READ", 3); // stop timer
    printf("# jr_read_binary_tables returns status %d\n", status);
    return status;
} // jr_read_binary_tables

