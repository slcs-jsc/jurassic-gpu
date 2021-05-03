#pragma once

#include <stdint.h> // uint64_t

  extern inline
  uint64_t jr_simple_string_hash(char const *string) {
    // this simple hash function maps only 9 out of 216000 english words onto the same integer number
      uint64_t hash = 5381; // 0x1505
      char const *c = string;
      while (*c != 0) {
         hash = hash * 33 + (*c);
         ++c;
      } // while
      return hash;
  } // simple_string_hash
