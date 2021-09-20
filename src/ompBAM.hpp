#ifndef _ompBAM
#define _ompBAM

#include <fstream>    // std::ifstream

// [[Rcpp::depends(zlibbioc)]]
#include <zlib.h>     // For BAM decompression
#include <zconf.h>

#include <string>    
#include <cstring>    // To compare between strings
#include <vector>     // For vector types

#ifdef _OPENMP
  #include <omp.h>    // For OpenMP
#endif

// #include "Rcpp.h"     // For Rcout
// The above is declared in includedefine.h

#include "pbam_defs.hpp"
#include "pbam1_t.hpp"
#include "pbam_in.hpp"

inline void ompBAM_version() {
  std::string version = "0.99.0";
  Rcout << "Compiled using ompBAM version " << version 
    << " for NxtIRFcore\n";
}

inline void is_compiled_with_OpenMP() {
  #ifdef _OPENMP
    Rcout << "Compiled with OpenMP; " << omp_get_max_threads()
      << " threads available\n";
  #else
    Rcout << "Compiled without OpenMP support\n";
  #endif
}

#endif

