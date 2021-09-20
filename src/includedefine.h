#ifndef INCLUDEDEFINE_DEF
#define INCLUDEDEFINE_DEF

// #pragma GCC diagnostic ignored "-Wpedantic"

#include <cstring>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <sys/types.h>
#include <vector>
#include <map>
#include <algorithm> // std::sort
#include <functional> // std::function
#include <stdexcept>
#include <zlib.h>
#include <zconf.h>
#include <math.h>
#include <chrono>
#include <thread>

using namespace std;

// Declare GALAXY in make file so that these sources can be compiled into
// executable file from within Galaxy / Linux

#ifndef GALAXY
  #include <Rcpp.h>
  using namespace Rcpp;
  #include <progress.hpp>
  // [[Rcpp::depends(RcppProgress)]]
  
  #include "ompBAM.hpp"
#else
  #define Rcout cout
  #define Rcpp::Rcout cout
  
  #include "ompBAM.hpp"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

// Tell Linux not to cache the large files:
//__asm__(".symver memcpy,memcpy@GLIBC_2.2.5");

#define DEF_lineLengthMax 10000
#define DEF_adaptLengthMax 500

// A class storing the refID, chrom name and length
class chr_entry {
  public:
    unsigned int refID;
    std::string chr_name;
    int32_t chr_len;
    
    chr_entry(unsigned int a, std::string b, int32_t c) {
      refID = a;
      chr_name = b;
      chr_len = c;
    };
};

// Sort a vector of chr_entry by chr_name in alphabetical order
inline bool operator< (const chr_entry& lhs, const chr_entry& rhs){
  return lhs.chr_name < rhs.chr_name;
}

#endif