/* includedefine.h Shared includes and declarations for IRFinder/NxtIRF

Copyright (C) 2021 Alex Chit Hei Wong
Copyright (C) 2016 William Ritchie
  - original: https://github.com/williamritchie/IRFinder/tree/IRFinder-1.3.1)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.  */

#ifndef INCLUDEDEFINE_DEF
#define INCLUDEDEFINE_DEF

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
#include <math.h>
#include <chrono>
#include <thread>

// [[Rcpp::depends(zlibbioc)]]
#include <zlib.h>
#include <zconf.h>

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

union stream_uint64 {
  char c[8];
  uint64_t u;
};
union stream_uint32 {
  char c[4];
  uint32_t u;
};
union stream_int32 {
  char c[4];
  int32_t i;
};
union stream_uint16 {
  char c[2];
  uint16_t u;
};

#endif