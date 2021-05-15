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

//__asm__(".symver memcpy,memcpy@GLIBC_2.2.5");


using namespace std;

	#ifndef GALAXY
		#include <Rcpp.h>
		using namespace Rcpp;
		#include <progress.hpp>
		// [[Rcpp::depends(RcppProgress)]]
	#else
		#define Rcout cout
	#endif

#define DEF_lineLengthMax 10000
#define DEF_adaptLengthMax 500

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

// Sort a vector of chr_entry by chr_name in alphabetical order
inline bool operator< (const chr_entry& lhs, const chr_entry& rhs){
  return lhs.chr_name < rhs.chr_name;
}


#endif