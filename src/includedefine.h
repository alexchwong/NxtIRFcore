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

class chr_index {
  public:
    unsigned int refID;
    std::string chr_name;
    int32_t chr_len;
    chr_index(unsigned int a, std::string b, int32_t c) {
      refID = a;
      chr_name = b;
      chr_len = c;
    };
};

// Sort a vector of chr_index by chr_name in alphabetical order
inline bool operator< (const chr_index& lhs, const chr_index& rhs){
  return lhs.chr_name < rhs.chr_name;
}


#endif