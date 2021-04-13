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

//__asm__(".symver memcpy,memcpy@GLIBC_2.2.5");


using namespace std;

	#ifndef GALAXY
		#include "RcppArmadillo.h"
		using namespace Rcpp;
		#include <progress.hpp>
		// [[Rcpp::depends(RcppArmadillo)]]
		// [[Rcpp::depends(RcppProgress)]]
	#else
		#define Rcout cout
	#endif

#define DEF_lineLengthMax 10000
#define DEF_adaptLengthMax 500

#endif