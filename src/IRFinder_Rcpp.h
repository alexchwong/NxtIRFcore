#ifndef CODE_IRFINDER_RCPP
#define CODE_IRFINDER_RCPP

// Avoid including full Rcpp.h library as this increases compile time dramatically

// Also imports RcppProgress

// Declare GALAXY in make file so that these sources can be compiled into
// executable file from within Galaxy / Linux
#ifndef GALAXY
  #include <Rcpp.h>
  using namespace Rcpp;
  
  // [[Rcpp::depends(RcppProgress)]]
  #include <progress.hpp>
  
#else
  #define Rcout cout
  #define Rcpp::Rcout cout
#endif

#endif