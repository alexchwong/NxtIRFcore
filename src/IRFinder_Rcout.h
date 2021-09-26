#ifndef CODE_IRFINDER_RCOUT
#define CODE_IRFINDER_RCOUT

// For Rcout
#ifndef GALAXY
  #include <RcppCommon.h>
  using namespace Rcpp;
#else
  #define Rcout cout
  #define Rcpp::Rcout cout
#endif

#endif