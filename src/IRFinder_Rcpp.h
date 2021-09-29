#ifndef CODE_IRFINDER_RCPP
#define CODE_IRFINDER_RCPP

// Declare RNXTIRF in Makevars
// If not declared, this will produce an executable that does not depend on Rcpp
#ifdef RNXTIRF
  #include <Rcpp.h>
  using namespace Rcpp;
  
  // [[Rcpp::depends(RcppProgress)]]
  #include <progress.hpp>
  
  #define cout Rcpp::Rcout
#else
  #include <iostream>  
#endif

#endif