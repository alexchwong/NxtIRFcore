#ifndef _CODE_COVSTREAMS
#define _CODE_COVSTREAMS

#include "includedefine.h"
#include "IRFinder_Rcpp.h"
// #include "IRFinder_Rcout.h" // For Rcout

#include <zlib.h>
#include <zconf.h>

#include "pbam_defs.hpp"

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