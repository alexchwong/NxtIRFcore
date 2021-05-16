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

#ifdef _OPENMP
#include <omp.h>
#endif

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

union bam_header {
  char c[8];
  struct {
    char magic[4];
    int32_t l_text;
  } magic;
};

// TODO -- are structs best hidden inside the class? Does doing so push them into namespace of the class only?
struct bam_read_core {
  union {
    char c_block_size[4];
    uint32_t block_size;
  };    
  union {
    char c[32];
    struct {
    // uint32_t block_size;
    int32_t refID;
    int32_t pos;
    uint8_t l_read_name;
    uint8_t mapq;
    uint16_t bin;
    uint16_t n_cigar_op;
    uint16_t flag;
    uint32_t l_seq;
    int32_t next_refID;
    int32_t next_pos;
    int32_t tlen;
    } core; // anonymous struct is now named.
  };
  char read_name[256];
  union {
    char cigar_buffer[2000];
    uint32_t cigar[500];
  };
};



#endif