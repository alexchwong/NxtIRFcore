#ifndef CODE_IRFINDER
#define CODE_IRFINDER

#include "includedefine.h"

// Full Rcpp functionality
#include "IRFinder_Rcpp.h"

#include "BAM2blocks.h"

#include "covReader.h"
#include "covWriter.h"

#include "ReadBlockProcessor.h"
#include "ReadBlockProcessor_CoverageBlocks.h"  // includes FragmentsMap

#include "GZReader.h"
#include "GZWriter.h"

int Has_OpenMP();
int Set_Threads(int n_threads);

bool IRF_Check_Cov(std::string s_in);

#ifndef GALAXY
// Rcpp-only functions

  List IRF_RLE_From_Cov(
    std::string s_in, std::string seqname, int start, int end, int strand
  );

  List IRF_RLEList_From_Cov(std::string s_in, int strand);

  List IRF_gunzip_DF(std::string s_in, StringVector s_header_begin);
  
#endif

int IRF_gunzip(std::string s_in, std::string s_out);

int ReadChrAlias(std::istringstream &IN,
    std::vector<std::string> &ref_names, 
    std::vector<std::string> &ref_alias,
    std::vector<uint32_t> &ref_lengths
);

int IRF_ref(std::string &reference_file, 
    std::vector<std::string> &ref_names, 
    std::vector<std::string> &ref_alias,
    std::vector<uint32_t> &ref_lengths,
    CoverageBlocksIRFinder &CB_template, 
    SpansPoint &SP_template, 
    FragmentsInROI &ROI_template,
    JunctionCount &JC_template, 
    bool verbose
);

int IRF_core(std::string const &bam_file, 
    std::string const &s_output_txt, std::string const &s_output_cov,
    std::vector<std::string> &ref_names, 
    std::vector<std::string> &ref_alias,
    std::vector<uint32_t> &ref_lengths,
    CoverageBlocksIRFinder const &CB_template, 
    SpansPoint const &SP_template, 
    FragmentsInROI const &ROI_template,
    JunctionCount const &JC_template,
    bool const verbose,
    int n_threads = 1
);

#ifndef GALAXY
  int IRF_main(
      std::string bam_file, std::string reference_file, std::string output_file, 
      bool verbose = true, int n_threads = 1
  );

  int IRF_main_multi(
      std::string reference_file, StringVector bam_files, StringVector output_files,
      int max_threads = 1, bool verbose = true
  );

  int IRF_GenerateMappabilityRegions(
      std::string bam_file, std::string output_file, 
      int threshold, int includeCov, bool verbose,
      int n_threads
  );

#else
  int IRF_main(
      std::string bam_file, std::string reference_file, std::string s_output_txt,
      std::string s_output_cov, int n_threads = 1
  );

  int IRF_GenerateMappabilityRegions(
      std::string bam_file, std::string s_output_txt, 
      int threshold, std::string s_output_cov
  );

  int main(int argc, char * argv[]);
#endif

#endif