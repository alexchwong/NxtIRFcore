#include "ReadBlockProcessor.h" // No need to include GZWriter.h as this is included in ReadBlockProcessor.h
#include "BAM2blocks.h"
#include "FastaReader.h"
#include "includedefine.h"


std::string GenerateReadError(char * input_read, unsigned int read_len, unsigned int error_pos,
  unsigned int direction, unsigned int error_seed);

bool checkDNA(char * input_read, unsigned int read_len);

int IRF_GenerateMappabilityReads(std::string genome_file, std::string out_fa,
	int read_len, int read_stride, int error_pos);

#ifndef GALAXY
  int IRF_GenerateMappabilityRegions(std::string bam_file, std::string output_file, int threshold, int includeCov = 0);
#else
  int IRF_GenerateMappabilityRegions(std::string bam_file, std::string s_output_txt, int threshold, std::string s_output_cov = "");	
#endif
