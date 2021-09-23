/* Mappability.h Generates synthetic reads; maps low-mappability regions

Copyright (C) 2021 Alex Chit Hei Wong

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.  */

#ifndef CODE_MAPPABILITY
#define CODE_MAPPABILITY

#include "ReadBlockProcessor_FragmentsMap.h"
#include "GZWriter.h"
#include "BAM2blocks.h"
#include "FastaReader.h"
#include "includedefine.h"

int Set_Threads(int n_threads);

std::string GenerateReadError(
    char * input_read, 
    const unsigned int read_len, 
    const unsigned int error_pos,
    const unsigned int direction, 
    const unsigned int error_seed
);

bool checkDNA(char * input_read, unsigned int read_len);

int IRF_GenerateMappabilityReads(
  std::string genome_file, std::string out_fa,
	int read_len, int read_stride, int error_pos
);

#ifndef GALAXY
int IRF_GenerateMappabilityRegions(
  std::string bam_file, std::string output_file, 
  int threshold, int includeCov = 0, bool verbose = true,
  int n_threads = 1
);
#else
int IRF_GenerateMappabilityRegions(
  std::string bam_file, std::string s_output_txt, 
  int threshold, std::string s_output_cov = ""
);	
#endif

#endif