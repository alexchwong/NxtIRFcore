/* Mappability.cpp Generates synthetic reads; maps low-mappability regions

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

#include "Mappability.h"

std::string GenerateReadError(
    char * input_read, 
    const unsigned int read_len, 
    const unsigned int error_pos,
    const unsigned int direction, 
    const size_t error_seed
) {
  
  // Copy https://github.com/williamritchie/IRFinder/blob/master/bin/util/generateReadsError.pl

  char * new_read_inter = new char[read_len + 1];
  new_read_inter[read_len] = '\0';
  memcpy(&new_read_inter[0], input_read, read_len);  

  char error_nuc = '\0';  // set this as something to avoid warning at compile
  if(error_seed % 3 == 0) {
    switch(new_read_inter[error_pos - 1]) {
      case 'A':
        error_nuc = 'G'; break;
      case 'C':
        error_nuc = 'A'; break;
      case 'G':
        error_nuc = 'T'; break;
      case 'T':
        error_nuc = 'C'; break;
      case 'a':
        error_nuc = 'g'; break;
      case 'c':
        error_nuc = 'a'; break;
      case 'g':
        error_nuc = 't'; break;
      case 't':
        error_nuc = 'c'; break;
      default:
        error_nuc = 'N';
    }
  } else if(error_seed % 3 == 1) {
    switch(new_read_inter[error_pos - 1]) {
      case 'A':
        error_nuc = 'T'; break;
      case 'C':
        error_nuc = 'G'; break;
      case 'G':
        error_nuc = 'C'; break;
      case 'T':
        error_nuc = 'A'; break;
      case 'a':
        error_nuc = 't'; break;
      case 'c':
        error_nuc = 'g'; break;
      case 'g':
        error_nuc = 'c'; break;
      case 't':
        error_nuc = 'a'; break;
      default:
        error_nuc = 'N';
    }
  } else {
    switch(new_read_inter[error_pos - 1]) {
      case 'A':
        error_nuc = 'C'; break;
      case 'C':
        error_nuc = 'T'; break;
      case 'G':
        error_nuc = 'A'; break;
      case 'T':
        error_nuc = 'G'; break;
      case 'a':
        error_nuc = 'c'; break;
      case 'c':
        error_nuc = 't'; break;
      case 'g':
        error_nuc = 'a'; break;
      case 't':
        error_nuc = 'g'; break;
      default:
        error_nuc = 'N';
    }
  }
  
  memcpy(&new_read_inter[error_pos - 1], &error_nuc, 1);
  
  char * new_read = new char[read_len + 1];
  new_read[read_len] = '\0';
  if(direction == 0) {
    memcpy(&new_read[0], new_read_inter, read_len);  
  } else {
    for(unsigned int i = 0; i < read_len; i++) {
      switch(new_read_inter[i]) {   
        case 'A':
          new_read[read_len - i - 1] = 'T'; break;
        case 'T':
          new_read[read_len - i - 1] = 'A'; break;
        case 'G':
          new_read[read_len - i - 1] = 'C'; break;
        case 'C':
          new_read[read_len - i - 1] = 'G'; break;
        case 'a':
          new_read[read_len - i - 1] = 't'; break;
        case 't':
          new_read[read_len - i - 1] = 'a'; break;
        case 'g':
          new_read[read_len - i - 1] = 'c'; break;
        case 'c':
          new_read[read_len - i - 1] = 'g'; break;
        default :
          new_read[read_len - i - 1] = 'N';
      }         
    }
  }

  string return_str = string(new_read);
  delete[] new_read;
  return(return_str);
}

// Replicate old PERL script; return true if N's constitute less than half of length
bool checkDNA(char * input_read, unsigned int read_len) {
  unsigned int numN = 0;
  for(unsigned int i = 0; i < read_len; i++) {
    if(input_read[i]!='A' && input_read[i]!='T' && input_read[i]!='G' && input_read[i]!='C' &&
      input_read[i]!='a' && input_read[i]!='t' && input_read[i]!='g' && input_read[i]!='c') {
      numN++;
    }
  }
  return(numN < read_len / 2);
}

// [[Rcpp::export]]
int IRF_GenerateMappabilityReads(
  std::string genome_file, std::string out_fa,
	int read_len, int read_stride, int error_pos
) {
  
  std::ifstream inGenome;
  inGenome.open(genome_file, std::ifstream::in);
  
  // Allows writing to standard output if filename is '-'
  int is_stdout = 0;
  std::ofstream outFA;
  if(out_fa == "-") {
    is_stdout = 1;
  } else {
    outFA.open(out_fa, std::ios::binary);    
  }
      
  unsigned int direction = 0;
  char * read = new char[read_len + 1];
  size_t num_reads = 0;
  
  string chr;
  string sequence;

  FastaReader inFA;
  inFA.SetInputHandle(&inGenome);
  inFA.Profile();

#ifndef GALAXY
  Progress p(inFA.total_size, (is_stdout == 0));
#endif

  while(!inGenome.eof() && !inGenome.fail()) {

    inFA.ReadSeq();
    sequence = inFA.sequence;
    chr = inFA.seqname;
    char * buffer = new char[sequence.length() + 1];
    std::strcpy (buffer, sequence.c_str());

    size_t seq_progress = 0;
    for(
        unsigned int bufferPos = 1; 
        bufferPos < sequence.length() - read_len + 1; 
        bufferPos += read_stride
    ) {
      memcpy(read, &buffer[bufferPos - 1], read_len);
      num_reads += 1;
      if(checkDNA(read, read_len)) {       
        std::string write_seq = GenerateReadError(
          read, read_len, error_pos, direction, num_reads
        ) ;

        if(is_stdout == 1) {
          Rcout << (direction == 0 ? ">RF!" : ">RR!") << chr << "!" << std::to_string(bufferPos)
            << '\n' << write_seq << '\n';  
        } else {
          outFA << (direction == 0 ? ">RF!" : ">RR!") << chr << "!" << std::to_string(bufferPos)
            << '\n' << write_seq << '\n'; 
        }
        direction = (direction == 0 ? 1 : 0);
        
        // update progress bar
#ifndef GALAXY
        if(num_reads % 100000 == 0) {
          p.increment(bufferPos - seq_progress);
          seq_progress = bufferPos;
        }
#endif
      }
    }
#ifndef GALAXY
    p.increment(sequence.length() - seq_progress);
#endif
    delete[] buffer;
  }
  delete[] read;
  
  inGenome.close();
  if(is_stdout == 0) {
    outFA.flush();
    outFA.close();
  }
  Rcout << num_reads << " synthetic reads generated\n";
  return(0);
}
