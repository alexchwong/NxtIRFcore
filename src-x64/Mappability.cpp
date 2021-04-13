#include "Mappability.h"
#include <cassert>

std::string GenerateReadError(char * input_read, unsigned int read_len, unsigned int error_pos,
  unsigned int direction, unsigned int error_seed) {
  
  char * new_read = new char[read_len + 1];
  new_read[read_len] = '\0';
  if(direction == 0) {
    memcpy(&new_read[0], input_read, read_len);  
  } else {
    for(unsigned int i = 0; i < read_len; i++) {
      switch(input_read[i])
      {   
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
        new_read[read_len - i - 1] = input_read[i];
      }         
    }
  }

  char error_nuc = '\0';  // set this as something to avoid warning at compile
  switch(error_seed % 2) {
  case 0:
    switch(new_read[error_pos - 1]) {
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
    }
  case 1:
    switch(new_read[error_pos - 1]) {
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
    }    
  case 2:
    switch(new_read[error_pos - 1]) {
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
    }    
  }
  memcpy(&new_read[error_pos - 1], &error_nuc, 1);
  
  string return_str = string(new_read);
  delete[] new_read;
  return(return_str);
}

bool checkDNA(char * input_read, unsigned int read_len) {
  for(unsigned int i = 0; i < read_len; i++) {
    if(input_read[i]!='A' && input_read[i]!='T' && input_read[i]!='G' && input_read[i]!='C' &&
      input_read[i]!='a' && input_read[i]!='t' && input_read[i]!='g' && input_read[i]!='c') {
      return false;
    }
  }
  return true;
}

// [[Rcpp::export]]
int IRF_GenerateMappabilityReads(std::string genome_file, std::string out_fa,
	int read_len, int read_stride, int error_pos) {
  
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
  unsigned int seed = 0;
  
  string chr;
  string sequence;

  FastaReader inFA;
  inFA.SetInputHandle(&inGenome);
  
  while(!inGenome.eof() && !inGenome.fail()) {

    inFA.ReadSeq();
    sequence = inFA.sequence;
    chr = inFA.seqname;
    char * buffer = new char[sequence.length() + 1];
    std::strcpy (buffer, sequence.c_str());
    
    for(unsigned int bufferPos = 1; (bufferPos < sequence.length() - read_len - 1); bufferPos += read_stride) {
      memcpy(read, &buffer[bufferPos - 1], read_len);
      if(checkDNA(read, read_len)) {
        // std::string write_name;
        // write_name = (direction == 0 ? ">RF!" : ">RR!");
        // write_name.append(chr);
        // write_name.append("!");
        // write_name.append(std::to_string(bufferPos));

        
        std::string write_seq = GenerateReadError(read, read_len, error_pos, direction, seed) ;
        // outGZ.writeline(write_seq);
        if(is_stdout == 1) {
          // Rcout << write_name << '\n' << write_seq << '\n';
          Rcout << (direction == 0 ? ">RF!" : ">RR!") << chr << "!" << std::to_string(bufferPos)
            << '\n' << write_seq << '\n';  
        } else {
          // outFA << write_name << '\n' << write_seq << '\n';     
          outFA << (direction == 0 ? ">RF!" : ">RR!") << chr << "!" << std::to_string(bufferPos)
            << '\n' << write_seq << '\n'; 
        }
        seed += 1;
        direction = (direction == 0 ? 1 : 0);
      }
      if((seed % 100000 == 0) & (seed > 0)) {
        if(is_stdout == 0) {
          Rcout << "Processed " << bufferPos << " coord of chrom:" << chr << '\n';
        }
      }
    }
    delete[] buffer;
  }
  delete[] read;
  
  inGenome.close();
  if(is_stdout == 0) {
    outFA.flush();
    outFA.close();
  }
  return(0);
}

#ifndef GALAXY
// [[Rcpp::export]]
int IRF_GenerateMappabilityRegions(std::string bam_file, std::string output_file, int threshold, int includeCov){
  
  std::string s_output_txt = output_file + ".txt";
  std::string s_output_cov = output_file + ".cov";
#else
int IRF_GenerateMappabilityRegions(std::string bam_file, std::string s_output_txt, int threshold, std::string s_output_cov){	
#endif
  std::string s_inBAM = bam_file;
  
  FragmentsMap oFragMap;
  
  BAM2blocks BB;
  
  BB.registerCallbackChrMappingChange( std::bind(&FragmentsMap::ChrMapUpdate, &oFragMap, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&FragmentsMap::ProcessBlocks, &oFragMap, std::placeholders::_1) );
  
  BAMReader inbam;
  std::ifstream inbam_stream;
  if(bam_file == "-") {
    inbam.SetInputHandle(&std::cin);        
  } else {
    inbam_stream.open(s_inBAM, std::ifstream::binary);
    inbam.SetInputHandle(&inbam_stream);    
  }
  
  BB.openFile(&inbam);
  
  std::string BBreport;
  BB.processAll(BBreport);
  
  std::ofstream outFragsMap;
  outFragsMap.open(s_output_txt, std::ifstream::out);
  oFragMap.WriteOutput(&outFragsMap, BB.chr_names, BB.chr_lens, threshold);
  outFragsMap.flush(); outFragsMap.close();

#ifndef GALAXY
  if(includeCov == 1) {
#else
  if(!s_output_cov.empty()) {
#endif
    std::ofstream ofCOV;
    ofCOV.open(s_output_cov, std::ofstream::binary);
     
    covFile outCOV;
    outCOV.SetOutputHandle(&ofCOV);
    
    oFragMap.WriteBinary(&outCOV, BB.chr_names, BB.chr_lens);
    ofCOV.close();    
  }
  
  return(0);
}
