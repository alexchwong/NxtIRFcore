#include "Mappability.h"
#include <cassert>

int Set_Threads(int n_threads) {
#ifdef _OPENMP
  int use_threads = 1;
	if(n_threads > 0 && n_threads <= omp_get_thread_limit()) {
    use_threads = n_threads;
	} else {
		use_threads = omp_get_thread_limit();
		if(use_threads < 1) {
			use_threads = 1;
		}
	}
	omp_set_num_threads(use_threads);
  return(use_threads);
#else
	return(1);
#endif
}

std::string GenerateReadError(
    char * input_read, 
    const unsigned int read_len, 
    const unsigned int error_pos,
    const unsigned int direction, 
    const unsigned int error_seed
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
  inFA.Profile();
  
  Progress p(inFA.chr_names.size(), (is_stdout == 0));
  
  while(!inGenome.eof() && !inGenome.fail()) {
    p.increment(1);
    inFA.ReadSeq();
    sequence = inFA.sequence;
    chr = inFA.seqname;
    char * buffer = new char[sequence.length() + 1];
    std::strcpy (buffer, sequence.c_str());

    for(unsigned int bufferPos = 1; (bufferPos < sequence.length() - read_len + 1); bufferPos += read_stride) {
      memcpy(read, &buffer[bufferPos - 1], read_len);
      seed += 1;
      if(checkDNA(read, read_len)) {       
        std::string write_seq = GenerateReadError(read, read_len, error_pos, direction, seed) ;

        if(is_stdout == 1) {
          Rcout << (direction == 0 ? ">RF!" : ">RR!") << chr << "!" << std::to_string(bufferPos)
            << '\n' << write_seq << '\n';  
        } else {
          outFA << (direction == 0 ? ">RF!" : ">RR!") << chr << "!" << std::to_string(bufferPos)
            << '\n' << write_seq << '\n'; 
        }
        direction = (direction == 0 ? 1 : 0);        
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
int IRF_GenerateMappabilityRegions(std::string bam_file, std::string output_file, int threshold, int includeCov, bool verbose, int n_threads = 1){
  
  std::string s_output_txt = output_file + ".txt";
  std::string s_output_cov = output_file + ".cov";
#else
int IRF_GenerateMappabilityRegions(std::string bam_file, std::string s_output_txt, int threshold, std::string s_output_cov){	
	bool verbose = true;
#endif

  // int use_threads = Set_Threads(n_threads);

  unsigned int n_threads_to_use = (unsigned int)n_threads;   // Should be sorted out in calling function
 
  std::string myLine;
	if(verbose) Rcout << "Calculating Mappability Exclusions from aligned synthetic reads in BAM file " << bam_file << "\n";
  
  
  // std::ifstream inbam_stream;   inbam_stream.open(bam_file, std::ios::in | std::ios::binary);
  pbam_in inbam((size_t)1000000000, (size_t)2000000000, 5);
  // inbam.SetInputHandle(&inbam_stream, n_threads_to_use);
    inbam.openFile(bam_file, n_threads_to_use);

  // Assign children:
  std::vector<FragmentsMap*> oFM;
  std::vector<BAM2blocks*> BBchild;

  for(unsigned int i = 0; i < n_threads_to_use; i++) {
    oFM.push_back(new FragmentsMap);
    BBchild.push_back(new BAM2blocks);

    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&FragmentsMap::ChrMapUpdate, &(*oFM.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&FragmentsMap::ProcessBlocks, &(*oFM.at(i)), std::placeholders::_1) );

    BBchild.at(i)->openFile(&inbam);
  }
  
  // BAM processing loop
  Progress p(inbam.GetFileSize(), verbose);
  // Rcout << "Total blocks: " << n_bgzf_blocks << '\n';
  // unsigned int blocks_read_total = 0;
  // int ret = 0;

  while(0 == inbam.fillReads() && !p.check_abort()) {
    p.increment(inbam.IncProgress());
    
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads_to_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      BBchild.at(i)->processAll(i);
    }


  }

  if(p.check_abort()) {
    // interrupted:
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      delete oFM.at(i);
      delete BBchild.at(i);
    }
    return(-1);
  }

  // inbam_stream.close();
  // Rcout << "BAM processing finished\n";
  
  if(n_threads_to_use > 1) {
    if(verbose) Rcout << "Compiling data from threads\n";
  // Combine BB's and process spares
    for(unsigned int i = 1; i < n_threads_to_use; i++) {
      BBchild.at(0)->processSpares(*BBchild.at(i));
      delete BBchild.at(i);
    }
  // Combine objects:
    for(unsigned int i = 1; i < n_threads_to_use; i++) {
      oFM.at(0)->Combine(*oFM.at(i));

      delete oFM.at(i);
    }
  }

#ifndef GALAXY
  if(includeCov == 1) {
#else
  if(!s_output_cov.empty()) {
#endif
   // Write Coverage Binary file:
    std::ofstream ofCOV;                          ofCOV.open(s_output_cov, std::ofstream::binary);  
    covFile outCOV;                               outCOV.SetOutputHandle(&ofCOV);
    oFM.at(0)->WriteBinary(&outCOV, verbose);     ofCOV.close();
  }
  
  std::ofstream outFragsMap;
  outFragsMap.open(s_output_txt, std::ifstream::out);
	
  oFM.at(0)->WriteOutput(&outFragsMap, threshold, verbose);
  outFragsMap.flush(); outFragsMap.close();

  // destroy objects:
  // for(unsigned int i = 0; i < n_threads_to_use; i++) {
    delete oFM.at(0);
    delete BBchild.at(0);
  // }

  return(0);
}
