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
  inFA.Profile();
  
  Progress p(inFA.chr_names.size(), (is_stdout == 0));
  
  while(!inGenome.eof() && !inGenome.fail()) {
    p.increment(1);
    inFA.ReadSeq();
    sequence = inFA.sequence;
    chr = inFA.seqname;
    char * buffer = new char[sequence.length() + 1];
    std::strcpy (buffer, sequence.c_str());
		// Rcout << "Processing chromosome: " << chr << '\n';
		// Progress p(sequence.length(), (is_stdout == 0));

		// unsigned int bufferPos_prev = 1;
    for(unsigned int bufferPos = 1; (bufferPos < sequence.length() - read_len - 1); bufferPos += read_stride) {
      memcpy(read, &buffer[bufferPos - 1], read_len);
      if(checkDNA(read, read_len)) {       
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
      // if((seed % 100000 == 0) & (seed > 0)) {
        // if(is_stdout == 0) {
          // Rcout << "Processed " << bufferPos << " coord of chrom:" << chr << '\n';
					// p.increment((unsigned long)(bufferPos - bufferPos_prev));	
					// bufferPos_prev = bufferPos;
				// }
      // }
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

  int use_threads = Set_Threads(n_threads);

  unsigned int n_threads_to_use = (unsigned int)n_threads;   // Should be sorted out in calling function
 
  std::string myLine;
	if(verbose) Rcout << "Calculating Mappability Exclusions from aligned synthetic reads in BAM file " << bam_file << "\n";
  
  
  std::ifstream inbam_stream;   inbam_stream.open(bam_file, std::ios::in | std::ios::binary);
  BAMReader_Multi inbam;        inbam.SetInputHandle(&inbam_stream); // Rcout << "BAMReader_Multi handle set\n";  
  
  BAM2blocks BB;  
  if(verbose) Rcout << "Identifying BGZF blocks in BAM file\n";
  unsigned int n_bgzf_blocks = BB.openFile(&inbam, verbose, n_threads_to_use);
  // This step writes chrs to BB, and BB obtains bgzf block positions for each worker
  if(n_bgzf_blocks == 0) {
    Rcout << "Error occurred profiling BAM file\n";
    return(-1);
  }
  // Assign children:
  std::vector<FragmentsMap*> oFM;
  std::vector<BAM2blocks*> BBchild;
  std::vector<BAMReader_Multi*> BRchild;

  for(unsigned int i = 0; i < n_threads_to_use; i++) {
    oFM.push_back(new FragmentsMap);
    BBchild.push_back(new BAM2blocks);
    BRchild.push_back(new BAMReader_Multi);

    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&FragmentsMap::ChrMapUpdate, &(*oFM.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&FragmentsMap::ProcessBlocks, &(*oFM.at(i)), std::placeholders::_1) );

    // Assign task:
    uint64_t begin_bgzf; unsigned int begin_pos;
    uint64_t end_bgzf; unsigned int end_pos;
    BB.ProvideTask(i, begin_bgzf, begin_pos, end_bgzf, end_pos);
    
    BRchild.at(i)->AssignTask(&inbam_stream, begin_bgzf, begin_pos, end_bgzf, end_pos);
    BRchild.at(i)->SetAutoLoad(false);
    
    BBchild.at(i)->AttachReader(BRchild.at(i));
    BBchild.at(i)->TransferChrs(BB);
  }
  
  // BAM processing loop
  Progress p(n_bgzf_blocks, verbose);
  // Rcout << "Total blocks: " << n_bgzf_blocks << '\n';
  unsigned int blocks_read_total = 0;
  int ret = 0;
#ifdef _OPENMP
  #pragma omp parallel for
  for(unsigned int i = 0; i < n_threads_to_use; i++) {
    unsigned int n_blocks_read = 1;
    int ret2 = 0;
    while(!BRchild.at(i)->eob() && !p.check_abort() && n_blocks_read > 0 && ret == 0) {
      #pragma omp critical
      n_blocks_read = (unsigned int)BRchild.at(i)->read_from_file(100);
      
      if(n_blocks_read > 0) {
        BRchild.at(i)->decompress();
        ret2 = BBchild.at(i)->processAll();
        
        if(ret2 == -1) {
          #pragma omp critical
          ret = -1;    // abort if broken reads detected
        }
        
        #pragma omp atomic
        blocks_read_total += n_blocks_read;
        
        #pragma omp critical
        p.increment(n_blocks_read);
      }
      // Rcout << "Blocks read: " << n_blocks_read << '\n';
    }
  }
#else
  for(unsigned int i = 0; i < n_threads_to_use; i++) {
    while(!BRchild.at(i)->eob()) {
      
      int n_blocks_read = BRchild.at(i)->read_from_file(100);
      BRchild.at(i)->decompress();
      BBchild.at(i)->processAll();
      
      p.increment(n_blocks_read);
    }
  }
#endif
  if(p.check_abort()) {
    // interrupted:
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      delete oFM.at(i);
      delete BRchild.at(i);
      delete BBchild.at(i);
    }
    return(-1);
  }

  int final_inc = max((int)(n_bgzf_blocks - blocks_read_total), 1);
  p.increment(final_inc);

  inbam_stream.close();
  // Rcout << "BAM processing finished\n";
  
  if(n_threads_to_use > 1) {
    if(verbose) Rcout << "Compiling data from threads\n";
  // Combine BB's and process spares
    for(unsigned int i = 1; i < n_threads_to_use; i++) {
      BBchild.at(0)->processSpares(*BBchild.at(i));
      delete BBchild.at(i);
      delete BRchild.at(i);
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
    delete BRchild.at(0);
    delete BBchild.at(0);
  // }

  return(0);
}
