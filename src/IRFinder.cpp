#include "includedefine.h"

#include "ReadBlockProcessor.h"
#include "ReadBlockProcessor_CoverageBlocks.h"
#include "BAM2blocks.h"
#include "GZReader.h"
#include "Mappability.h"

const char refEOF[5] =
		"\x20\x45\x4f\x46";

#ifndef GALAXY

    // [[Rcpp::export]]
    int Has_OpenMP() {
    #ifdef _OPENMP
      return omp_get_max_threads();
    #else
      return 0;
    #endif
    }

// [[Rcpp::export]]
bool IRF_Check_Cov(std::string s_in) {
	// Checks if given file is a valid COV file
	
  std::ifstream inCov_stream;
  inCov_stream.open(s_in, std::ifstream::binary);	

  covFile inCov;
  inCov.SetInputHandle(&inCov_stream);

  if(inCov.fail()){
		inCov_stream.close();
    return(false);
  }

  int ret = inCov.ReadHeader();
  if(ret == -1){
		inCov_stream.close();	
    return(false);
  }	
	
  inCov_stream.close();	
	return(true);
}

// [[Rcpp::export]]
List IRF_RLE_From_Cov(std::string s_in, std::string seqname, int start, int end, int strand) {
// Returns an RLE covering the region described above
// s_in: The coverage file
// strand: 0 = -, 1 = +, 2 = *
  
  List NULL_RLE = List::create(
    _["values"] = 0,
    _["lengths"] = 0 
  );
  
  if(start > end || start < 0){
    return(NULL_RLE);
  }

  std::ifstream inCov_stream;
  inCov_stream.open(s_in, std::ifstream::binary);
  
  covFile inCov;
  inCov.SetInputHandle(&inCov_stream);

  if(inCov.fail()){
		inCov_stream.close();
    return(NULL_RLE);
  }
	
  int ret = inCov.ReadHeader();
  if(ret < 0){
		Rcout << s_in << " appears to not be valid COV file... exiting\n";
		inCov_stream.close();	
    return(NULL_RLE);
  }
  
  // Find corresponding seqname
  unsigned int ref_index = 0;
  std::vector<chr_entry> chrs;
  inCov.GetChrs(chrs);
  while(0 != seqname.compare(0, seqname.size(), chrs.at(ref_index).chr_name)) {
    if(ref_index == chrs.size()) break;
    ref_index++;
  }
  if(ref_index == chrs.size()) {
    inCov_stream.close();	
    return(NULL_RLE);
  }
  
  // auto it_chr = std::find(inCov.chr_names.begin(), inCov.chr_names.end(), seqname);
  // if(it_chr == inCov.chr_names.end()) {
		// inCov_stream.close();	
    // return(NULL_RLE);
  // } else {
    // ref_index = distance(inCov.chr_names.begin(), it_chr);
  // }
  // end = 0 implies fetch whole chromosome
  int eff_end = 0;
  if(end == 0) {
    eff_end = chrs.at(ref_index).chr_len;
  } else {
    eff_end = end;
  }
  
  std::vector<int> values;
  std::vector<unsigned int> lengths;
  // Push first value
  values.push_back(0);
  lengths.push_back((unsigned int)start);
  
  inCov.FetchRLE(seqname, (uint32_t)start, (uint32_t)eff_end, strand, &values, &lengths);

  inCov_stream.close();
  // Push last value
  if((uint32_t)eff_end < (uint32_t)chrs.at(ref_index).chr_len) {
    values.push_back(0);
    lengths.push_back((uint32_t)chrs.at(ref_index).chr_len - eff_end);
  }
    
  List RLE = List::create(
    _["values"] = values,
    _["lengths"] = lengths 
  );

  return(RLE);
}

// [[Rcpp::export]]
List IRF_RLEList_From_Cov(std::string s_in, int strand) {
  // Returns an RLEList
  // s_in: The coverage file
  // strand: 0 = -, 1 = +, 2 = *
  
  List NULL_RLE = List::create(
    _["values"] = 0,
    _["lengths"] = 0 
  );
  
  List RLEList;
  
  std::ifstream inCov_stream;
  inCov_stream.open(s_in, std::ifstream::binary);
  
  covFile inCov;
  inCov.SetInputHandle(&inCov_stream);

  if(inCov.fail()){
		inCov_stream.close();
    return(NULL_RLE);
  }
  
  int ret = inCov.ReadHeader();
  if(ret == -1){
		Rcout << s_in << " appears to not be valid COV file... exiting";
		inCov_stream.close();	
    return(NULL_RLE);
  }
  
  std::vector<chr_entry> chrs;
  inCov.GetChrs(chrs);

  for (unsigned int i = 0; i < chrs.size(); i++) {
    uint32_t eff_end = chrs.at(i).chr_len;
    uint32_t start = 0;
    
    std::vector<int> values;
    std::vector<unsigned int> lengths;

    inCov.FetchRLE(chrs.at(i).chr_name, (uint32_t)start, (uint32_t)eff_end, strand, &values, &lengths);
    
    List RLE = List::create(
      _["values"] = values,
      _["lengths"] = lengths 
    );
    RLEList.push_back(RLE, chrs.at(i).chr_name);
  }

  inCov_stream.close();
  
  return(RLEList);
}

// [[Rcpp::export]]
int IRF_gunzip(std::string s_in, std::string s_out) {
  
  GZReader gz_in;
  int ret = gz_in.LoadGZ(s_in, true);
  if(ret != 0) return(ret);
	
  std::ofstream out;
  out.open(s_out, std::ofstream::binary);
  std::string myLine;
  
  while(!gz_in.iss.eof()) {
    getline(gz_in.iss, myLine, '\n');
    out << myLine << "\n";
  }
  out.flush(); out.close();
  
  return(0);
}

// [[Rcpp::export]]
List IRF_gunzip_DF(std::string s_in, StringVector s_header_begin) {
  List Final_final_list;
  
  GZReader gz_in;
  int ret = gz_in.LoadGZ(s_in, false, true);
  if(ret != 0) return(Final_final_list);
	
  // std::ofstream out;
  // out.open(s_out, std::ifstream::out);
  
  // Look for first line of data to return
  std::string myLine;
  std::string myEntry;
  unsigned int q = 0;
  char delim = '\n';
	bool check_line = true;
	
  for(int z = 0; z < s_header_begin.size(); z++) {
    std::string header = string(s_header_begin(z));
    std::vector<std::string> columns;
    while(!gz_in.eof()) {
      gz_in.getline(myLine, delim); q++;
      myLine.erase( std::remove(myLine.begin(), myLine.end(), '\r'), myLine.end() ); // remove \r 
      
			if(check_line == true) {
				if(strncmp(myLine.c_str(), header.c_str(), header.size()) == 0) {
					// read columns
					std::istringstream column_iss;
					column_iss.str(myLine);
					
					while(!column_iss.eof() && !column_iss.fail()) {
						getline(column_iss, myEntry, '\t');
						columns.push_back(myEntry);
					}
					
					break;
				} else {
					check_line = false;
				}
			} else {
				// screen for empty lines - then reactivate check_line
				if (myLine.length() == 0) {
					check_line = true;
				}
			}

    }

    // use a map of string vectors
    std::map< std::string, std::vector<std::string> > column_data;
      
    while(!gz_in.eof()) {
      gz_in.getline(myLine, delim); q++;
      myLine.erase( std::remove(myLine.begin(), myLine.end(), '\r'), myLine.end() ); // remove \r 
      if (myLine.length() == 0) {
        break;  // End at detection of empty line
      }
      std::istringstream entry_iss;
      entry_iss.str(myLine);
      unsigned int j = 0;
      while(!entry_iss.eof() && !entry_iss.fail() && j < columns.size()) {
        getline(entry_iss, myEntry, '\t');
        column_data[columns.at(j)].push_back(myEntry);
        j++;
      }
      if(j > columns.size()) {
        Rcout << "Detecting extra rows at line" << q << '\n';
        // ignore for now
      } else if(j != columns.size()) {
        Rcout << "Missing entries detected at line" << q << '\n';
        // attempt to repair by putting blank entries
        for(unsigned int k = j; k < columns.size(); k++) {
          column_data[columns.at(k)].push_back("");
        }
      }
    }
    List final_list;
    for(unsigned int i = 0; i < columns.size(); i++) {
      final_list.push_back(column_data[columns.at(i)], columns.at(i));
    }
    Final_final_list.push_back(final_list, header);
  }
	gz_in.closeGZ();
  return(Final_final_list);
}

#else
	// galaxy
#endif

// IRFinder reference reader:
int IRF_ref(std::string &reference_file, 
    CoverageBlocksIRFinder &CB_template, 
    SpansPoint &SP_template, 
    FragmentsInROI &ROI_template,
    JunctionCount &JC_template, 
    bool verbose
) { 
  GZReader * gz_in = new GZReader;
  int ret = gz_in->LoadGZ(reference_file, true);
  if(ret != 0) return(-1);
  
  std::string myLine;
  std::string myBuffer;
  
  getline(gz_in->iss, myLine, '#');    // discard first >
  getline(gz_in->iss, myLine, '\n');   // ignore file names for now
  getline(gz_in->iss, myBuffer, '#');  // this is the data block for ref-cover.bed

  std::istringstream inCoverageBlocks;
  inCoverageBlocks.str(myBuffer);
  CB_template.loadRef(inCoverageBlocks);
  
    getline(gz_in->iss, myLine, '\n');
    getline(gz_in->iss, myBuffer, '#');

  SP_template.setSpanLength(5,4);
  std::istringstream inSpansPoint;
  inSpansPoint.str(myBuffer);
  SP_template.loadRef(inSpansPoint);

  getline(gz_in->iss, myLine, '\n');
  getline(gz_in->iss, myBuffer, '#');

  std::istringstream inFragmentsInROI;
  inFragmentsInROI.str(myBuffer);
  ROI_template.loadRef(inFragmentsInROI);

  getline(gz_in->iss, myLine, '\n');
  getline(gz_in->iss, myBuffer, '#');

  std::istringstream inJuncCount;
  inJuncCount.str(myBuffer);
  JC_template.loadRef(inJuncCount);

// Ensure valid reference termination:
  getline(gz_in->iss, myLine, '\n');
  delete gz_in;
  
  if(strncmp(myLine.c_str(), refEOF, 4) != 0) {
    Rcout << "Invalid IRFinder reference detected\n";
    return(-1);
  }
  return(0);
}

// IRFinder core:
int IRF_core(std::string const &bam_file, 
    std::string const &s_output_txt, std::string const &s_output_cov,
    CoverageBlocksIRFinder const &CB_template, 
    SpansPoint const &SP_template, 
    FragmentsInROI const &ROI_template,
    JunctionCount const &JC_template,
    bool const verbose,
    int n_threads = 1
) {
  unsigned int n_threads_to_use = (unsigned int)n_threads;   // Should be sorted out in calling function
 
  std::string myLine;
	if(verbose) Rcout << "Processing BAM file " << bam_file << "\n";
  
  
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
  std::vector<CoverageBlocksIRFinder*> oCB;
  std::vector<SpansPoint*> oSP;
  std::vector<FragmentsInROI*> oROI;
  std::vector<FragmentsInChr*> oChr;
  std::vector<JunctionCount*> oJC;
  std::vector<FragmentsMap*> oFM;
  std::vector<BAM2blocks*> BBchild;
  std::vector<BAMReader_Multi*> BRchild;

  for(unsigned int i = 0; i < n_threads_to_use; i++) {
    oCB.push_back(new CoverageBlocksIRFinder(CB_template));
    oSP.push_back(new SpansPoint(SP_template));
    oROI.push_back(new FragmentsInROI(ROI_template));
    oChr.push_back(new FragmentsInChr);
    oJC.push_back(new JunctionCount(JC_template));
    oFM.push_back(new FragmentsMap);
    BBchild.push_back(new BAM2blocks);
    BRchild.push_back(new BAMReader_Multi);

    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&JunctionCount::ChrMapUpdate, &(*oJC.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&JunctionCount::ProcessBlocks, &(*oJC.at(i)), std::placeholders::_1) );
    
    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&FragmentsInChr::ChrMapUpdate, &(*oChr.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&FragmentsInChr::ProcessBlocks, &(*oChr.at(i)), std::placeholders::_1) );
    
    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&SpansPoint::ChrMapUpdate, &(*oSP.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&SpansPoint::ProcessBlocks, &(*oSP.at(i)), std::placeholders::_1) );
        
    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&FragmentsInROI::ChrMapUpdate, &(*oROI.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&FragmentsInROI::ProcessBlocks, &(*oROI.at(i)), std::placeholders::_1) );
    
    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&CoverageBlocks::ChrMapUpdate, &(*oCB.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&CoverageBlocks::ProcessBlocks, &(*oCB.at(i)), std::placeholders::_1) );

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
    unsigned int n_blocks_read = 0;
    while(!BRchild.at(i)->eob() && !p.check_abort() && ret == 0) {
      n_blocks_read = (unsigned int)BRchild.at(i)->read_from_file(100);
      BRchild.at(i)->decompress();
      ret = BBchild.at(i)->processAll();
      
      blocks_read_total += n_blocks_read;
      p.increment(n_blocks_read);
    }
    
  }
#endif
  if(p.check_abort()) {
    // interrupted:
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      delete oJC.at(i);
      delete oChr.at(i);
      delete oSP.at(i);
      delete oROI.at(i);
      delete oCB.at(i);
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
      oJC.at(0)->Combine(*oJC.at(i));
      oChr.at(0)->Combine(*oChr.at(i));
      oSP.at(0)->Combine(*oSP.at(i));
      oROI.at(0)->Combine(*oROI.at(i));
      oCB.at(0)->Combine(*oCB.at(i));
      oFM.at(0)->Combine(*oFM.at(i));
      
      delete oJC.at(i);
      delete oChr.at(i);
      delete oSP.at(i);
      delete oROI.at(i);
      delete oCB.at(i);
      delete oFM.at(i);
    }
  }

  // Write Coverage Binary file:
  std::ofstream ofCOV;                          ofCOV.open(s_output_cov, std::ofstream::binary);  
  covWriter outCOV;                             outCOV.SetOutputHandle(&ofCOV);
  oFM.at(0)->WriteBinary(&outCOV, verbose);     ofCOV.close();

// Write output to file:  
	if(verbose) Rcout << "Writing output file\n";

  std::ofstream out;                            out.open(s_output_txt, std::ios::binary);  // Open binary file
  GZWriter outGZ;                               outGZ.SetOutputHandle(&out); // GZ compression

// Write stats here:
  BBchild.at(0)->WriteOutput(myLine);
  outGZ.writeline("BAM_report\tValue"); outGZ.writestring(myLine); outGZ.writeline("");

  int directionality = oJC.at(0)->Directional(myLine);
  outGZ.writeline("Directionality\tValue"); outGZ.writestring(myLine); outGZ.writeline("");

  // Generate output but save this to strings:
  std::string myLine_ROI;
  std::string myLine_JC;
  std::string myLine_SP;
  std::string myLine_Chr;
  std::string myLine_ND;
  std::string myLine_Dir;
  std::string myLine_QC;
  
  // Rcout << "Writing ROIs\n";
  oROI.at(0)->WriteOutput(myLine_ROI, myLine_QC);
  // Rcout << "Writing Juncs\n";
	oJC.at(0)->WriteOutput(myLine_JC, myLine_QC);
  // Rcout << "Writing Spans\n";
	oSP.at(0)->WriteOutput(myLine_SP, myLine_QC);
  // Rcout << "Writing Chrs\n";
	oChr.at(0)->WriteOutput(myLine_Chr, myLine_QC);
  // Rcout << "Writing CoverageBlocks\n";
	oCB.at(0)->WriteOutput(myLine_ND, myLine_QC, *oJC.at(0), *oSP.at(0), *oFM.at(0), n_threads_to_use);
  if (directionality != 0) {
    // Rcout << "Writing Stranded CoverageBlocks\n";
    oCB.at(0)->WriteOutput(myLine_Dir, myLine_QC, *oJC.at(0), *oSP.at(0), *oFM.at(0), n_threads_to_use, directionality); // Directional.
	}

  outGZ.writeline("QC\tValue"); outGZ.writestring(myLine_QC); outGZ.writeline("");
	
  outGZ.writeline("ROIname\ttotal_hits\tpositive_strand_hits\tnegative_strand_hits");
  outGZ.writestring(myLine_ROI); outGZ.writeline("");
  
  outGZ.writeline("JC_seqname\tstart\tend\tstrand\ttotal\tpos\tneg");
  outGZ.writestring(myLine_JC); outGZ.writeline("");
  
  outGZ.writeline("SP_seqname\tcoord\ttotal\tpos\tneg");
  outGZ.writestring(myLine_SP); outGZ.writeline("");
  
  outGZ.writeline("ChrCoverage_seqname\ttotal\tpos\tneg");
  outGZ.writestring(myLine_Chr); outGZ.writeline("");
  
  outGZ.writestring(myLine_ND); outGZ.writeline("");
  
  if (directionality != 0) {
    outGZ.writestring(myLine_Dir); outGZ.writeline("");
  }
  outGZ.flush(true);
  out.flush(); out.close();
  
  // destroy objects:

  delete oJC.at(0);
  delete oChr.at(0);
  delete oSP.at(0);
  delete oROI.at(0);
  delete oCB.at(0);
  delete oFM.at(0);
  delete BRchild.at(0);
  delete BBchild.at(0);

  return(0);
}



#ifndef GALAXY
// [[Rcpp::export]]
int IRF_main(std::string bam_file, std::string reference_file, std::string output_file, bool verbose = true, int n_threads = 1){
  
  std::string s_output_txt = output_file + ".txt.gz";
  std::string s_output_cov = output_file + ".cov";
#else
int IRF_main(std::string bam_file, std::string reference_file, std::string s_output_txt, std::string s_output_cov, int n_threads = 1){
	
  bool verbose = true;
#endif
  
  int use_threads = Set_Threads(n_threads);
  
  std::string s_bam = bam_file;
  std::string s_ref = reference_file;
		
  if(verbose) {
    Rcout << "Running IRFinder on " << s_bam;
    if(Has_OpenMP() != 0) Rcout << " with OpenMP ";
    Rcout << "using " << use_threads << " threads"
      << "\n" << "Reference: " << s_ref << "\n"
      << "Output file: " << s_output_txt << "\t" << s_output_cov << "\n\n"
      << "Reading reference file\n";
  }

  CoverageBlocksIRFinder * CB_template = new CoverageBlocksIRFinder;
  SpansPoint * SP_template = new SpansPoint;
  FragmentsInROI * ROI_template = new FragmentsInROI;
  JunctionCount * JC_template = new JunctionCount;
  
  int ret = 0;
  
  ret = IRF_ref(s_ref, *CB_template, *SP_template, *ROI_template, *JC_template, verbose);
  if(ret != 0) {
    Rcout << "Reading Reference file failed. Check if IRFinder.ref.gz exists and is a valid NxtIRF-generated IRFinder reference\n";
    return(ret);
  }
  // main:
  ret = IRF_core(s_bam, s_output_txt, s_output_cov,
    *CB_template, *SP_template, *ROI_template, *JC_template, verbose, use_threads);
    
  if(ret != 0) Rcout << "Process interrupted running IRFinder on " << s_bam << '\n';
  
  delete CB_template;
  delete SP_template;
  delete ROI_template;
  delete JC_template;
  return(ret);
}

#ifndef GALAXY
// [[Rcpp::export]]
int IRF_main_multi(std::string reference_file, StringVector bam_files, StringVector output_files, int max_threads, bool verbose = true){
	
	int use_threads = Set_Threads(max_threads);

	if(bam_files.size() != output_files.size() || bam_files.size() < 1) {
		Rcout << "bam_files and output_files are of different sizes\n";
		return(1);	
	}
	
	std::vector< std::string > v_bam;
	std::vector< std::string > v_out;
  for(int z = 0; z < bam_files.size(); z++) {
		v_bam.push_back(string(bam_files(z)));
		v_out.push_back(string(output_files(z)));
	}

  std::string s_ref = reference_file;
  Rcout << "Reading reference file\n";
  
  CoverageBlocksIRFinder * CB_template = new CoverageBlocksIRFinder;
  SpansPoint * SP_template = new SpansPoint;
  FragmentsInROI * ROI_template = new FragmentsInROI;
  JunctionCount * JC_template = new JunctionCount;
  
  int ret = 0;
  
  ret = IRF_ref(s_ref, *CB_template, *SP_template, *ROI_template, *JC_template, false);
  if(ret != 0) {
    Rcout << "Reading Reference file failed. Check if IRFinder.ref.gz exists and is a valid NxtIRF-generated IRFinder reference\n";
    return(ret);
  }

	Rcout << "Running IRFinder with OpenMP using " << use_threads << " threads\n";

  for(unsigned int z = 0; z < v_bam.size(); z++) {
    std::string s_bam = v_bam.at(z);
		std::string s_output_txt = v_out.at(z) + ".txt.gz";
		std::string s_output_cov = v_out.at(z) + ".cov";
		// Rcout << "Processing " << s_bam << " with output " << v_out.at(z) << '\n';
    
    int ret2 = IRF_core(s_bam, s_output_txt, s_output_cov,
      *CB_template, *SP_template, *ROI_template, *JC_template, verbose, use_threads);
    if(ret2 != 0) {
      Rcout << "Process interrupted running IRFinder on " << s_bam << '\n';
      delete CB_template;
      delete SP_template;
      delete ROI_template;
      delete JC_template;
      return(ret);
    } else {
      Rcout << s_bam << " processed\n";
    }
	}

  delete CB_template;
  delete SP_template;
  delete ROI_template;
  delete JC_template;
  return(0);
}

#else
// Galaxy main
int main(int argc, char * argv[]) {
	// Usage:
    // irfinder_galaxy main sample.bam IRFinder.ref.gz OutputHeader
    // irfinder_galaxy gen_map_reads genome.fa reads_to_map.fa 70 10
    // irfinder_galaxy process_mappability_bam mappedreads.bam mappability.bed

  if(std::string(argv[1]) == "gen_map_reads") {
      std::string s_genome = argv[2];
      std::string s_output = argv[3];
      int read_len = atoi(argv[4]);
      int read_stride = atoi(argv[5]);
      int read_error = atoi(argv[4]) / 2;
      IRF_GenerateMappabilityReads(s_genome, s_output, read_len, read_stride, read_error);
      exit(0);
  } else if(std::string(argv[1]) == "process_mappability_bam") {
      std::string s_bam = argv[2];
      std::string s_output = argv[3];
      int threshold = atoi(argv[4]);
      if(argc == 6) {
        std::string s_cov = argv[5];
        IRF_GenerateMappabilityRegions(s_bam, s_output, threshold, s_cov);
        exit(0);          
      } else {
        IRF_GenerateMappabilityRegions(s_bam, s_output, threshold);
        exit(0);
      }
  } else if(std::string(argv[1]) == "main") {
      std::string s_bam = argv[2];
      std::string s_ref = argv[3];
      std::string s_output_txt = argv[4];		
      std::string s_output_cov = argv[5];		
  IRF_main(s_bam, s_ref, s_output_txt, s_output_cov);
  } else {
    Rcout << "Usage:\n\t"
      << argv[0] <<  " main samplename.bam IRFinder.ref.gz samplename.txt.gz samplename.cov\n"
      << argv[0] <<  " main samplename.bam IRFinder.ref.gz samplename.txt.gz samplename.cov\n"
      << argv[0] <<  " irfinder_galaxy gen_map_reads genome.fa reads_to_map.fa 70 10\n"
      << argv[0] <<  " irfinder_galaxy process_mappability_bam mappedreads.bam mappability.bed {mappability.cov}";   
  }
}
	
#endif