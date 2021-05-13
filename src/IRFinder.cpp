#include "ReadBlockProcessor.h"
#include "ReadBlockProcessor_CoverageBlocks.h"
#include "BAM2blocks.h"
#include "GZReader.h"
#include "Mappability.h"

#include "includedefine.h"

#ifndef GALAXY

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
int Has_OpenMP() {
#ifdef _OPENMP
	return omp_get_max_threads();
#else
	return 0;
#endif
}

const char refEOF[5] =
		"\x20\x45\x4f\x46";

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
  int ref_index;
  auto it_chr = std::find(inCov.chr_names.begin(), inCov.chr_names.end(), seqname);
  if(it_chr == inCov.chr_names.end()) {
		inCov_stream.close();	
    return(NULL_RLE);
  } else {
    ref_index = distance(inCov.chr_names.begin(), it_chr);
  }
  // end = 0 implies fetch whole chromosome
  int eff_end = 0;
  if(end == 0) {
    eff_end = (int)inCov.chr_lens.at(ref_index);
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
  if((uint32_t)eff_end < inCov.chr_lens.at(ref_index)) {
    values.push_back(0);
    lengths.push_back(inCov.chr_lens.at(ref_index) - eff_end);
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
  
  for (unsigned int i = 0; i < inCov.chr_names.size(); i++) {
    uint32_t eff_end = inCov.chr_lens.at(i);
    uint32_t start = 0;
    
    std::vector<int> values;
    std::vector<unsigned int> lengths;

    inCov.FetchRLE(inCov.chr_names.at(i), (uint32_t)start, (uint32_t)eff_end, strand, &values, &lengths);
    
    List RLE = List::create(
      _["values"] = values,
      _["lengths"] = lengths 
    );
    RLEList.push_back(RLE, inCov.chr_names.at(i));
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

#ifndef GALAXY
// [[Rcpp::export]]
int IRF_main(std::string bam_file, std::string reference_file, std::string output_file, bool verbose = true){
  
  std::string s_output_txt = output_file + ".txt.gz";
  std::string s_output_cov = output_file + ".cov";
#else
int IRF_main(std::string bam_file, std::string reference_file, std::string s_output_txt, std::string s_output_cov){	
	bool verbose = true;
#endif

  std::string s_bam = bam_file;
  
  std::string s_ref = reference_file;
		
		if(verbose) {
			Rcout << "Running IRFinder on " << s_bam << "\nReference: " << reference_file << "\n";
			Rcout << "Output file: " << s_output_txt << "\t" << s_output_cov << "\n\n";

			Rcout << "Reading reference file\n";
		}

    
    GZReader gz_in;
    int ret = gz_in.LoadGZ(reference_file, true);
		if(ret != 0) return(-1);
		
    std::string myLine;
    std::string myBuffer;
    
    getline(gz_in.iss, myLine, '#');    // discard first >
    getline(gz_in.iss, myLine, '\n');   // ignore file names for now
		
    getline(gz_in.iss, myBuffer, '#');  // this is the data block for ref-cover.bed

  CoverageBlocksIRFinder oCoverageBlocks;
  std::istringstream inCoverageBlocks;
  inCoverageBlocks.str(myBuffer);
  oCoverageBlocks.loadRef(inCoverageBlocks);

    getline(gz_in.iss, myLine, '\n');
    getline(gz_in.iss, myBuffer, '#');

  SpansPoint oSpansPoint;
  oSpansPoint.setSpanLength(5,4);
  std::istringstream inSpansPoint;
  inSpansPoint.str(myBuffer);
  oSpansPoint.loadRef(inSpansPoint);

    getline(gz_in.iss, myLine, '\n');
    getline(gz_in.iss, myBuffer, '#');
  
  FragmentsInROI oFragmentsInROI;
  FragmentsInChr oFragmentsInChr;

    std::istringstream inFragmentsInROI;
    inFragmentsInROI.str(myBuffer);
    oFragmentsInROI.loadRef(inFragmentsInROI);

    getline(gz_in.iss, myLine, '\n');
    getline(gz_in.iss, myBuffer, '#');

  JunctionCount oJuncCount;
  std::istringstream inJuncCount;
  inJuncCount.str(myBuffer);
  oJuncCount.loadRef(inJuncCount);

	// Ensure valid reference termination:
		getline(gz_in.iss, myLine, '\n');    
		if(strncmp(myLine.c_str(), refEOF, 4) != 0) {
			Rcout << "Invalid IRFinder reference detected\n";
			return(0);
		}
  
  FragmentsMap oFragMap;
  
  BAM2blocks BB;
  
  BB.registerCallbackChrMappingChange( std::bind(&JunctionCount::ChrMapUpdate, &oJuncCount, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&JunctionCount::ProcessBlocks, &oJuncCount, std::placeholders::_1) );
  
  BB.registerCallbackChrMappingChange( std::bind(&FragmentsInChr::ChrMapUpdate, &oFragmentsInChr, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&FragmentsInChr::ProcessBlocks, &oFragmentsInChr, std::placeholders::_1) );
  
  BB.registerCallbackChrMappingChange( std::bind(&SpansPoint::ChrMapUpdate, &oSpansPoint, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&SpansPoint::ProcessBlocks, &oSpansPoint, std::placeholders::_1) );
      
  BB.registerCallbackChrMappingChange( std::bind(&FragmentsInROI::ChrMapUpdate, &oFragmentsInROI, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&FragmentsInROI::ProcessBlocks, &oFragmentsInROI, std::placeholders::_1) );
  
  BB.registerCallbackChrMappingChange( std::bind(&CoverageBlocks::ChrMapUpdate, &oCoverageBlocks, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&CoverageBlocks::ProcessBlocks, &oCoverageBlocks, std::placeholders::_1) );

  BB.registerCallbackChrMappingChange( std::bind(&FragmentsMap::ChrMapUpdate, &oFragMap, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&FragmentsMap::ProcessBlocks, &oFragMap, std::placeholders::_1) );

	if(verbose) {  
		Rcout << "Processing BAM file\n";
  }  
  BAMReader inbam;
  std::ifstream inbam_stream;
  inbam_stream.open(s_bam, std::ios::in | std::ios::binary);
  inbam.SetInputHandle(&inbam_stream);
  
  BB.openFile(&inbam); // This file needs to be a decompressed BAM. (setup via fifo / or expect already decompressed via stdin).
  BB.processAll(myLine, verbose);
	// oFragMap.sort_and_collapse_final(verbose);

  // Write Coverage Binary file:
  
  std::ofstream ofCOV;
  ofCOV.open(s_output_cov, std::ofstream::binary);
   
  covFile outCOV;
  outCOV.SetOutputHandle(&ofCOV);
  
  oFragMap.WriteBinary(&outCOV, BB.chr_names, BB.chr_lens, verbose);
  ofCOV.close();

// Write output to file:  
	if(verbose) {  
		Rcout << "Writing output file\n";
	}
  std::ofstream out;
  out.open(s_output_txt, std::ios::binary);

// GZ compression:
  GZWriter outGZ;
  outGZ.SetOutputHandle(&out);

// Write stats here:

  outGZ.writeline("BAM_report\tValue");
  outGZ.writestring(myLine);
  outGZ.writeline("");

  int directionality = oJuncCount.Directional(myLine);
  outGZ.writeline("Directionality\tValue");
  outGZ.writestring(myLine);
  outGZ.writeline("");

// Generate output but save this to strings:
std::string myLine_ROI;
std::string myLine_JC;
std::string myLine_SP;
std::string myLine_Chr;
std::string myLine_ND;
std::string myLine_Dir;
std::string myLine_QC;

  oFragmentsInROI.WriteOutput(myLine_ROI, myLine_QC);
	oJuncCount.WriteOutput(myLine_JC, myLine_QC);
	oSpansPoint.WriteOutput(myLine_SP, myLine_QC);
	oFragmentsInChr.WriteOutput(myLine_Chr, myLine_QC);
	oCoverageBlocks.WriteOutput(myLine_ND, myLine_QC, oJuncCount, oSpansPoint, oFragMap);
  if (directionality != 0) {
    oCoverageBlocks.WriteOutput(myLine_Dir, myLine_QC, oJuncCount, oSpansPoint, oFragMap, directionality); // Directional.
	}

  outGZ.writeline("QC\tValue");
  outGZ.writestring(myLine_QC);
  outGZ.writeline("");
	
  outGZ.writeline("ROIname\ttotal_hits\tpositive_strand_hits\tnegative_strand_hits");
  outGZ.writestring(myLine_ROI);
  outGZ.writeline("");
  
  outGZ.writeline("JC_seqname\tstart\tend\tstrand\ttotal\tpos\tneg");
  outGZ.writestring(myLine_JC);
  outGZ.writeline("");
  
  outGZ.writeline("SP_seqname\tcoord\ttotal\tpos\tneg");
  outGZ.writestring(myLine_SP);
  outGZ.writeline("");
  
  outGZ.writeline("ChrCoverage_seqname\ttotal\tpos\tneg");
  outGZ.writestring(myLine_Chr);
  outGZ.writeline("");
  
  outGZ.writestring(myLine_ND);
  outGZ.writeline("");
  
  if (directionality != 0) {
    outGZ.writestring(myLine_Dir);
    outGZ.writeline("");
  }
  outGZ.flush(true);
  out.flush(); out.close();
  

  
  return(0);
}

#ifndef GALAXY

#ifndef _OPENMP
int IRF_main_multithreaded(std::string reference_file, StringVector bam_files, StringVector output_files, int max_threads){
	Rcout << "NxtIRF was built without OpenMP; exiting...";
	return(1);
}
#else
// [[Rcpp::export]]
int IRF_main_multithreaded(std::string reference_file, StringVector bam_files, StringVector output_files, int max_threads){
	
	int use_threads = 0;
	if(max_threads > 0 && max_threads <= omp_get_thread_limit()) {
		if(max_threads > bam_files.size()) {
			use_threads = bam_files.size();
		} else {
			use_threads = max_threads;
		}
	} else {
		use_threads = omp_get_thread_limit();
		if(use_threads < 1) {
			use_threads = 1;
		}
	}
	omp_set_num_threads(use_threads);

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

	GZReader gz_in;
	int ret = gz_in.LoadGZ(reference_file, true);
	if(ret != 0) return(-1);

	Rcout << "Running IRFinder with OpenMP using " << use_threads << " threads\n";
	#pragma omp parallel for
  for(unsigned int z = 0; z < v_bam.size(); z++) {
    std::string s_bam = v_bam.at(z);
    std::string output_file = v_out.at(z);

		std::string s_output_txt = output_file + ".txt.gz";
		std::string s_output_cov = output_file + ".cov";
		
    std::string myLine;
    std::string myBuffer;
    
		std::istringstream iss(gz_in.iss.str());
		
    getline(iss, myLine, '>');    // discard first >
    getline(iss, myLine, '\n');   // ignore file names for now
    getline(iss, myBuffer, '>');  // this is the data block for ref-cover.bed

		CoverageBlocksIRFinder oCoverageBlocks;
		std::istringstream inCoverageBlocks;
		inCoverageBlocks.str(myBuffer);
		oCoverageBlocks.loadRef(inCoverageBlocks);

		getline(iss, myLine, '\n');
		getline(iss, myBuffer, '>');

		SpansPoint oSpansPoint;
		oSpansPoint.setSpanLength(5,4);
		std::istringstream inSpansPoint;
		inSpansPoint.str(myBuffer);
		oSpansPoint.loadRef(inSpansPoint);

		getline(iss, myLine, '\n');
		getline(iss, myBuffer, '>');
		
		FragmentsInROI oFragmentsInROI;
		FragmentsInChr oFragmentsInChr;

		std::istringstream inFragmentsInROI;
		inFragmentsInROI.str(myBuffer);
		oFragmentsInROI.loadRef(inFragmentsInROI);

		getline(iss, myLine, '\n');
		getline(iss, myBuffer, '>');

		JunctionCount oJuncCount;
		std::istringstream inJuncCount;
		inJuncCount.str(myBuffer);
		oJuncCount.loadRef(inJuncCount);
		
		FragmentsMap oFragMap;
		
		BAM2blocks BB;
  
		BB.registerCallbackChrMappingChange( std::bind(&JunctionCount::ChrMapUpdate, &oJuncCount, std::placeholders::_1) );
		BB.registerCallbackProcessBlocks( std::bind(&JunctionCount::ProcessBlocks, &oJuncCount, std::placeholders::_1) );
		
		BB.registerCallbackChrMappingChange( std::bind(&FragmentsInChr::ChrMapUpdate, &oFragmentsInChr, std::placeholders::_1) );
		BB.registerCallbackProcessBlocks( std::bind(&FragmentsInChr::ProcessBlocks, &oFragmentsInChr, std::placeholders::_1) );
		
		BB.registerCallbackChrMappingChange( std::bind(&SpansPoint::ChrMapUpdate, &oSpansPoint, std::placeholders::_1) );
		BB.registerCallbackProcessBlocks( std::bind(&SpansPoint::ProcessBlocks, &oSpansPoint, std::placeholders::_1) );
		
		BB.registerCallbackChrMappingChange( std::bind(&FragmentsInROI::ChrMapUpdate, &oFragmentsInROI, std::placeholders::_1) );
		BB.registerCallbackProcessBlocks( std::bind(&FragmentsInROI::ProcessBlocks, &oFragmentsInROI, std::placeholders::_1) );
		
		BB.registerCallbackChrMappingChange( std::bind(&CoverageBlocks::ChrMapUpdate, &oCoverageBlocks, std::placeholders::_1) );
		BB.registerCallbackProcessBlocks( std::bind(&CoverageBlocks::ProcessBlocks, &oCoverageBlocks, std::placeholders::_1) );

		BB.registerCallbackChrMappingChange( std::bind(&FragmentsMap::ChrMapUpdate, &oFragMap, std::placeholders::_1) );
		BB.registerCallbackProcessBlocks( std::bind(&FragmentsMap::ProcessBlocks, &oFragMap, std::placeholders::_1) );
		
		BAMReader inbam;
		std::ifstream inbam_stream;
		inbam_stream.open(s_bam, std::ios::in | std::ios::binary);
		inbam.SetInputHandle(&inbam_stream);
		
		Rcout << "Processing " << s_bam << "\n";
		
		BB.openFile(&inbam); // This file needs to be a decompressed BAM. (setup via fifo / or expect already decompressed via stdin).
		BB.processAll(myLine, false);

		// Write Coverage Binary file:
		
		std::ofstream ofCOV;
		ofCOV.open(s_output_cov, std::ofstream::binary);
		 
		covFile outCOV;
		outCOV.SetOutputHandle(&ofCOV);
		
		oFragMap.WriteBinary(&outCOV, BB.chr_names, BB.chr_lens);
		ofCOV.close();	
		
		std::ofstream out;
		out.open(s_output_txt, std::ios::binary);

	// GZ compression:
		GZWriter outGZ;
		outGZ.SetOutputHandle(&out);

	// Write stats here:

		outGZ.writeline("BAM_report\tValue");
		outGZ.writestring(myLine);
		outGZ.writeline("");

		int directionality = oJuncCount.Directional(myLine);
		outGZ.writeline("Directionality\tValue");
		outGZ.writestring(myLine);
		outGZ.writeline("");

		// Generate output but save this to strings:
		std::string myLine_ROI;
		std::string myLine_JC;
		std::string myLine_SP;
		std::string myLine_Chr;
		std::string myLine_ND;
		std::string myLine_Dir;
		std::string myLine_QC;

		oFragmentsInROI.WriteOutput(myLine_ROI, myLine_QC);
		oJuncCount.WriteOutput(myLine_JC, myLine_QC);
		oSpansPoint.WriteOutput(myLine_SP, myLine_QC);
		oFragmentsInChr.WriteOutput(myLine_Chr, myLine_QC);
		oCoverageBlocks.WriteOutput(myLine_ND, myLine_QC, oJuncCount, oSpansPoint, oFragMap);
		if (directionality != 0) {
			oCoverageBlocks.WriteOutput(myLine_Dir, myLine_QC, oJuncCount, oSpansPoint, oFragMap, directionality); // Directional.
		}

		outGZ.writeline("QC\tValue");
		outGZ.writestring(myLine_QC);
		outGZ.writeline("");
		
		outGZ.writeline("ROIname\ttotal_hits\tpositive_strand_hits\tnegative_strand_hits");
		outGZ.writestring(myLine_ROI);
		outGZ.writeline("");
		
		outGZ.writeline("JC_seqname\tstart\tend\tstrand\ttotal\tpos\tneg");
		outGZ.writestring(myLine_JC);
		outGZ.writeline("");
		
		outGZ.writeline("SP_seqname\tcoord\ttotal\tpos\tneg");
		outGZ.writestring(myLine_SP);
		outGZ.writeline("");
		
		outGZ.writeline("ChrCoverage_seqname\ttotal\tpos\tneg");
		outGZ.writestring(myLine_Chr);
		outGZ.writeline("");
		
		outGZ.writestring(myLine_ND);
		outGZ.writeline("");
		
		if (directionality != 0) {
			outGZ.writestring(myLine_Dir);
			outGZ.writeline("");
		}
		outGZ.flush(true);
		out.flush(); out.close();
		
	
	}

	return(0);
}
#endif

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