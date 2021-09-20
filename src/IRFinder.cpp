/* IRFinder.cpp Main IRFinder exported functions

Copyright (C) 2021 Alex Chit Hei Wong
Copyright (C) 2016 William Ritchie
  - original: https://github.com/williamritchie/IRFinder/tree/IRFinder-1.3.1)

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

int ReadChrAlias(std::istringstream &IN,
    std::vector<std::string> &ref_names, 
    std::vector<std::string> &ref_alias,
    std::vector<uint32_t> &ref_lengths
) {
  ref_names.clear();
  ref_alias.clear();
  
  std::string myLine;
  myLine.reserve(1000);
  std::string myChr;
  myChr.reserve(100);
  std::string myAlias;
  myAlias.reserve(100);
  std::string myLength;
  myLength.reserve(100);
  
  while(!IN.eof() && !IN.fail()) {
    getline(IN, myLine, '\n');
    if (IN.eof() || IN.fail()) {
      if (myLine.length() == 0) {
        // This line is empty - just a blank line at the end of the file.
        // Checking at this stage allows correct handling of files both with and without a trailing \n after the last record.
        break;
      }else{
        // Error line in input, ignore.
        break;
      }
    }
    std::istringstream lineStream;
    lineStream.str(myLine);
    getline(lineStream, myChr, '\t');
    getline(lineStream, myLength, '\t');
    getline(lineStream, myAlias, '\t');
    if(myChr.size() > 0) {
      ref_names.push_back(myChr);
      ref_lengths.push_back((uint32_t)stoul(myLength));
      ref_alias.push_back(myAlias);      
    }
  }
  // Rcout << "Debug:" << ref_names.size() << " chromosome aliases loaded\n";
  return(0);
}

// IRFinder reference reader:
int IRF_ref(std::string &reference_file, 
    std::vector<std::string> &ref_names, 
    std::vector<std::string> &ref_alias,
    std::vector<uint32_t> &ref_lengths,
    CoverageBlocksIRFinder &CB_template, 
    SpansPoint &SP_template, 
    FragmentsInROI &ROI_template,
    JunctionCount &JC_template, 
    bool verbose
) { 
  GZReader * gz_in = new GZReader;
  int ret = gz_in->LoadGZ(reference_file, true);
  if(ret != 0) return(-1);
  
  // Allows IRFinder reference blocks to be read in any order
  std::string headerCover ("ref-cover.bed");
  std::string headerSpans ("ref-read-continues.ref");
  std::string headerROI ("ref-ROI.bed");
  std::string headerSJ ("ref-sj.ref");
  std::string headerChr ("ref-chrs.ref");
  std::string headerEOF ("EOF");
  
  bool doneCover = false;
  bool doneSpans = false;
  bool doneROI = false;
  bool doneSJ = false;
  bool doneChrs = false;
  
  std::string myLine;
  std::string myBuffer;
  
  getline(gz_in->iss, myLine, '#');    // discard anything before the first hash
  getline(gz_in->iss, myLine, '\n');   // Get block name
  
  // Check non-empty IRF ref block name
  if(myLine.size() == 0) {
    Rcout << "Invalid IRFinder reference detected\n";
    return(-1);
  }

  while(myLine.find(headerEOF)==std::string::npos) {
    getline(gz_in->iss, myBuffer, '#');  // this is the data block
    
    // Check that data block is not empty.
    // if(myBuffer.size() == 0) {
      // Rcout << "Invalid IRFinder reference detected\n";
      // return(-1);
    // }
    
    if(myLine.find(headerCover)!=std::string::npos && !doneCover) {
      std::istringstream inCoverageBlocks;
      inCoverageBlocks.str(myBuffer);
      CB_template.loadRef(inCoverageBlocks);
      // Rcout << "doneCover\n";
      doneCover = true;
    } else if(myLine.find(headerSpans)!=std::string::npos && !doneSpans) {
      SP_template.setSpanLength(5,4);
      std::istringstream inSpansPoint;
      inSpansPoint.str(myBuffer);
      SP_template.loadRef(inSpansPoint);
      // Rcout << "doneSpans\n";
      doneSpans = true;
    } else if(myLine.find(headerROI)!=std::string::npos && !doneROI) {
      std::istringstream inFragmentsInROI;
      inFragmentsInROI.str(myBuffer);
      ROI_template.loadRef(inFragmentsInROI);
      // Rcout << "doneROI\n";
      doneROI = true;
    } else if(myLine.find(headerSJ)!=std::string::npos && !doneSJ) {
      std::istringstream inJuncCount;
      inJuncCount.str(myBuffer);
      JC_template.loadRef(inJuncCount);
      // Rcout << "doneSJ\n";
      doneSJ = true;
    } else if(myLine.find(headerChr)!=std::string::npos && !doneChrs) {
      std::istringstream inChrAlias;
      inChrAlias.str(myBuffer);
      ReadChrAlias(inChrAlias, ref_names, ref_alias, ref_lengths);
      // Rcout << "doneChrs\n";
      doneChrs = true;
    } else {
      Rcout << "Error: Invalid IRFinder reference block detected\n";
      return(-1);
    }

    // Get next data block name
    getline(gz_in->iss, myLine, '\n');
  }
  
  if(!doneCover || !doneSpans || !doneROI || !doneSJ) {
    Rcout << "Error: Incomplete IRFinder reference detected\n";
    return(-1);
  }
  return(0);
}

// IRFinder core:
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
) {
  unsigned int n_threads_to_use = (unsigned int)n_threads;   // Should be sorted out in calling function
 
  std::string myLine;
	if(verbose) Rcout << "Processing BAM file " << bam_file << "\n";
  
  
  // std::ifstream inbam_stream;   inbam_stream.open(bam_file, std::ios::in | std::ios::binary);
  // BAMReader_Multi inbam;        inbam.SetInputHandle(&inbam_stream); // Rcout << "BAMReader_Multi handle set\n";  
  pbam_in inbam((size_t)5e8, (size_t)1e9, 5);
  // inbam.SetInputHandle(&inbam_stream, n_threads_to_use);
  inbam.openFile(bam_file, n_threads_to_use);
  
  // Abort here if BAM corrupt
  std::vector<std::string> s_chr_names;
  std::vector<uint32_t> u32_chr_lens;
  int chrcount = inbam.obtainChrs(s_chr_names, u32_chr_lens);
  if(chrcount < 1) {
    Rcout << bam_file << " - contains no chromosomes mapped\n";
    return(-1);
  }
  
  // Compile here a list of chromosomes; use BAM chromosomes for order
  // Add reference-only chromosomes at the end
  std::vector<std::string> bam_chr_name;
  std::vector<uint32_t> bam_chr_len;
  for(unsigned int i = 0; i < s_chr_names.size(); i++) {
    for(unsigned int j = 0; j < ref_alias.size(); j++) {
      if( 0==strncmp(
            ref_alias.at(j).c_str(), 
            s_chr_names.at(i).c_str(), 
            s_chr_names.at(i).size()
          ) && s_chr_names.at(i).size() == ref_alias.at(j).size()
      ) {
        bam_chr_name.push_back(ref_names.at(j));
        bam_chr_len.push_back(u32_chr_lens.at(i));
        break;
      }
    }
    if(i == bam_chr_name.size()) {
      bam_chr_name.push_back(s_chr_names.at(i));
      bam_chr_len.push_back(u32_chr_lens.at(i));
    }
  }
  // Now fill in reference chromosomes not in BAM:
  for(unsigned int i = 0; i < ref_names.size(); i++) {
    auto it = std::find(bam_chr_name.begin(), bam_chr_name.end(), ref_names.at(i));
    if(it == bam_chr_name.end()) {
      bam_chr_name.push_back(ref_names.at(i));
      bam_chr_len.push_back(ref_lengths.at(i));      
    }
  }
  
  std::vector<CoverageBlocksIRFinder*> oCB;
  std::vector<SpansPoint*> oSP;
  std::vector<FragmentsInROI*> oROI;
  std::vector<FragmentsInChr*> oChr;
  std::vector<JunctionCount*> oJC;
  std::vector<FragmentsMap*> oFM;
  std::vector<BAM2blocks*> BBchild;

  for(unsigned int i = 0; i < n_threads_to_use; i++) {
    oCB.push_back(new CoverageBlocksIRFinder(CB_template));
    oSP.push_back(new SpansPoint(SP_template));
    oROI.push_back(new FragmentsInROI(ROI_template));
    oChr.push_back(new FragmentsInChr);
    oJC.push_back(new JunctionCount(JC_template));
    oFM.push_back(new FragmentsMap);
    BBchild.push_back(new BAM2blocks(ref_names, ref_lengths));

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
      delete oJC.at(i);
      delete oChr.at(i);
      delete oSP.at(i);
      delete oROI.at(i);
      delete oCB.at(i);
      delete oFM.at(i);
      delete BBchild.at(i);
    }
    return(-1);
  }

  // inbam.closeFile();
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

  if(verbose) Rcout << "Writing COV file\n";

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
  
  std::vector<std::string> ref_names;
  std::vector<std::string> ref_alias;
  std::vector<uint32_t> ref_lengths;
  
  int ret = 0;
  
  ret = IRF_ref(s_ref, ref_names, ref_alias, ref_lengths,
    *CB_template, *SP_template, *ROI_template, *JC_template, verbose
  );
  if(ret != 0) {
    Rcout << "Reading Reference file failed. Check if IRFinder.ref.gz exists and is a valid NxtIRF-generated IRFinder reference\n";
    return(ret);
  }
  // main:
  ret = IRF_core(s_bam, s_output_txt, s_output_cov,
    ref_names, ref_alias, ref_lengths,
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
  
  std::vector<std::string> ref_names;
  std::vector<std::string> ref_alias;
  std::vector<uint32_t> ref_lengths;
  
  int ret = 0;
  
  ret = IRF_ref(s_ref, ref_names, ref_alias, ref_lengths,
    *CB_template, *SP_template, *ROI_template, *JC_template, false);
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
      ref_names, ref_alias, ref_lengths,
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