// WARNING: code is little endian only!

#include "BAM2blocks.h"
#include "includedefine.h"
// using namespace std;

BAM2blocks::BAM2blocks() {
	oBlocks = FragmentBlocks(); //Right syntax to call the default constructor on an object variable, declared but not initialised?

	cShortPairs = 0;
	cIntersectPairs = 0;
	cLongPairs = 0;
	cSingleReads = 0;
	cPairedReads = 0;
	cErrorReads = 0;
	cSkippedReads = 0;
	cChimericReads = 0;
}

// OK.
void BAM2blocks::readBamHeader() {
  char buffer[1000];
  std::string chrName;

  bam_header bamhead;

  IN->read(bamhead.c, BAM_HEADER_BYTES);

  char * headertext = new char[bamhead.magic.l_text+1];
  IN->read(headertext, bamhead.magic.l_text);
  samHeader = string(headertext, bamhead.magic.l_text);
  delete[] headertext;
  
  stream_int32 i32;
  IN->read(i32.c ,4);
  unsigned int n_chr = i32.i;


  for (unsigned int i = 0; i < n_chr; i++) {
    IN->read(i32.c ,4);
    IN->read(buffer , i32.i);
    chrName = string(buffer, i32.i-1);
    chr_names.push_back(chrName);

    IN->read(i32.c ,4);
    chr_lens.push_back(i32.i);
    chrs.push_back(chr_entry(i, chrName, i32.i));
  }
  std::sort(chrs.begin(), chrs.end());
  
	for (auto & callback : callbacksChrMappingChange ) {
		callback(chrs);
	}
  
}

void BAM2blocks::cigar2block(uint32_t * cigar, uint16_t n_cigar_op, std::vector<int> &starts, std::vector<int> &lens, int &ret_genome_len) {
  bool inBlock = true;
  int relpos = 0;
  int curblock = 0;
  starts.resize(1);  // Is this expensive or not -- does this call destroy on further items, or is it a single op, adjusting the end? If expensive we can revert to earlier behaviour where we keep track of how many blocks, just overwriting relevant parts of the vector.
  lens.resize(1);
  starts[curblock] = 0;
  lens[curblock] = 0;

  for (; n_cigar_op > 0; n_cigar_op--) {
    if (inBlock) {
      switch (*cigar & 15) {
        case 0: case 2: case 7: case 8:
          // increment len of last block
          lens[curblock] += (*cigar >> 4);
          relpos += (*cigar >> 4);
          break;
        case 3:
          curblock++;
          relpos += (*cigar >> 4);
          // extend arrays. 
          starts.push_back(relpos);
          lens.push_back(0);
          inBlock = false;
          break;
      }
    }else{
      switch (*cigar & 15) {
        case 0: case 2: case 7: case 8:
          lens[curblock] = (*cigar >> 4);
          relpos += (*cigar >> 4);
          inBlock = true;
          break;
        case 3:
          // push start of next further out
          relpos += (*cigar >> 4);
          starts[curblock] = relpos;
          break;
        }
    }
    cigar++;
  }
  ret_genome_len = relpos;
//  *ret_blocks = curblock+1;  // Unnecessary if we are using vectors in the expected manner - changing their length as needed.
}


//OK - translated - doesn't call the callbacks yet though.
unsigned int BAM2blocks::processPair(bam_read_core * read1, bam_read_core * read2) {
  // R1 is to the left of R2 (or equal starts).
  int r1_genome_len;
  //int r1_blocks;
  int r2_genome_len;

  string debugstate;

  //int r2_blocks;
  //char dir;

  if (read1->core.flag & 0x40) {
    //this is first of pair.
    if (read1->core.flag & 0x10) {
      oBlocks.direction = 0;
    }else{
      oBlocks.direction = 1;
    }
  }else{
    if (read1->core.flag & 0x20) {
      oBlocks.direction = 0;
    }else{
      oBlocks.direction = 1;
    }
  }


  cigar2block(read1->cigar, read1->core.n_cigar_op, oBlocks.rStarts[0], oBlocks.rLens[0], r1_genome_len);
  cigar2block(read2->cigar, read2->core.n_cigar_op, oBlocks.rStarts[1], oBlocks.rLens[1], r2_genome_len);

  if (read1->core.pos + r1_genome_len < read2->core.pos) {
    cLongPairs++;
    //reads do not intersect
    oBlocks.readCount = 2;
    debugstate.append( "-Long-");
  }else if (read1->core.pos + r1_genome_len >= read2->core.pos + r2_genome_len){
    cShortPairs++;
    // Read 2 is a short read & read 1 fully contains it (or perhaps just a trimmed read with two exactly complementary reads remaining).
    oBlocks.readCount = 1;
    debugstate.append( "-Short-");    
  }else{
    debugstate.append( "-Intersect-");
    cIntersectPairs++;
    bool goodPair = true;
    oBlocks.readCount = 1;
    // We have two reads that intersect - construct just one complete fragment.

// Guaranteed assumptions:
//   Read 1 starts to the left of Read 2.
//   Read 2 end extends beyond the end of Read 1 end.
    int r1pos = read1->core.pos;
    int r2pos = read2->core.pos;
    for (unsigned int i = 0; i < oBlocks.rStarts[0].size(); i++) {
        if (r1pos + oBlocks.rStarts[0][i] + oBlocks.rLens[0][i] >= r2pos) {
          if (r1pos + oBlocks.rStarts[0][i] <= r2pos) {
            oBlocks.rLens[0][i] = r2pos - r1pos - oBlocks.rStarts[0][i] + oBlocks.rLens[1][0];
            //r1_blocks = i + r2_blocks;
            oBlocks.rStarts[0].resize(i + oBlocks.rStarts[1].size());
            oBlocks.rLens[0].resize(i + oBlocks.rStarts[1].size());
            // Maybe this can be optimised by using push_back below instead of running resize.
            for (unsigned int j = 1; j < oBlocks.rStarts[1].size(); j++) {
              i++;
              oBlocks.rLens[0][i] = oBlocks.rLens[1][j];
              oBlocks.rStarts[0][i] = oBlocks.rStarts[1][j] + r2pos - r1pos;
            }
            r1_genome_len = r2pos - r1pos + r2_genome_len;
            break;
          }else{
            //cerr << "Fault with this synthetic read, outputting each of the overlapping reads as singles: " << read1->read_name << endl;
            // This error is not worth reporting. The current version of STAR outputs a good number of these, concordance would be nice, but it is better to get at least one read illustrating the splice junction.
            goodPair = false;
            oBlocks.readCount = 2;
          }
        }
    }

    if (!goodPair) {
	    oBlocks.readCount = 2;
    }
  }
	oBlocks.chr_id = read1->core.refID;
	oBlocks.readStart[0] = read1->core.pos;
	oBlocks.readEnd[0] = read1->core.pos + r1_genome_len;
	oBlocks.readName.resize(read1->core.l_read_name - 1);
	oBlocks.readName.replace(0, read1->core.l_read_name - 1, read1->read_name, read1->core.l_read_name - 1); // is this memory/speed efficient?

	unsigned int totalBlockLen = 0;
	for (auto blockLen: oBlocks.rLens[0]) {
		totalBlockLen += blockLen;
	}
	if (oBlocks.readCount > 1) {
		oBlocks.readStart[1] = read2->core.pos;
		oBlocks.readEnd[1] = read2->core.pos + r2_genome_len;
		for (auto blockLen: oBlocks.rLens[1]) {
			totalBlockLen += blockLen;
		}
	}
	//DEBUG:
	oBlocks.readName.append(debugstate);
	oBlocks.readName.append(to_string(oBlocks.readCount));
// TODO - restructure -- we could instead do the manipulation from 2 reads-> 1 synthetic in a non-const callback.
//        not required until that future flexibility is needed if part of the framework is repurposed.
	for (auto & callback : callbacksProcessBlocks ) {
		callback(oBlocks);
	}
	return totalBlockLen;
}


unsigned int BAM2blocks::processSingle(bam_read_core * read1) {
  int r1_genome_len;

  string debugstate;

  if (read1->core.flag & 0x10) {
    oBlocks.direction = 0;
  }else{
    oBlocks.direction = 1;
  }

  cigar2block(read1->cigar, read1->core.n_cigar_op, oBlocks.rStarts[0], oBlocks.rLens[0], r1_genome_len);
  oBlocks.readCount = 1;

	oBlocks.chr_id = read1->core.refID;
	oBlocks.readStart[0] = read1->core.pos;
	oBlocks.readEnd[0] = read1->core.pos + r1_genome_len;
	oBlocks.readName.resize(read1->core.l_read_name - 1);
	oBlocks.readName.replace(0, read1->core.l_read_name - 1, read1->read_name, read1->core.l_read_name - 1); // is this memory/speed efficient?
	//DEBUG:
	oBlocks.readName.append(debugstate);
	oBlocks.readName.append(to_string(oBlocks.readCount));
	//cout << "process pair - callbacks" << endl;  
	for (auto & callback : callbacksProcessBlocks ) {
		callback(oBlocks);
	}
	unsigned int totalBlockLen = 0;
	for (auto blockLen: oBlocks.rLens[0]) {
		totalBlockLen += blockLen;
	}
	return totalBlockLen;
}



int BAM2blocks::processAll(std::string& output, bool verbose) {

	unsigned long long totalNucleotides = 0;
	unsigned long j = 0;
	unsigned int idx = 0;
  int ret = 0;
	// int pair = 0;
	//int bytesread = 0;
	std::ostringstream oss;
#ifndef GALAXY
  uint64_t prev_bam_pos = IN->tellg();	// For progress bar
	Progress p(IN->GetLength(), verbose);
#endif

	// Use map pointer spare_reads:
  std::map< std::string, bam_read_core* > * spare_reads;
	std::map< std::string, bam_read_core* > * new_spare_reads;  
  
	spare_reads = new std::map< std::string, bam_read_core* >;
	int32_t spare_reads_count = 0;
	
	while(1) {
		ret = IN->read(reads[idx].c_block_size, 4);  // Should return 4 if all 4 bytes are read
		if (IN->eob()) {
      cErrorReads = spare_reads_count;
			oss << "Total reads processed\t" << j-1 << '\n';
			oss << "Total nucleotides\t" << totalNucleotides << '\n';
			oss << "Total singles processed\t" << cSingleReads << '\n';
			oss << "Total pairs processed\t" << cShortPairs+cIntersectPairs+cLongPairs << '\n';
			oss << "Short pairs\t" << cShortPairs << '\n';
			oss << "Intersect pairs\t" << cIntersectPairs << '\n';
			oss << "Long pairs\t" << cLongPairs << '\n';
			oss << "Skipped reads\t" << cSkippedReads << '\n';
			oss << "Chimeric reads\t" << cChimericReads << '\n';
			oss << "Error / Unpaired reads\t" << cErrorReads << '\n';
			oss << "Error detected on line\t" << "NA" << '\n';
      output = oss.str();
			return(0);   
		} else if(IN->fail() || (ret != 4)) {
      cErrorReads = spare_reads_count;
			// cerr << "Input error at line:" << j << ", return error (" << ret << ")\n";
			// cerr << "Characters read on last read call:" << IN->gcount() << '\n';
			oss << "Total reads processed\t" << j-1 << '\n';
			oss << "Total nucleotides\t" << totalNucleotides << '\n';
			oss << "Total singles processed\t" << cSingleReads << '\n';
			oss << "Total pairs processed\t" << cShortPairs+cIntersectPairs+cLongPairs << '\n';
			oss << "Short pairs\t" << cShortPairs << '\n';
			oss << "Intersect pairs\t" << cIntersectPairs << '\n';
			oss << "Long pairs\t" << cLongPairs << '\n';
			oss << "Skipped reads\t" << cSkippedReads << '\n';
			oss << "Chimeric reads\t" << cChimericReads << '\n';
			oss << "Error / Unpaired reads\t" << cErrorReads << '\n';
			oss << "Error detected on line\t" << j << '\n';
      output = oss.str();
			return(1);
			//This is possibly also just about the end of the file (say an extra null byte).
			//IN->gcount() knows how many characters were actually read last time.
		}
    if(reads[idx].block_size > BAM_READ_CORE_BYTES - 4) {
      j++;
      ret = IN->read(reads[idx].c, BAM_READ_CORE_BYTES - 4);
      ret = IN->read(reads[idx].read_name, reads[idx].core.l_read_name);
      ret = IN->read(reads[idx].cigar_buffer, reads[idx].core.n_cigar_op*4);    
      // debugs
      ret = IN->ignore(reads[idx].block_size - BAM_READ_CORE_BYTES + 4 - reads[idx].core.l_read_name - (reads[idx].core.n_cigar_op*4));

      if (reads[idx].core.flag & 0x904) {
        // If is an unmapped / secondary / supplementary alignment -- discard/overwrite
        cSkippedReads ++;
      }else if (! (reads[idx].core.flag & 0x1)) {
        // If is a single read -- process it as a single -- then discard/overwrite
        cSingleReads ++;
        totalNucleotides += processSingle(&reads[idx]);
      }else{
        if(idx == 0 && spare_reads->size() == 0) {
          // If BAM is sorted by read name, then we don't need read size, simply use old system
          idx++;
        } else if(idx == 1 && spare_reads->size() == 0 && 
            reads[0].core.l_read_name == reads[1].core.l_read_name &&
            (0 == strncmp(reads[0].read_name, reads[1].read_name, reads[1].core.l_read_name))) {
          cPairedReads ++;
          totalNucleotides += processPair(&reads[0], &reads[1]);
          idx = 0;
        } else {
          // Likely a coordinate sorted BAM file:
          for(unsigned int k = 0; k <= idx; k++) {
            std::string read_name = string(reads[k].read_name);
            auto it_read = spare_reads->find(read_name);
            
            if(it_read != spare_reads->end()){
              cPairedReads ++;
              if (reads[k].core.refID != it_read->second->core.refID) {
                cChimericReads += 1;
              } else {
                if (reads[k].core.pos <= it_read->second->core.pos) {
                  //cout << "procesPair call1" << endl;        
                  totalNucleotides += processPair(&reads[k], &(*(it_read->second)));
                  delete (it_read->second);
                  spare_reads->erase(read_name);
                  spare_reads_count -= 1;
                }else{
                  //cout << "procesPair call2" << endl;                
                  totalNucleotides += processPair(&(*(it_read->second)), &reads[k]);
                  delete (it_read->second);
                  spare_reads->erase(read_name);
                  spare_reads_count -= 1;
                }
              }
            } else {
              bam_read_core * store_read = new bam_read_core;
              *(store_read) = reads[k];
              spare_reads->insert({read_name, store_read});
              spare_reads_count += 1;
            }
          }
          idx = 0;
        }
      }
    }
#ifndef GALAXY
		if(verbose) {
			if ( (cPairedReads + cSingleReads) % 1000000 == 0 ) {
				p.increment((unsigned long)(IN->tellg() - prev_bam_pos));		
				prev_bam_pos = IN->tellg();				
			}
		}
#endif
		if ( (cPairedReads + cSingleReads) % 1000000 == 0 ) {
			// Clean map by swapping for a new one
			new_spare_reads = new std::map< std::string, bam_read_core* >;
			new_spare_reads->insert(spare_reads->begin(), spare_reads->end());
			spare_reads->swap(*new_spare_reads);
			delete new_spare_reads;
		}
	}
	delete spare_reads;
	return(0);
}



void BAM2blocks::openFile(BAMReader_Multi * _IN) {
	IN = _IN;
  readBamHeader(); // readBamHeader needs to call the ChrMappingChange callbacks.
}

void BAM2blocks::registerCallbackChrMappingChange( std::function<void(const std::vector<chr_entry> &)> callback ) {
	callbacksChrMappingChange.push_back(callback);
}

void BAM2blocks::registerCallbackProcessBlocks( std::function<void(const FragmentBlocks &)> callback ) {	
	callbacksProcessBlocks.push_back(callback);
}
