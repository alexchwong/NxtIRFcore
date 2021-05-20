// WARNING: code is little endian only!

#include "BAM2blocks.h"
#include "includedefine.h"
// using namespace std;

BAM2blocks::BAM2blocks() {
	oBlocks = FragmentBlocks(); //Right syntax to call the default constructor on an object variable, declared but not initialised?

  cReadsProcessed = 0;
  totalNucleotides = 0;
	cShortPairs = 0;
	cIntersectPairs = 0;
	cLongPairs = 0;
	cSingleReads = 0;
	cPairedReads = 0;
	cErrorReads = 0;
	cSkippedReads = 0;
	cChimericReads = 0;
  
  spare_reads = new std::map< std::string, bam_read_core* >;
}

BAM2blocks::~BAM2blocks() {
  for(auto it = spare_reads->begin(); it != spare_reads->end(); it++) {
    delete it->second;
  }
  delete spare_reads;
}

unsigned int BAM2blocks::openFile(BAMReader_Multi * _IN, unsigned int n_workers) {
  // Change for multi-threading:
  // Only BBchild processes BAM reads from file input
  // BAM2blocks (parent) reads header, analyses file, and delegates tasks to children
  // - so if running 1 thread, have 1 BAM2blocks and 1 BBchild
  
  // All BBchild(s) share the same file handle
  // So in IRFinder main, only 1 child reads file buffer at any one time
  // Decompression and read processing runs in parallel
	IN = _IN;
  return(readBamHeader(block_begins, read_offsets, n_workers)); // readBamHeader needs to call the ChrMappingChange callbacks.
}

void BAM2blocks::AttachReader(BAMReader_Multi * _IN) {
	IN = _IN;
}

unsigned int BAM2blocks::readBamHeader(
    std::vector<uint64_t> &block_begins, 
    std::vector<unsigned int> &read_offsets,
    unsigned int n_workers
) {
  unsigned int n_bgzf = IN->readBamHeader(block_begins, read_offsets, n_workers);
  IN->fillChrs(chrs);
  return(n_bgzf);
}

void BAM2blocks::ProvideTask(unsigned int thread_number, 
        uint64_t &begin_bgzf, unsigned int &begin_pos,
        uint64_t &end_bgzf, unsigned int &end_pos
) {
  begin_bgzf = block_begins.at(thread_number);
  begin_pos = read_offsets.at(thread_number);
  end_bgzf = block_begins.at(thread_number + 1);
  end_pos = read_offsets.at(thread_number + 1);
}

void BAM2blocks::TransferChrs(BAM2blocks& other) {
  for(unsigned int i = 0; i < other.chrs.size(); i++) {
    chrs.push_back(other.chrs.at(i));
  }
	for (auto & callback : callbacksChrMappingChange ) {
		callback(chrs);
	}
  chrs_prepped = true;
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
  }else if (read1->core.pos + r1_genome_len >= read2->core.pos + r2_genome_len &&
      read1->core.pos < read2->core.pos){
    // if read1->core.pos == read2->core.pos, then order matters:
    // if r1_genome_len > r2_genome_len, will be short read, and intersect read if not
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

int BAM2blocks::WriteOutput(std::string& output) {
  std::ostringstream oss;
  cErrorReads = spare_reads->size();
  oss << "Total reads processed\t" << cReadsProcessed << '\n';
  oss << "Total nucleotides\t" << totalNucleotides << '\n';
  oss << "Total singles processed\t" << cSingleReads << '\n';
  oss << "Total pairs processed\t" << cShortPairs+cIntersectPairs+cLongPairs << '\n';
  oss << "Short pairs\t" << cShortPairs << '\n';
  oss << "Intersect pairs\t" << cIntersectPairs << '\n';
  oss << "Long pairs\t" << cLongPairs << '\n';
  oss << "Skipped reads\t" << cSkippedReads << '\n';
  oss << "Chimeric reads\t" << cChimericReads << '\n';
  oss << "Error / Unpaired reads\t" << cErrorReads << '\n';
  output = oss.str();
  return(0);
}

bam_read_core *  BAM2blocks::SupplyRead(std::string& read_name) {
  // Supplies the pointer to the last spare read
  // When called, transfer ownership of read to the parent BB
  if(spare_reads->size() == 0) return(NULL);
  auto it = spare_reads->begin();
  read_name = it->first;
  bam_read_core * read = it->second;
  spare_reads->erase(it);
  return(read);
}

int BAM2blocks::processSpares(BAM2blocks& other) {
  // Combines two BB's, and processes any matching paired reads
  cReadsProcessed += other.cReadsProcessed;
  totalNucleotides += other.totalNucleotides;
    
  cShortPairs += other.cShortPairs;
  cIntersectPairs += other.cIntersectPairs;
  cLongPairs += other.cLongPairs;
  cSingleReads += other.cSingleReads;
  cPairedReads += other.cPairedReads;
  cErrorReads += other.cErrorReads;
  cSkippedReads += other.cSkippedReads;
  cChimericReads += other.cChimericReads;
  
  while(1) {
    bam_read_core * spare_read;
    std::string read_name;
    spare_read = other.SupplyRead(read_name);
    
    if(!spare_read) {
      break;
    }
    
    auto it_read = spare_reads->find(read_name);
    if(it_read != spare_reads->end()){
      cPairedReads ++;
      if (spare_read->core.refID != it_read->second->core.refID) {
        cChimericReads += 1;
      } else {
        // Rcout << "Read matched: " << read_name << '\n';
        if (spare_read->core.pos <= it_read->second->core.pos) {    
          totalNucleotides += processPair(&(*spare_read), &(*(it_read->second)));
        } else{              
          totalNucleotides += processPair(&(*(it_read->second)), &(*spare_read));
        }
        cReadsProcessed+=2;
        delete (it_read->second);
        spare_reads->erase(read_name);
        delete spare_read;
      }
    } else {
      spare_reads->insert({read_name, spare_read});
    }
  }
  
  return(0);
}

int BAM2blocks::processAll() {
  // Reads from BAMReader until finished; do not create output
	unsigned int idx = 0;
  int ret = 0;
  bool any_reads_processed = false;
	// Use map pointer spare_reads:
	std::map< std::string, bam_read_core* > * new_spare_reads;  

	while(1) {
		ret = IN->read(reads[idx].c_block_size, 4);  // Should return 4 if all 4 bytes are read

		if (IN->eob() || IN->fail() || (ret != 4)) {
      // Bank the spare read:
      // Rcout << "BAM2blocks::processAll ret == " << ret << '\n';
      if(idx == 1 && spare_reads->size() == 0) {
        std::string read_name = string(reads[0].read_name);
        bam_read_core * store_read = new bam_read_core;
        *(store_read) = reads[0];
        spare_reads->insert({read_name, store_read});
      }
      cErrorReads = spare_reads->size();
      if(!any_reads_processed) return(1);
			return(0);   // This will happen if read fails - i.e. end of loaded buffer
		} else {
      any_reads_processed = true;
    }
    if(reads[idx].block_size > BAM_READ_CORE_BYTES - 4) {
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
        cReadsProcessed++;
      }else{
        if(idx == 0 && spare_reads->size() == 0) {
          // If BAM is sorted by read name, then we don't need read size, simply use old system
          idx++;
        } else if(idx == 1 && spare_reads->size() == 0 && 
            reads[0].core.l_read_name == reads[1].core.l_read_name &&
            (0 == strncmp(reads[0].read_name, reads[1].read_name, reads[1].core.l_read_name))) {
          cPairedReads ++;
          totalNucleotides += processPair(&reads[0], &reads[1]);
          cReadsProcessed+=2;
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
                }else{
                  //cout << "procesPair call2" << endl;                
                  totalNucleotides += processPair(&(*(it_read->second)), &reads[k]);
                }
                cReadsProcessed+=2;
                delete (it_read->second);
                spare_reads->erase(read_name);
              }
            } else {
              bam_read_core * store_read = new bam_read_core;
              *(store_read) = reads[k];
              spare_reads->insert({read_name, store_read});
            }
          }
          idx = 0;
        }
      }
    }

		if ( (cPairedReads + cSingleReads) % 1000000 == 0 ) {
			// Clean map by swapping for a new one
			new_spare_reads = new std::map< std::string, bam_read_core* >;
			new_spare_reads->insert(spare_reads->begin(), spare_reads->end());
			spare_reads->swap(*new_spare_reads);
			delete new_spare_reads;
		}
	}
	return(0);
}

void BAM2blocks::registerCallbackChrMappingChange( std::function<void(const std::vector<chr_entry> &)> callback ) {
	callbacksChrMappingChange.push_back(callback);
}

void BAM2blocks::registerCallbackProcessBlocks( std::function<void(const FragmentBlocks &)> callback ) {	
	callbacksProcessBlocks.push_back(callback);
}
