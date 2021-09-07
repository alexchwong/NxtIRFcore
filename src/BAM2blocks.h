#ifndef CODE_BAM2BLOCKS
#define CODE_BAM2BLOCKS

#include "FragmentBlocks.h"

// #include "BAMReader_Multi.h"


/* Little Endian .. for big endian each group of 4 bytes needs to be reversed before individual members are accessed. */
// std c11 allows anonymous struct/union. -Wall may give a warning as non-portable to older c++ standards.



class BAM2blocks {
    static const int BAM_HEADER_BYTES = 8;
    static const int BAM_READ_CORE_BYTES = 36;
    static const int BAM_READ_CORE_MAX_CIGAR = 2000;

    FragmentBlocks oBlocks;

    std::vector< std::function<void(const std::vector<chr_entry> &)> > callbacksChrMappingChange;
    std::vector< std::function<void(const FragmentBlocks &)> > callbacksProcessBlocks;

    void cigar2block(uint32_t * cigar, uint16_t n_cigar_op, std::vector<int> &starts, std::vector<int> &lens, int &ret_genome_len);

    unsigned int processPair(pbam1_t * read1, pbam1_t * read2);
    unsigned int processSingle(pbam1_t * read1, bool mappability_mode = false);

  	unsigned int readBamHeader(
      std::vector<uint64_t> &block_begins, 
      std::vector<unsigned int> &read_offsets, bool verbose = false,
      unsigned int n_workers = 1);  // implied by openFile. So perhaps should be private.

    // Statistics.
    unsigned long cReadsProcessed;
    unsigned long long totalNucleotides;
    
    unsigned long cShortPairs;
    unsigned long cIntersectPairs;
    unsigned long cLongPairs;
    unsigned long cSingleReads;
    unsigned long cPairedReads;
    unsigned long cErrorReads;
    unsigned long cSkippedReads;
    unsigned long cChimericReads;

    // bool error_detected;

    pbam1_t reads[2];
    pbam_in * IN;
    
    std::vector<chr_entry> chrs;
    std::vector<std::string> BB_ref_names;
    std::vector<std::string> BB_ref_alias;

    std::vector<uint64_t> block_begins;
    std::vector<unsigned int> read_offsets;

    std::map< std::string, pbam1_t* > * spare_reads;
    pbam1_t * SupplyRead(std::string& read_name);    
    int realizeSpareReads();
  public:
  	BAM2blocks();
  	BAM2blocks(
      std::vector<std::string> & ref_names, 
      std::vector<std::string> & ref_alias
    );
    ~BAM2blocks();
  	unsigned int openFile(pbam_in * _IN);

  	int processAll(unsigned int thread_number = 0, bool mappability_mode = false);
  	int processSpares(BAM2blocks& other);

  	int WriteOutput(std::string& output);

    void registerCallbackChrMappingChange( std::function<void(const std::vector<chr_entry> &)> callback );
    void registerCallbackProcessBlocks( std::function<void(const FragmentBlocks &)> callback );

};


#endif
