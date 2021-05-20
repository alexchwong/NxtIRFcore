#ifndef CODE_BAM2BLOCKS
#define CODE_BAM2BLOCKS

#include "FragmentBlocks.h"

#include "BAMReader_Multi.h"


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

    unsigned int processPair(bam_read_core * read1, bam_read_core * read2);
    unsigned int processSingle(bam_read_core * read1);

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

    bool error_detected;

    bam_read_core reads[2];
    BAMReader_Multi * IN;  
    std::vector<chr_entry> chrs;
    bool chrs_prepped = false;

    std::vector<uint64_t> block_begins;
    std::vector<unsigned int> read_offsets;

    std::map< std::string, bam_read_core* > * spare_reads;    
  public:
  	BAM2blocks(); ~BAM2blocks();
  	unsigned int openFile(BAMReader_Multi * _IN, bool verbose = false,
      unsigned int n_workers = 1);
    void AttachReader(BAMReader_Multi * _IN);
    
    void ProvideTask(unsigned int thread_number, 
        uint64_t &begin_bgzf, unsigned int &begin_pos,
        uint64_t &end_bgzf, unsigned int &end_pos     
    );
    void TransferChrs(BAM2blocks& other);
  	int processAll();
  	int processSpares(BAM2blocks& other);

  	int WriteOutput(std::string& output);

    bam_read_core * SupplyRead(std::string& read_name);

    void registerCallbackChrMappingChange( std::function<void(const std::vector<chr_entry> &)> callback );
    void registerCallbackProcessBlocks( std::function<void(const FragmentBlocks &)> callback );

};


#endif
