#ifndef CODE_BAM2BLOCKS
#define CODE_BAM2BLOCKS

#include "FragmentBlocks.h"

#include "BAMReader.h"


/* Little Endian .. for big endian each group of 4 bytes needs to be reversed before individual members are accessed. */
// std c11 allows anonymous struct/union. -Wall may give a warning as non-portable to older c++ standards.



class BAM2blocks {
	static const int BAM_HEADER_BYTES = 8;
	static const int BAM_READ_CORE_BYTES = 36;
	static const int BAM_READ_CORE_MAX_CIGAR = 2000;

	FragmentBlocks oBlocks;

	std::vector< std::function<void(const std::vector<chr_entry> &)> > callbacksChrMappingChange;
	std::vector< std::function<void(const FragmentBlocks &)> > callbacksProcessBlocks;

	// Statistics.
	unsigned long cShortPairs;
	unsigned long cIntersectPairs;
	unsigned long cLongPairs;
	unsigned long cSingleReads;
	unsigned long cPairedReads;
	unsigned long cErrorReads;
	unsigned long cSkippedReads;
	unsigned long cChimericReads;


	bam_read_core reads[2];

	BAMReader * IN;

	void cigar2block(uint32_t * cigar, uint16_t n_cigar_op, std::vector<int> &starts, std::vector<int> &lens, int &ret_genome_len);

	unsigned int processPair(bam_read_core * read1, bam_read_core * read2);
	unsigned int processSingle(bam_read_core * read1);

  public:
  	BAM2blocks();
  	void openFile(BAMReader * _IN);
  	void readBamHeader();  // implied by openFile. So perhaps should be private.
  	int processAll(std::string& output, bool threaded = false);

    void registerCallbackChrMappingChange( std::function<void(const std::vector<chr_entry> &)> callback );
    void registerCallbackProcessBlocks( std::function<void(const FragmentBlocks &)> callback );

    std::string samHeader;
    std::vector<std::string> chr_names;   //tab terminated chromosome names.
    std::vector<int32_t> chr_lens;	//length of each chromosome (not used when reading, used if optionally outputting an altered BAM file)

    std::vector<chr_entry> chrs;
};


#endif
