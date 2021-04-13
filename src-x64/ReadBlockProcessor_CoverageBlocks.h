#ifndef CODE_READBLOCKPROCESSOR_COVERAGEBLOCKS
#define CODE_READBLOCKPROCESSOR_COVERAGEBLOCKS

#include "CoverageBlock.h"
#include "ReadBlockProcessor.h"
#include "FragmentBlocks.h"

struct BEDrecord {
	std::string chrName;
	std::string name;
	unsigned int start;
	unsigned int end;
	bool direction;
	
	std::vector<std::pair<unsigned int,unsigned int>> blocks;
};


class CoverageBlocks : public ReadBlockProcessor {
	//Store the Blocked BED record for each ROI/intron. This won't be referred to again until the end.
	//XX Create the temporary vectors (per Chr) which simply list the blocks sequentially as read.
	//XX Sort the temporary vectors
	//XX Build the final vectors of "blocks of interest"
	//xx Delete the temporary vectors
	//xx Create the parallel vectors with counter objects. (do these as a batch at the end, once vector size is known - for best memory layout)
	//xx Process fragments against the counter structure. (have I already written a class/object for this?)
	
	//Produce summary statistical output for each Blocked BED record, using the counter structure.

	private:

		// Coverage depth data-structures.
		std::map<string, std::vector<CoverageBlock>> chrName_CoverageBlocks;
		// Shortcut pointers to depth data-structures.
		std::vector<std::vector<CoverageBlock>*> chrID_CoverageBlocks;

		// TODO: what is optimal for speed & memory usage?
//		static const unsigned int coverage_block_max_length = 5000;
		static const unsigned int coverage_block_max_length = 500;

	protected:
		std::vector<BEDrecord> BEDrecords;


	public:
        ~CoverageBlocks();
		void ProcessBlocks(const FragmentBlocks &fragblock);
		void ChrMapUpdate(const std::vector<string> &chrmap);
		void loadRef(std::istringstream &IN);
		int WriteOutput(std::string& output) const;
		
		void fillHist(std::map<unsigned int,unsigned int> &hist, const std::string &chrName, const std::vector<std::pair<unsigned int,unsigned int>> &blocks) const;
		void fillHist(std::map<unsigned int,unsigned int> &hist, const std::string &chrName, const std::vector<std::pair<unsigned int,unsigned int>> &blocks, bool direction) const;

		double meanFromHist(const std::map<unsigned int,unsigned int> &hist) const;
		double coverageFromHist(const std::map<unsigned int,unsigned int> &hist) const;
		double percentileFromHist(const std::map<unsigned int,unsigned int> &hist, unsigned int percentile) const;
		double trimmedMeanFromHist(const std::map<unsigned int,unsigned int> &hist, unsigned int centerPercent) const;
};

class CoverageBlocksIRFinder : public CoverageBlocks {
	public:
		int WriteOutput(std::string& output, std::string& QC, const JunctionCount &JC, const SpansPoint &SP, int directionality = 0) const;
};


#endif
