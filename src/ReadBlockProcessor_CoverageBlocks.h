#ifndef CODE_READBLOCKPROCESSOR_COVERAGEBLOCKS
#define CODE_READBLOCKPROCESSOR_COVERAGEBLOCKS

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

  
	protected:
		std::vector<BEDrecord> BEDrecords;

	public:
		void ProcessBlocks(const FragmentBlocks &fragblock);
		void ChrMapUpdate(const std::vector<chr_entry> &chrmap);
		void loadRef(std::istringstream &IN);
		int WriteOutput(std::string& output, const FragmentsMap &FM) const;
		
	  void fillHist(std::map<unsigned int,unsigned int> &hist, const unsigned int &refID, const std::vector<std::pair<unsigned int,unsigned int>> &blocks, const FragmentsMap &FM, bool debug = false) const;
		void fillHist(std::map<unsigned int,unsigned int> &hist, const unsigned int &refID, const std::vector<std::pair<unsigned int,unsigned int>> &blocks, bool direction, const FragmentsMap &FM, bool debug = false) const;

		double meanFromHist(const std::map<unsigned int,unsigned int> &hist) const;
		double coverageFromHist(const std::map<unsigned int,unsigned int> &hist) const;
		double percentileFromHist(const std::map<unsigned int,unsigned int> &hist, unsigned int percentile) const;
		double trimmedMeanFromHist(const std::map<unsigned int,unsigned int> &hist, unsigned int centerPercent, bool debug = false) const;

    vector<chr_entry> chrs;
};

class CoverageBlocksIRFinder : public CoverageBlocks {
	public:
		void Combine(CoverageBlocksIRFinder &child);
		int WriteOutput(std::string& output, std::string& QC, const JunctionCount &JC, const SpansPoint &SP, const FragmentsMap &FM, int directionality = 0) const;
};


#endif
