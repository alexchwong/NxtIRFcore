#ifndef CODE_READBLOCKPROCESSOR
#define CODE_READBLOCKPROCESSOR

#include "FragmentBlocks.h"
#include "GZWriter.h"
#include "covFile.h"

/*
The code can be finished faster if we force a requirement that all input files are coordinate sorted by the start of each block.
ie: sort -k2,2n (for BED files).
Chromosome sorted or not won't matter, as these get split into different vectors in all cases.
*/



class ReadBlockProcessor {
	public:
		virtual void ProcessBlocks(const FragmentBlocks &fragblock) = 0;
		virtual void ChrMapUpdate(const std::vector<chr_entry> &chrmap) = 0; //Maybe some of these funcs shouldn't be pure virtual - overloadable if needed, but default often ok.
};

class JunctionCount : public ReadBlockProcessor {
	private:
		std::map<string, std::map<std::pair<unsigned int,unsigned int>,unsigned int[3]>> chrName_junc_count;
		std::vector<std::map<std::pair<unsigned int,unsigned int>,unsigned int[3]>*> chrID_junc_count;
		//unsigned int[3] - 0, neg strand count; 1, pos strand count; 2 = expected direction from ref: 0=unknown, 1=neg, 2=pos.

		std::map<string, std::map<unsigned int,unsigned int[2]>> chrName_juncLeft_count;
		std::vector<std::map<unsigned int,unsigned int[2]>*> chrID_juncLeft_count;

		std::map<string, std::map<unsigned int,unsigned int[2]>> chrName_juncRight_count;
		std::vector<std::map<unsigned int,unsigned int[2]>*> chrID_juncRight_count;
		  //chrID_... stores a fast access pointer to the appropriate structure in chrName_... 
	public:

		void ProcessBlocks(const FragmentBlocks &fragblock);
		void ChrMapUpdate(const std::vector<chr_entry> &chrmap);
		int WriteOutput(std::string& output, std::string& QC) const;
		void loadRef(std::istringstream &IN); //loadRef is optional, it allows directional detection to determine not just non-dir vs dir, but also which direction.

		int Directional(std::string& output) const;
		
		unsigned int lookup(std::string ChrName, unsigned int left, unsigned int right, bool direction) const;
		unsigned int lookup(std::string ChrName, unsigned int left, unsigned int right) const;
		unsigned int lookupLeft(std::string ChrName, unsigned int left, bool direction) const;
		unsigned int lookupLeft(std::string ChrName, unsigned int left) const;
		unsigned int lookupRight(std::string ChrName, unsigned int right, bool direction) const;
		unsigned int lookupRight(std::string ChrName, unsigned int right) const;

// Ideally we would read the XS junction strand attribute from the BAM if we want to count junctions from non-directional sequencing.
//   that will require BAM2blocks to be informed it should read the optional attributes looking for that attrib in that case.
// -- or we can just ignore direction -- the splice start/end information effectively determines the XS info (by ref to the reference)
};


class SpansPoint : public ReadBlockProcessor {
	private:
		std::map<string, std::vector<unsigned int>> chrName_pos;
		std::map<string, std::vector<unsigned int>> chrName_count[2];
		std::vector<std::vector<unsigned int>*> chrID_pos;
		std::vector<std::vector<unsigned int>*> chrID_count[2];
		char overhangLeft;
		char overhangRight;
		char overhangTotal;
		//chrID_... stores a fast access pointer to the appropriate structure in chrName_... 
	public:
		void setSpanLength(unsigned int overhang_left, unsigned int overhang_right);
		void loadRef(std::istringstream &IN);
		void ProcessBlocks(const FragmentBlocks &fragblock);
		void ChrMapUpdate(const std::vector<chr_entry> &chrmap);
		//void SetOutputStream(std::ostream *os);
		int WriteOutput(std::string& output, std::string& QC) const;
		unsigned int lookup(std::string ChrName, unsigned int pos, bool direction) const;
		unsigned int lookup(std::string ChrName, unsigned int pos) const;
};

class FragmentsInChr : public ReadBlockProcessor {
	// Counts the number of fragments in each Chromosome. (for both + & - strands).
	private:
		std::map<string, std::vector<unsigned int>> chrName_count; //only expecting 2 items in our vector.
		std::vector<std::vector<unsigned int>*> chrID_count;
	public:
		void ProcessBlocks(const FragmentBlocks &blocks);
		void ChrMapUpdate(const std::vector<chr_entry> &chrmap);
		int WriteOutput(std::string& output, std::string& QC) const;		
};


class FragmentsInROI : public ReadBlockProcessor {
	// Counts the number of fragments fully contained within a ROI.
	//   the ROIs may not overlap. Direction ignored for overlap detect.
	private:
		std::map<string, unsigned long> RegionID_counter[2];
 
		std::map<string, std::vector<std::pair<unsigned int,unsigned int>>> chrName_ROI;
		std::map<string, std::vector<unsigned long*>> chrName_count[2];

		std::vector<std::vector<std::pair<unsigned int,unsigned int>>*> chrID_ROI;
		std::vector<std::vector<unsigned long*>*> chrID_count[2];

		// Perhaps we want to store some text relating to each record too? Easy to do if the input is pre-sorted (at least within each Chr).
		//   if pre-sorted, it may be easier to check for no overlapping blocks on read .. or can do this immediately after read with a single nested-walk.
		std::map<string, std::vector<string>> chrName_ROI_text;
	public:
		void ProcessBlocks(const FragmentBlocks &blocks);
		void ChrMapUpdate(const std::vector<chr_entry> &chrmap);
		void loadRef(std::istringstream &IN);
		int WriteOutput(std::string& output, std::string& QC) const;		
};

class FragmentsMap : public ReadBlockProcessor {
  // Counts mappability.
private:
  // 0 = -, 1 = +, 2 = both
  std::vector< std::vector< std::pair<unsigned int, int> > > chrName_vec_final[3];
  std::vector< std::vector< std::pair<unsigned int, int> > > chrName_vec_new[3];
  std::vector< std::vector< std::pair<unsigned int, int> > > temp_chrName_vec_new[3];

  uint32_t frag_count = 0;
	int sort_and_collapse_temp();

	bool final_is_sorted = false;
  
  vector<chr_entry> chrs;
public:
  int sort_and_collapse_final(bool verbose);

  void ProcessBlocks(const FragmentBlocks &blocks);
  void ChrMapUpdate(const std::vector<chr_entry> &chrmap);
  int WriteOutput(std::ostream *os, int threshold = 4, bool verbose = false) ;
  int WriteBinary(covFile *os, bool verbose = false) ;

  void updateCoverageHist(std::map<unsigned int,unsigned int> &hist, unsigned int start, unsigned int end, unsigned int dir, const unsigned int &refID, bool debug = false) const;
};


/*
class CoverageBlocks : public ReadBlockProcessor { ... }
// In it's own file -- bigger code.
*/

#endif
