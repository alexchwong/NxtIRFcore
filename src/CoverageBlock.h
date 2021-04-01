#ifndef CODE_COVERAGEBLOCK
#define CODE_COVERAGEBLOCK

#include "includedefine.h"

class start_stops {
	public:
		unsigned char start[2];
		unsigned char end[2];

		start_stops() {
			start[0]=0;
			start[1]=0;
			end[0]=0;
			end[1]=0;
		};
};

class start_stopsL {
	public:
		unsigned int start[2];
		unsigned int end[2];

		start_stopsL() {
			start[0]=0;
			start[1]=0;
			end[0]=0;
			end[1]=0;
		};
		start_stopsL(const start_stops &copy) {
			start[0]=copy.start[0];
			start[1]=copy.start[1];
			end[0]=copy.end[0];
			end[1]=copy.end[1];
		};

};


class CoverageBlock {
	private:
	    unsigned int blockStart;
	    unsigned int blockEnd;
	    unsigned int firstDepth[2];

		std::vector<start_stops>* blockExtents;
		std::vector<start_stopsL>* blockExtentsL;

		inline unsigned int vectorLen() {
			return (blockEnd - blockStart - 1);
		};
	public:
		CoverageBlock(unsigned int start, unsigned int end);
        ~CoverageBlock();
		void RecordCover(unsigned int start, unsigned int end, bool dir);
		//RetrieveCover(..);
		void print(std::ostream& os) const;
		
		//First form, non-directional. Second form, directional with "dir" specifiying whether sense or anti-sense.
		void updateCoverageHist(std::map<unsigned int,unsigned int> &hist, unsigned int start, unsigned int end) const;
		void updateCoverageHist(std::map<unsigned int,unsigned int> &hist, unsigned int start, unsigned int end, bool dir) const;

		inline bool posIsAfterStart(const unsigned int &compareval) const {
		  return (compareval > blockStart);
		};

		// http://www.learncpp.com/cpp-tutorial/94-overloading-the-comparison-operators/
		// http://en.cppreference.com/w/cpp/language/operator_comparison
		inline bool operator<(const CoverageBlock &b) const {
			return (blockEnd < b.blockEnd);
		};
		inline bool operator<(const unsigned int &b) const {
			return (blockEnd < b);  //a is the object.
		};
		friend inline bool operator<(const unsigned int &a, const CoverageBlock &b) {
			return (a < b.blockEnd);  //a is a uint.
		};
};

std::ostream& operator<<( std::ostream& os, const CoverageBlock& cb);

#endif
