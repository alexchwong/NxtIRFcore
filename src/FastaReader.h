#include "includedefine.h"

class FastaReader {
private:
  istream * IN;
  
  bool FirstSeq;
public:
  std::string seqname;
  std::string sequence;
  
  void SetInputHandle(std::istream *in_stream);
  void Profile();

  bool ReadSeq();

	std::vector<std::string> chr_names;   //tab terminated chromosome names.
	std::vector<int32_t> chr_lens;	//length of each chromosome (not used when reading, used if optionally outputting an altered BAM file)
};