#include "includedefine.h"

class FastaReader {
private:
  istream * IN;
  
  bool FirstSeq;
public:
  std::string seqname;
  std::string sequence;
  
  void SetInputHandle(std::istream *in_stream);
  bool ReadSeq();
};