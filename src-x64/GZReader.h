// GZ File reader

#include <stdexcept>
#include "includedefine.h"

#define CHUNK_gz 262144

class GZReader {
private:
  gzFile gz_in;
  int GetBuffer();

public:
  GZReader();
  ~GZReader();
  int LoadGZ(std::string s_filename, bool asStream = false, bool lazy = false);
  int getline(std::string & s_myLine, const char delim);

  void read(char * dest, const unsigned long len);
  void ignore(const unsigned long len);
  bool eof();
  
  std::istringstream iss;
  char * buffer;
  unsigned long bufferLen;
  unsigned long bufferPos;
};
