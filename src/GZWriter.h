// GZ File writer

#include "includedefine.h"

#define CHUNK_gz 262144

class GZWriter {
private:
  ostream * OUT;
//  gzFile gz_out;
  
  char compressed_buffer[CHUNK_gz];
  
  char buffer[CHUNK_gz];
  unsigned int bufferPos = 0;
  
public:
  void SetOutputHandle(std::ostream *out_stream);

  int writebuffer(const char * src, unsigned int len);
  int writeline(const std::string& s_src);
  int writestring(const std::string& s_src);
  int flush(bool final = false);
};
