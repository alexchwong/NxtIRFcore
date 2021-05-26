#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include "includedefine.h"

class BAMReader {
	private:
    // char compressed_buffer[65536];
    // char buffer[65536];
    char * compressed_buffer;
    char * buffer;

    unsigned long bufferPos;
    unsigned long bufferMax;
    
    istream * IN;
    int IS_EOF;
    int IS_FAIL;
    size_t IS_LENGTH;
        
  public:
    BAMReader();
    ~BAMReader();
    bool eof();
    bool fail();
    uint64_t tellg();
    streamsize gcount();
    
    int read(char * dest, unsigned int len);
    int ignore(unsigned int len);
    int LoadBuffer();
    
    size_t GetLength() { return(IS_LENGTH); };
    
    void SetInputHandle(std::istream *in_stream);
    
};