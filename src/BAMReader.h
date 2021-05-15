#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include "includedefine.h"

class BAMReader {
	private:
    char compressed_buffer[65536];
    char buffer[65536];
    unsigned long bufferPos;
    unsigned long bufferMax;
    
    istream * IN;
    int IS_EOF;
    int IS_FAIL;

        
 		static const int bamEOFlength = 28;
		static const char bamEOF[bamEOFlength+1];

		static const int bamGzipHeadLength = 16;  // +2 a uint16 with the full block length.
		static const char bamGzipHead[bamGzipHeadLength+1];

  public:
    BAMReader();
    bool eof();
    bool fail();
    uint64_t tellg();
    streamsize gcount();

    
    int read(char * dest, unsigned int len);
    int ignore(unsigned int len);
    int LoadBuffer();
    
    void SetInputHandle(std::istream *in_stream);
    
    size_t IS_LENGTH;
};
