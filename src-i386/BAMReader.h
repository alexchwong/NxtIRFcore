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
        
        // size_t EOF_POS;
        
 		static const int bamEOFlength = 28;
		static const char bamEOF[bamEOFlength+1];

		static const int bamGzipHeadLength = 16;  // +2 a uint16 with the full block length.
		static const char bamGzipHead[bamGzipHeadLength+1];

		union stream_uint32 {
			char c[4];
			uint32_t u;
		};
		union stream_uint16 {
			char c[2];
			uint16_t u;
		};
        
        // bool eof_check();
       
    public:
        BAMReader();
        bool eof();
        bool fail();
        uint64_t tellg();
        streamsize gcount();

        size_t IS_LENGTH;
        
        int read(char * dest, unsigned int len);
        int ignore(unsigned int len);
        int LoadBuffer();
        
        void SetInputHandle(std::istream *in_stream);
};
