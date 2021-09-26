/* covReader.h Reads COV files (BAM coverage)

Copyright (C) 2021 Alex Chit Hei Wong

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.  */

#ifndef CODE_COVREADER
#define CODE_COVREADER

#include "covCommon.h"

class covReader {
	private:
    // Buffers
    char * compressed_buffer;
    char * buffer;

    unsigned long bufferPos;  // Position of decompressed buffer
    unsigned long bufferMax;  // Size of decompressed buffer
    
    uint32_t index_begin;     // File position of first byte of COV index
    uint32_t body_begin;      // File position of first byte of COV body
    
    istream * IN;
    
    int IS_EOF;               // Set to 1 if istream hits eof()
    int IS_FAIL;              // Set to 1 if istream hits fail()
    
    size_t IS_LENGTH;         // Size of COV file
    size_t EOF_POS;           // Position of first byte of BGZF EOF block
       
    std::vector<std::string> chr_names;   // seqnames
    std::vector<uint32_t> chr_lens;       // chromosome lengths

  public:
    covReader();
    ~covReader();
    void SetInputHandle(std::istream *in_stream);
    
    int ReadBuffer();
    int read(char * dest, unsigned int len);
    int ignore(unsigned int len);
    bool eof();
    bool fail();

    // Input functions
    
    int ReadHeader();
    int GetChrs(std::vector<chr_entry> &chrs);
    int FetchPos(const std::string seqname, const uint32_t start, const int strand,
      uint64_t * file_offset, uint32_t * block_start);
    int FetchRLE(const std::string seqname, 
      const uint32_t start, const uint32_t end, const int strand,
      std::vector<int> * values, std::vector<unsigned int> * lengths
    );
};

#endif