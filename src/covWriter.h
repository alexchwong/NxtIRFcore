/* covWriter.h Writes COV files (BAM coverage)

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

#ifndef CODE_COVWRITER
#define CODE_COVWRITER

#include "includedefine.h"

static const unsigned int BGZF_max = 65536 - 18 - 8;

// Body buffer for CovWriter
class buffer_out_chunk {
  private:
    char * buffer;
    char * compressed_buffer;
    
    // current write position of buffer
    unsigned int buffer_pos = 0;
    // number of bytes that need to be compressed
    unsigned int buffer_size = 0;   
    // number of bytes compressed that need to be written out
    unsigned int compressed_size = 0;   
    
  public:
    buffer_out_chunk();
    ~buffer_out_chunk();


    unsigned int write(char * src, unsigned int len);
    int WriteToFile(ostream * OUT);
    int Compress();
    
    unsigned int getBGZFSize() { return(compressed_size); };

    unsigned int SetPos(unsigned int pos) {
      if(pos >= BGZF_max) return(buffer_pos);
      buffer_pos = pos;
      if(pos > buffer_size) {
        buffer_size = pos;
      }
      return(pos);
    };

    unsigned int GetPos() { return(buffer_pos); };
    
    unsigned int write_to_pos(char * src, unsigned int len, unsigned int pos) {
      if(len + pos > BGZF_max) return(0);
      SetPos(pos);
      return(write(src, len));
    };
    
    bool IsAtCap(unsigned int len) {
      if(len + buffer_pos >= BGZF_max) return(true);
      return(false);
    }
};

class covWriter {
  private:
    ostream * OUT;
    
    // When chrs is set, initialize these:
    std::vector<chr_entry> chrs;
    
    // The buffers
    std::vector< std::vector<buffer_out_chunk> > body;          
    // The start coords of each bgzf
    std::vector< std::vector<uint32_t> > block_coord_starts;    

    int WriteEmptyEntry(unsigned int refID);
    int WriteHeaderToFile();
    int WriteIndexToFile();
  public:
    covWriter();
    ~covWriter();
    
    void SetOutputHandle(std::ostream *out_stream);
    
    int InitializeCOV(std::vector<chr_entry> chrs_to_copy);
  
    int WriteFragmentsMap(
      std::vector< std::pair<unsigned int, int> > * vec, 
      unsigned int chrID, unsigned int strand,
      unsigned int n_threads_to_use = 1
    );
    
    int WriteToFile();

};

#endif