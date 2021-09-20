/* covFile.h Implement COV format for BAM alignment coverage depth

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

#include "includedefine.h"

union stream_uint64 {
  char c[8];
  uint64_t u;
};
union stream_uint32 {
  char c[4];
  uint32_t u;
};
union stream_int32 {
  char c[4];
  int32_t i;
};
union stream_uint16 {
  char c[2];
  uint16_t u;
};

class covBuffer{
	private:
    char * compressed_buffer;
    char * buffer;

    unsigned long bufferPos;
    unsigned long bufferMax;
    uint64_t l_file_buffer;
    char * file_buffer;
    uint64_t file_bufferPos;
  public:
    covBuffer();
    ~covBuffer();
    bool BufferIsFull(unsigned int threshold = 8);
    int write(char * src, unsigned int len);
    int WriteBuffer();
    
    char * get_buffer_ptr() { return file_buffer; };
    uint64_t get_buffer_pos() { return file_bufferPos; };
};

class covFile {
	private:
    char * compressed_buffer;
    char * buffer;

    unsigned long bufferPos;
    unsigned long bufferMax;
    
    uint32_t index_begin;
    uint32_t body_begin;
    
    std::string mode;
    istream * IN;
    ostream * OUT;
    
    int IS_EOF;
    int IS_FAIL;
    
    size_t IS_LENGTH;
    size_t EOF_POS;
       
    unsigned int out_cur_seqID;
    char * chr_index;
    uint32_t chr_index_alloc;
    uint32_t chr_index_pos;
    uint32_t chr_coord;
    
    std::vector<std::string> chr_names;
    std::vector<unsigned int> chr_lens;

    covBuffer body;
  public:
    covFile();
    ~covFile();
    void SetOutputHandle(std::ostream *out_stream);
    void SetInputHandle(std::istream *in_stream);
    
    int ReadBuffer();
    int read(char * dest, unsigned int len);
    int write(char * src, unsigned int len);
    int WriteBuffer();
    int ignore(unsigned int len);
    bool eof();
    bool fail();

    // Input functions
    
    int ReadHeader();
    int FetchPos(const std::string seqname, const uint32_t start, const int strand,
      uint64_t * file_offset, uint32_t * block_start);
    int FetchRLE(const std::string seqname, const uint32_t start, const uint32_t end, const int strand,
                 std::vector<int> * values, std::vector<unsigned int> * lengths);    

    // Output functions
    
    int FlushBody();
    int WriteHeader(std::vector<std::string> s_chr, std::vector<int32_t> u_lens);
    int WriteEntry(unsigned int seqID, int value, unsigned int length);

    int GetChrs(std::vector<chr_entry> &chrs) {
      if(chr_names.size() > 0) {
        for(unsigned int i = 0; i < chr_names.size(); i++) {
          chrs.push_back(chr_entry(i, chr_names.at(i), chr_lens.at(i)));
        }
      }
      return(0);
    };

};

class buffer_out_chunk {
  private:
    static const int BUFFER_OUT_CAP = 65536 - 18 - 8;


    
    char * buffer;
    char * compressed_buffer;
    
    unsigned int buffer_pos = 0;
    unsigned int buffer_size = 0;   // number of bytes that need to be compressed
    
    unsigned int compressed_size = 0;   // number of bytes compressed that need to be written out
    
  public:
    buffer_out_chunk();
    ~buffer_out_chunk();

    unsigned int getBGZFSize() { return(compressed_size); };
    unsigned int write(char * src, unsigned int len);

    int Compress();
    int WriteToFile(ostream * OUT);
    
    unsigned int SetPos(unsigned int pos) {
      if(pos >= BUFFER_OUT_CAP) return(buffer_pos);
      buffer_pos = pos;
      if(pos > buffer_size) {
        buffer_size = pos;
      }
      return(pos);
    };
    unsigned int GetPos() { return(buffer_pos); };
    
    unsigned int write_to_pos(char * src, unsigned int len, unsigned int pos) {
      if(len + pos > BUFFER_OUT_CAP) return(0);
      SetPos(pos);
      return(write(src, len));
    };
    
    bool IsAtCap(unsigned int len) {
      if(len + buffer_pos >= BUFFER_OUT_CAP) return(true);
      return(false);
    }
};

class covWriter {
  private:
    ostream * OUT;
    
    std::vector<chr_entry> chrs;
    
    // When chrs is set, initialize these:
    std::vector< std::vector<buffer_out_chunk> > body;    // The buffers
    std::vector< std::vector<uint32_t> > block_coord_starts;    // The start coords of each bgzf
  public:
    covWriter();
    ~covWriter();
    
    void SetOutputHandle(std::ostream *out_stream);
    
    int WriteHeader(std::vector<chr_entry> chrs_to_copy);
  
    int WriteFragmentsMap(std::vector< std::pair<unsigned int, int> > * vec, 
        unsigned int chrID, unsigned int strand);
    int WriteEmptyEntry(unsigned int refID);
    int WriteEmptyEntry(unsigned int chrID, unsigned int strand) {
      return(WriteEmptyEntry(chrID + chrs.size() * strand)); };
    
    int WriteHeaderToFile();
    int WriteIndexToFile();
    int WriteToFile();
};