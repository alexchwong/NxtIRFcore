#ifndef CODE_BAMREADER_MULTI
#define CODE_BAMREADER_MULTI

#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include "includedefine.h"

class buffer_chunk {
	private:
    size_t bgzf_pos;
    unsigned int max_buffer;
    unsigned int max_decompressed;
    unsigned int pos;
    unsigned int end_pos;
    char * buffer;
    char * decompressed_buffer;
    bool decompressed;
  public:
    buffer_chunk();
    ~buffer_chunk();
    int clear_buffer();
    
    int read_from_file(istream * IN);   // Reads a bgzf from this stream into buffer, then sets max_buffer to non-zero
    int decompress();                   // Decompresses buffer, sets max_decompressed to non-zero

    bool is_decompressed() { return(decompressed); };
    bool is_at_end() { return(decompressed && pos == max_decompressed); };
    bool is_eof_block() { 
      return(max_buffer == 10);
      // return(decompressed && max_decompressed == 0); 
    }
    
    unsigned int GetPos() { return(pos); };
    unsigned int GetRemainingBytes() { return(max_decompressed - pos); };
    unsigned int GetBGZFPos() { return(bgzf_pos); };
    unsigned int SetPos(unsigned int pos_to_set) {
      pos = pos_to_set;
      return(pos); 
    };
    unsigned int SetBGZFPos(size_t pos_to_set) { 
      bgzf_pos = pos_to_set;
      return(bgzf_pos); 
    };

    unsigned int SetEndPos(unsigned int pos_to_set) {
      end_pos = pos_to_set;
      return(end_pos);
    };
    
    unsigned int GetMaxBuffer() { return(max_buffer); };
    unsigned int GetMaxBufferDecompressed() { return(max_decompressed); };
    
    unsigned int read(char * dest, unsigned int len);  // returns the number of bytes actually read
    unsigned int peek(char * dest, unsigned int len);  // same as read but does not move cursor
    unsigned int ignore(unsigned int len);
};

class BAMReader_Multi {
	private:
    // Each thread will load n_bgzf BGZF blocks; 100 BGZF blocks ~ 6 Mb decompressed
    int n_bgzf = 1;  
    static const int BAM_HEADER_BYTES = 8;
    static const int BAM_READ_CORE_BYTES = 36;
    
    istream * IN;
    int IS_EOF;
    int IS_EOB;   // End of file AND buffer
    int IS_FAIL;
    size_t IS_LENGTH;
    size_t BAM_READS_BEGIN;
    
    size_t BAM_BLOCK_CURSOR;
    
    std::vector<buffer_chunk> buffer;

    unsigned int comp_buffer_count = 0;   // File reading will increase this count
    unsigned int buffer_count = 0;        // Multi-threaded decompress will increase this count
    unsigned int buffer_pos = 0;          // Reading past the current decompressed buffer will increase this count

    // Allows setting BAMReader_Multi with rules to read a subset of BAM:
    // These rules are controlled when relevant BGZF block is read from file
    bool auto_load_data = true;
    uint64_t begin_block_offset = 0;
    unsigned int begin_read_offset = 0;
    uint64_t end_block_offset = 0;
    unsigned int end_read_offset = 65536;

  public:
    BAMReader_Multi();
    BAMReader_Multi(uint64_t block_begin, unsigned int begin_offset,
      uint64_t block_end, unsigned int end_offset);
    ~BAMReader_Multi();

    void AssignTask(std::istream *in_stream, uint64_t block_begin, unsigned int begin_offset,
      uint64_t block_end, unsigned int end_offset);
    void SetAutoLoad(bool autoload) {auto_load_data = autoload;};
    
    void SetInputHandle(std::istream *in_stream);
    std::istream* GetFileHandle() { return(IN); };
    unsigned int readBamHeader(
      std::vector<uint64_t> &block_begins, 
      std::vector<unsigned int> &read_offsets, bool verbose,
      unsigned int n_workers = 1);
    void fillChrs(std::vector<chr_entry> &chrs);
    unsigned int ProfileBAM(std::vector<uint64_t> &block_begins, 
      std::vector<unsigned int> &read_offsets, bool verbose,
      unsigned int target_n_threads = 1);

    int getBGZFstarts(std::vector<uint64_t> & BGZF_begins, bool verbose = false);
#ifdef _OPENMP
    int getBGZFstarts_OpenMP(std::vector<uint64_t> & BGZF_begins, bool verbose = false, unsigned int n_threads = 1);
#endif
    int read_from_file(unsigned int n_blocks);
    int decompress(bool allow_openmp = false);
    
    unsigned int read(char * dest, unsigned int len);  // returns the number of bytes actually read
    unsigned int ignore(unsigned int len);
    unsigned int peek(char * dest, unsigned int len) ;
    
    bool isReadable();
    bool isReadableStrict();
    bool GotoNextRead(bool strict = false);
    bool eof();
    bool eob();
    bool fail() {return(IN->fail());};
    uint64_t tellg() {return((uint64_t)IN->tellg());};
    streamsize gcount() {return(IN->gcount());};
    size_t GetLength() { return(IS_LENGTH); };
    
    std::vector<chr_entry> chrs;
};

#endif