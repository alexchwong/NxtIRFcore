#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include "includedefine.h"

class buffer_chunk {
	private:
    unsigned int max_buffer;
    unsigned int max_decompressed;
    unsigned int pos;
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
    bool is_eof_block() { return(decompressed && max_decompressed == 0); }
    
    unsigned int GetMaxBuffer() { return(max_buffer); };
    unsigned int GetMaxBufferDecompressed() { return(max_decompressed); };
    
    unsigned int read(char * dest, unsigned int len);  // returns the number of bytes actually read
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
    
    int n_threads = 1;
 

    std::vector<buffer_chunk> buffer;

    unsigned int comp_buffer_count = 0;   // File reading will increase this count
    unsigned int buffer_count = 0;        // Multi-threaded decompress will increase this count
    unsigned int buffer_pos = 0;          // Reading past the current decompressed buffer will increase this count

    int read_from_file(unsigned int n_blocks);
    int decompress(unsigned int n_blocks);
  public:
    BAMReader_Multi();
    BAMReader_Multi(int threads_to_use);
    ~BAMReader_Multi();
    
    int SetThreads(int threads_to_use);
    void SetInputHandle(std::istream *in_stream);
    void readBamHeader();
    void fillChrs(std::vector<chr_entry> &chrs);
    void ProfileBAM(std::vector<uint64_t> &begin, std::vector<unsigned int> &first_read_offsets, int target_n_threads);
    
    unsigned int read(char * dest, unsigned int len);  // returns the number of bytes actually read
    unsigned int ignore(unsigned int len);

    bool eof();
    bool eob();
    bool fail();
    uint64_t tellg();
    streamsize gcount();

    size_t GetLength() { return(IS_LENGTH); };
    
    std::vector<chr_entry> chrs;
};
