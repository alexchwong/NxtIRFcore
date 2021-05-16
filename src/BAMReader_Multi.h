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
  public:
    buffer_chunk();
    ~buffer_chunk();
    int clear_buffer();
    
    int read_from_file(istream * IN);   // Reads a bgzf from this stream into buffer, then sets max_buffer to non-zero
    int decompress();                   // Decompresses buffer, sets max_decompressed to non-zero

    bool is_decompressed() { return(max_decompressed != 0); };
    bool is_at_end() { return(max_decompressed != 0 && pos == max_decompressed); };
    
    unsigned int GetMaxBuffer() { return(max_buffer); };
    unsigned int GetMaxBufferDecompressed() { return(max_decompressed); };
    
    unsigned int read(char * dest, unsigned int len);  // returns the number of bytes actually read
    unsigned int ignore(unsigned int len);
};

class BAMReader_Multi {
	private:
    static const int n_bgzf = 1;
    istream * IN;
    int IS_EOF;
    int IS_EOB;   // End of file AND buffer
    int IS_FAIL;
    size_t IS_LENGTH;

    int n_threads = 1;
    
    bam_header Header;

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

    unsigned int read(char * dest, unsigned int len);  // returns the number of bytes actually read
    unsigned int ignore(unsigned int len);

    bool eof();
    bool eob();
    bool fail();
    uint64_t tellg();
    streamsize gcount();

    size_t GetLength() { return(IS_LENGTH); };
};
