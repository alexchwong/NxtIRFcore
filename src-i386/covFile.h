#include "includedefine.h"

class covBuffer{
	private:
    static const int bamGzipHeadLength = 16;  // +2 a uint16 with the full block length.
    static const char bamGzipHead[bamGzipHeadLength+1];
    
    char compressed_buffer[65536];
    char buffer[65536];

    unsigned long bufferPos;
    unsigned long bufferMax;

	union stream_uint32 {
		char c[4];
		uint32_t u;
	};
	union stream_uint16 {
		char c[2];
		uint16_t u;
	};

  public:
    covBuffer();
    ~covBuffer();
    bool BufferIsFull(unsigned int threshold = 8);
    int write(char * src, unsigned int len);
    int WriteBuffer();
    
    uint64_t l_file_buffer;
    char * file_buffer;
    uint64_t file_bufferPos;
};

class covFile {
	private:
 		static const int bamEOFlength = 28;
		static const char bamEOF[bamEOFlength+1];

		static const int bamGzipHeadLength = 16;  // +2 a uint16 with the full block length.
		static const char bamGzipHead[bamGzipHeadLength+1];

    char compressed_buffer[65536];
    char buffer[65536];

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
        
	union stream_int32 {
		char c[4];
		int32_t i;
	};

	union stream_uint64 {
		char c[8];
		uint64_t u;
	};
	union stream_uint32 {
		char c[4];
		uint32_t u;
	};
	union stream_uint16 {
		char c[2];
		uint16_t u;
	};
        
    unsigned int out_cur_seqID;
    char * chr_index;
    uint32_t chr_index_alloc;
    uint32_t chr_index_pos;
    uint32_t chr_coord;
        
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

    std::vector<std::string> chr_names;
    std::vector<unsigned int> chr_lens;
    
    covBuffer body;
};
