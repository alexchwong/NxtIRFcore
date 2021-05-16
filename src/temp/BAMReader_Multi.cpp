#include <stddef.h>
#include "BAMReader_Multi.h"
#include <stdexcept>


buffer_chunk::buffer_chunk() {
  max_buffer = 0;
  max_decompressed = 0;
  pos = 0;
  buffer = NULL;
  decompressed_buffer = NULL;
}

buffer_chunk::~buffer_chunk() {
  if(buffer != NULL) free(buffer);
  if(decompressed_buffer != NULL) free(decompressed_buffer);
}

int buffer_chunk::read_from_file(istream * IN) {
  stream_uint16 u16;
  stream_uint32 u32;
  char GzipCheck[16];
  IN->read(GzipCheck, 16);

  if(IN->fail()) {
    // likely just EOF
    return(-1);
  } else if(IN->eof()) {
    IS_EOF = 1;
    return(1);
  }

  IN->read(u16.c, 2);
  max_buffer = u16.u + 1 - 2  - 16;
  
  buffer = (char*)malloc(max_buffer);
  IN->read(buffer, max_buffer);
  return(0);
}

int decompress() {
    stream_uint32 u32;

  if(max_decompressed == 0 && max_buffer > 10) {
    memcpy(u32.c, buffer + max_buffer - 8, 4);
    uint32_t crc_check = u32.u;
    
    max_decompressed = 65536;
    decompressed_buffer = (char*)malloc(max_decompressed);
    
    z_stream zs;
    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.msg = NULL;
    zs.next_in = (Bytef*)buffer;
    zs.avail_in = max_buffer;
    zs.next_out = (Bytef*)decompressed_buffer;
    zs.avail_out = max_decompressed;

    int ret = inflateInit2(&zs, -15);
    if(ret != Z_OK) {
        std::ostringstream oss;
        Rcout << "Exception during BAM decompression - inflateInit2() fail: (" << ret << ") ";
        return(ret);
    }
    ret = inflate(&zs, Z_FINISH);
    if(ret != Z_OK && ret != Z_STREAM_END) {
        std::ostringstream oss;
        Rcout << "Exception during BAM decompression - inflate() fail: (" << ret << ") ";
        return(ret);
    }
    ret = inflateEnd(&zs);
    
    max_decompressed -= zs.avail_out;
    // Don't really need to deallocate decompressed_buffer, as long as we know the real max

    uint32_t crc = crc32(crc32(0L, NULL, 0L), (Bytef*)decompressed_buffer, max_decompressed);
    if(crc_check != crc) {
        std::ostringstream oss;
        Rcout << "CRC fail during BAM decompression: (at " << IN->tellg() << " bytes) ";
        return(ret);
    }
    pos = 0;
    return(0);
  }
}

unsigned int buffer_chunk::read(char * dest, unsigned int len) {  
  if(!is_decompressed()) return(0);
  if(pos > max_decompressed) return(0);
  
  unsigned int bytes_to_read = min(max_decompressed, pos + len) - pos;
  memcpy(decompressed_buffer, decompressed_buffer + pos, bytes_to_read);
  return(bytes_to_read);
}

const char BAMReader_Multi::bamEOF[BAMReader_Multi::bamEOFlength+1] =
		"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00";
const char BAMReader_Multi::bamGzipHead[BAMReader_Multi::bamGzipHeadLength+1] = 
		"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00";

// Constructor
BAMReader_Multi::BAMReader_Multi() {
  
}

// Destructor
BAMReader_Multi::~BAMReader_Multi() {
  
}

int BAMReader_Multi::SetThreads(unsigned int &threads_to_use) {
	if(threads_to_use > 0 && threads_to_use <= omp_get_thread_limit()) {
    n_threads = threads_to_use;
	} else {
		n_threads = omp_get_thread_limit();
		if(n_threads < 1) {
			use_threads = 1;
		}
	}
  return(0);
}

void BAMReader_Multi::SetInputHandle(std::istream *in_stream) {
	IN = in_stream;
  //  get length of file:
  if(in_stream != &std::cin) {
    IN->seekg (0, std::ios_base::end);
    IS_LENGTH = IN->tellg();
    IN->seekg (0, std::ios_base::beg);    
  }  
}


unsigned int BAMReader_Multi::read(char * dest, unsigned int len) {  
  if(comp_buffer_count >= buffer_pos) {
    // Read from file
  } else if(buffer_count < comp_buffer_count) {
    // Decompress everything that's being read
  } else {
    // Read from current buffer
    
  }
}


bool BAMReader_Multi::eof() {
    if(IS_EOF == 1) {
        return (true);
    } else {
        if(IN->eof()) {
            IS_EOF = 1;
            return (true);
        } else {
            return (false);
        }
    }
}

uint64_t BAMReader_Multi::tellg() {
    return((uint64_t)IN->tellg());
}

bool BAMReader_Multi::fail() {
    return(IN->fail());
}


streamsize BAMReader_Multi::gcount() {
    return(IN->gcount());
}