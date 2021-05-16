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
  if(buffer) free(buffer);
  if(decompressed_buffer) free(decompressed_buffer);
}

int buffer_chunk::clear_buffer() {
  if(buffer) free(buffer);
  if(decompressed_buffer) free(decompressed_buffer);
  buffer = NULL;
  decompressed_buffer = NULL;
  max_buffer = 0;
  max_decompressed = 0;
  pos = 0;
  return(0);
}

int buffer_chunk::read_from_file(istream * IN) {
  stream_uint16 u16;

  char GzipCheck[16];
  IN->read(GzipCheck, 16);

  if(IN->fail()) {
    // likely just EOF
    return(-1);
  } else if(IN->eof()) {
    // IS_EOF = 1;
    return(1);
  }

  IN->read(u16.c, 2);
  max_buffer = u16.u + 1 - 2  - 16;
  Rcout << "BGZF block " << max_buffer << " bytes";
  buffer = (char*)malloc(max_buffer + 1);
  IN->read(buffer, max_buffer);
  Rcout << " read\n";

  return(0);
}

int buffer_chunk::decompress() {
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

    Rcout << "BGZF block " << max_decompressed << " bytes decompressed\n";

    uint32_t crc = crc32(crc32(0L, NULL, 0L), (Bytef*)decompressed_buffer, max_decompressed);
    if(crc_check != crc) {
        std::ostringstream oss;
        Rcout << "CRC fail during BAM decompression";
        return(ret);
    }
    pos = 0;
  }
  return(0);
}

unsigned int buffer_chunk::read(char * dest, unsigned int len) {  
  if(!is_decompressed()) return(0);
  if(pos >= max_decompressed) return(0);
  
  unsigned int bytes_to_read = min(max_decompressed, pos + len) - pos;
  memcpy(dest, decompressed_buffer + pos, bytes_to_read);
  pos += bytes_to_read;
  return(bytes_to_read);
}

unsigned int buffer_chunk::ignore(unsigned int len) {  
  if(!is_decompressed()) return(0);
  if(pos >= max_decompressed) return(0);

  unsigned int bytes_to_read = min(max_decompressed, pos + len) - pos;
  // memcpy(decompressed_buffer, decompressed_buffer + pos, bytes_to_read);
  pos += bytes_to_read;
  return(bytes_to_read);  
}

// Constructor
BAMReader_Multi::BAMReader_Multi() {
  IS_EOF = 0;
  IS_EOB = 0;   // End of file AND buffer
  IS_FAIL = 0;
  IS_LENGTH = 0;
  IN = NULL;
}

BAMReader_Multi::BAMReader_Multi(int threads_to_use) {
  IS_EOF = 0;
  IS_EOB = 0;   // End of file AND buffer
  IS_FAIL = 0;
  IS_LENGTH = 0;
  IN = NULL;
  SetThreads(threads_to_use);
}

// Destructor
BAMReader_Multi::~BAMReader_Multi() {
  if(buffer.size() > 0) {
    for(unsigned int i = 0; i < buffer.size(); i++) {
      buffer.at(i).clear_buffer();
    }
  }
}

int BAMReader_Multi::SetThreads(int threads_to_use) {
	if(threads_to_use > 0 && threads_to_use <= omp_get_thread_limit()) {
    n_threads = threads_to_use;
	} else {
		n_threads = omp_get_thread_limit();
		if(n_threads < 1) {
			n_threads = 1;
		}
	}
  omp_set_num_threads(n_threads);
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
  
  // Prime the pump:
  read_from_file(n_bgzf * n_threads);
  decompress(n_bgzf * n_threads);
}

int BAMReader_Multi::read_from_file(unsigned int n_blocks) {
  if(comp_buffer_count > 0) {
    if(buffer.at(comp_buffer_count - 1).GetMaxBuffer() == 10) {
      IS_EOF = 1;
      return(Z_STREAM_END);
    }
  }
  unsigned int i = 0;
  while(i < n_blocks) {
    buffer.push_back(buffer_chunk());
    int ret = buffer.at(comp_buffer_count).read_from_file(IN);
    if(ret != 0) return(ret);
    comp_buffer_count++;
    if(buffer.at(comp_buffer_count - 1).GetMaxBuffer() == 10) {
      IS_EOF = 1;
      return(Z_STREAM_END);
    }
    i++;
  }
  return(0);
}

int BAMReader_Multi::decompress(unsigned int n_blocks) {
  unsigned int end_blocks = min(comp_buffer_count, n_blocks + buffer_count);
  
  if(n_threads > 1) {
    #pragma omp parallel for
    for(unsigned int i = buffer_count; i < end_blocks; i++) {
      if(!buffer.at(i).is_decompressed()) buffer.at(i).decompress();
    }
  } else {
    for(unsigned int i = buffer_count; i < end_blocks; i++) {
      if(!buffer.at(i).is_decompressed()) buffer.at(i).decompress();
    }
  }
  buffer_count = end_blocks;
  Rcout << "BAMReader_Multi " << n_blocks << " decompressed\n";

  return(0);
}

unsigned int BAMReader_Multi::read(char * dest, unsigned int len) {  
  // Read from current buffer
  unsigned int cursor = 0;
  while(cursor < len) {
    if(buffer.at(buffer_pos).is_at_end()) {
      // destroy current buffer
      buffer.at(buffer_pos).clear_buffer();
      if(IS_EOF == 1) {
        IS_EOB = 1;
        return(cursor);
      }
      buffer_pos++;
    }
    if(buffer_pos == buffer.size() && IS_EOF != 1) {
      read_from_file(n_bgzf * n_threads);
      decompress(n_bgzf * n_threads);
    }
    cursor += buffer.at(buffer_pos).read(dest + cursor, len - cursor);
  }
  if(cursor < len) {
    IS_EOB = 1;
  } else if(IS_EOF == 1 && buffer.at(buffer_pos).is_at_end()) {
    buffer.at(buffer_pos).clear_buffer();
    IS_EOB = 1;
  }
  return(cursor);
}

unsigned int BAMReader_Multi::ignore(unsigned int len) {  
  // Read from current buffer
  unsigned int cursor = 0;
  while(cursor < len) {
    if(buffer.at(buffer_pos).is_at_end()) {
      // destroy current buffer
      buffer.at(buffer_pos).clear_buffer();
      buffer_pos++;
    }
    if(buffer_count < buffer_pos) {
      if(comp_buffer_count < buffer_pos) {
        if(IS_EOF == 1) break;
        read_from_file(10 * n_threads);
      }
      if(buffer_count < comp_buffer_count) {
        decompress(10 * n_threads);
      }
    }
    cursor += buffer.at(buffer_pos).ignore(len - cursor);
  }
  if(cursor < len) {
    IS_EOB = 1;
  } else if(IS_EOF == 1 && buffer.at(buffer_pos).is_at_end()) {
    buffer.at(buffer_pos).clear_buffer();
    IS_EOB = 1;
  }
  return(cursor);
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

// End of Buffer is when both file and buffer are exhausted
bool BAMReader_Multi::eob() {
  if(IS_EOB == 1) {
      return (true);
  } else {
    if(IS_EOF == 1 && buffer.at(buffer_pos).is_at_end()) {
      buffer.at(buffer_pos).clear_buffer();
      IS_EOB = 1;
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