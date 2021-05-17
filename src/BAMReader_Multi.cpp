#include <stddef.h>
#include "BAMReader_Multi.h"
#include <stdexcept>


buffer_chunk::buffer_chunk() {
  max_buffer = 0;
  max_decompressed = 0;
  pos = 0;
  buffer = NULL;
  decompressed_buffer = NULL;
  decompressed = false;
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
  decompressed = false;
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
  // Rcout << "BGZF block " << max_buffer << " bytes";
  buffer = (char*)malloc(max_buffer + 1);
  IN->read(buffer, max_buffer);
  // Rcout << " read\n";

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
    decompressed = true;
    // Don't really need to deallocate decompressed_buffer, as long as we know the real max

    // Rcout << "BGZF block " << max_decompressed << " bytes decompressed\n";

    uint32_t crc = crc32(crc32(0L, NULL, 0L), (Bytef*)decompressed_buffer, max_decompressed);
    if(crc_check != crc) {
        std::ostringstream oss;
        Rcout << "CRC fail during BAM decompression";
        return(ret);
    }
    pos = 0;
  } else {
    decompressed = true;
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
  BAM_READS_BEGIN = 0;
  IN = NULL;
}

BAMReader_Multi::BAMReader_Multi(int threads_to_use) {
  IS_EOF = 0;
  IS_EOB = 0;   // End of file AND buffer
  IS_FAIL = 0;
  IS_LENGTH = 0;
  BAM_READS_BEGIN = 0;
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
#ifndef _OPENMP
  if(threads_to_use > 1) Rcout << "OpenMP not built for this system... running in single thread\n";
  n_threads = 1
#else
	if(threads_to_use > 0 && threads_to_use <= omp_get_thread_limit()) {
    n_threads = threads_to_use;
	} else {
		n_threads = omp_get_thread_limit();
		if(n_threads < 1) {
			n_threads = 1;
		}
	}
  omp_set_num_threads(n_threads);
#endif  
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

// OK.
void BAMReader_Multi::readBamHeader() {
  char buffer[1000];
  std::string chrName;

  bam_header bamhead;
  read(bamhead.c, BAM_HEADER_BYTES);

  char * headertext = new char[bamhead.magic.l_text+1];
  read(headertext, bamhead.magic.l_text);
  
  std::string samHeader = string(headertext, bamhead.magic.l_text);
  delete[] headertext;
  
  stream_int32 i32;
  read(i32.c ,4);
  unsigned int n_chr = i32.i;

  for (unsigned int i = 0; i < n_chr; i++) {
    read(i32.c ,4);
    read(buffer , i32.i);
    chrName = string(buffer, i32.i-1);
    read(i32.c ,4);
    chrs.push_back(chr_entry(i, chrName, i32.i));
  }
  std::sort(chrs.begin(), chrs.end());
  // Rcout << "BAM cursor at" << tellg() << '\n';
  if(buffer_pos == comp_buffer_count) BAM_READS_BEGIN = tellg();
}

void BAMReader_Multi::fillChrs(std::vector<chr_entry> &chrs_dest) {
  for(unsigned int i = 0; i < chrs.size(); i++) {
    chrs_dest.push_back(chrs.at(i));
  }
}

void BAMReader_Multi::ProfileBAM(
    std::vector<uint64_t> &begin, std::vector<unsigned int> &first_read_offsets, int target_n_threads) {
  if(BAM_READS_BEGIN > 0) {
    std::vector<uint64_t> begins;
    std::vector<unsigned int> offsets;
    stream_uint32 u32;
    unsigned int bytes_read;
    
    bool end_of_buffer = false;
    bool break_at_read_head = false;
    bool break_at_read_body = false;
    
    
    begins.push_back(BAM_READS_BEGIN);
    offsets.push_back(0);
    IN->seekg (BAM_READS_BEGIN, std::ios_base::beg);
    
    while(!eof()) {
      buffer_chunk * temp_buffer = new buffer_chunk;
      temp_buffer->read_from_file(IN);
      temp_buffer->decompress();
      
      end_of_buffer = false;
      break_at_read_head = false;
      break_at_read_body = false;
      while(!end_of_buffer) {
        bytes_read = temp_buffer->read(u32.c, 4);
        if(bytes_read < 4) {
          end_of_buffer = true;
          break_at_read_head = true;
          break;
        }
        bytes_read = temp_buffer->ignore(u32.u - 4);
        if(bytes_read < u32.u - 4) {
          end_of_buffer = true;
          break_at_read_body = true;
          break;
        }
        if(temp_buffer->is_at_end()) {
          end_of_buffer = true;
        }
      }
      delete temp_buffer;
    }
  }
}

int BAMReader_Multi::read_from_file(unsigned int n_blocks) {
  unsigned int i = 0;
  buffer.resize(comp_buffer_count + n_blocks);
  while(i < n_blocks) {
    int ret = buffer.at(comp_buffer_count).read_from_file(IN);
    if(ret != 0) {
      if(ret == 1) {
        IS_EOF = 1;
      } else {
        Rcout << "Error reading file, error code: " << ret << '\n';
      }
      return(ret);
    }
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
#ifndef _OPENMP
#else
    #pragma omp parallel for
#endif
    for(unsigned int i = buffer_count; i < end_blocks; i++) {
      if(!buffer.at(i).is_decompressed()) buffer.at(i).decompress();
    }
  } else {
    for(unsigned int i = buffer_count; i < end_blocks; i++) {
      if(!buffer.at(i).is_decompressed()) buffer.at(i).decompress();
    }
  }
  buffer_count = end_blocks;
  // Rcout << "BAMReader_Multi " << n_blocks << " decompressed\n";

  return(0);
}

unsigned int BAMReader_Multi::read(char * dest, unsigned int len) {  
  // Read from current buffer
  unsigned int cursor = 0;
  if(IS_EOB == 1) return(cursor);
  while(cursor < len) {
    if(buffer_pos == comp_buffer_count && IS_EOF != 1) {
      read_from_file(n_bgzf * n_threads);
      decompress(n_bgzf * n_threads);
    } // reading will always start with reading buffer if current is empty
    if(!buffer.at(buffer_pos).is_eof_block()) {
      cursor += buffer.at(buffer_pos).read(dest + cursor, len - cursor);
    }
    if(buffer_pos < comp_buffer_count) {
      if(buffer.at(buffer_pos).is_at_end()) {
        buffer.at(buffer_pos).clear_buffer(); // destroy current buffer
        if(IS_EOF == 1 && buffer_pos == comp_buffer_count - 1) {
          IS_EOB = 1; // Rcout << "EOB reached\n";
          return(cursor);
        }
        buffer_pos++; // increment
      }
    } // reading will always end with end of buffer being cleared
  }
  return(cursor);
}

unsigned int BAMReader_Multi::ignore(unsigned int len) {  
  // Read from current buffer
  unsigned int cursor = 0;
  if(IS_EOB == 1) return(cursor);
  while(cursor < len) {
    if(buffer_pos == comp_buffer_count && IS_EOF != 1) {
      read_from_file(n_bgzf * n_threads);
      decompress(n_bgzf * n_threads);
    } // reading will always start with reading buffer if current is empty
    if(!buffer.at(buffer_pos).is_eof_block()) {
      cursor += buffer.at(buffer_pos).ignore(len - cursor);
    }
    if(buffer_pos < comp_buffer_count) {
      if(buffer.at(buffer_pos).is_at_end()) {
        buffer.at(buffer_pos).clear_buffer(); // destroy current buffer
        if(IS_EOF == 1 && buffer_pos == comp_buffer_count - 1) {
          IS_EOB = 1; // Rcout << "EOB reached\n";
          return(cursor);
        }
        buffer_pos++; // increment
      }
    } // reading will always end with end of buffer being cleared
  }
  return(cursor);
}

bool BAMReader_Multi::eof() {
  if(IS_EOF == 1) {
      return (true);
  } else {
    if(IN->eof()) {
      IS_EOF = 1; // Rcout << "EOF reached\n";
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
    if(IS_EOF == 1 && buffer_pos == comp_buffer_count - 1 && buffer.at(buffer_pos).is_at_end()) {
      buffer.at(buffer_pos).clear_buffer();
      IS_EOB = 1; // Rcout << "EOB reached\n";
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