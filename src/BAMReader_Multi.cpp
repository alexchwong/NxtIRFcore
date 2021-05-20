#include <stddef.h>
#include "BAMReader_Multi.h"
#include <stdexcept>

const char BAMReader_Multi::bamGzipHead[16+1] = 
		"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00";

buffer_chunk::buffer_chunk() {
  bgzf_pos = 0;
  max_buffer = 0;
  max_decompressed = 0;
  pos = 0;
  end_pos = 65536;
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
  end_pos = 65536;
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
        Rcout << "Exception during BAM decompression - inflateInit2() fail: (" << ret << ") "
          << ", bgzf pos = " << GetBGZFPos() << '\n';
        return(ret);
    }
    ret = inflate(&zs, Z_FINISH);
    if(ret != Z_OK && ret != Z_STREAM_END) {
        std::ostringstream oss;
        Rcout << "Exception during BAM decompression - inflate() fail: (" << ret << ") "
          << ", bgzf pos = " << GetBGZFPos() << '\n';
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
        Rcout << "CRC fail during BAM decompression" << ", bgzf pos = " << GetBGZFPos() << '\n';
        return(ret);
    }
    // pos = 0;
    if(end_pos < max_decompressed) max_decompressed = end_pos;
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
  if(len == 0) return(0);

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
  BAM_BLOCK_CURSOR = 0;
  IN = NULL;
}

BAMReader_Multi::BAMReader_Multi(uint64_t block_begin, unsigned int begin_offset,
      uint64_t block_end, unsigned int end_offset) {
  IS_EOF = 0;
  IS_EOB = 0;   // End of file AND buffer
  IS_FAIL = 0;
  IS_LENGTH = 0;
  BAM_READS_BEGIN = 0;
  BAM_BLOCK_CURSOR = 0;
  
  begin_block_offset = block_begin;
  begin_read_offset = begin_offset;
  end_block_offset = block_end;
  end_read_offset = end_offset;
  
  BAM_BLOCK_CURSOR = block_begin;
  IN = NULL;
}

void BAMReader_Multi::AssignTask(std::istream *in_stream, 
      uint64_t block_begin, unsigned int begin_offset,
      uint64_t block_end, unsigned int end_offset) {
  IS_EOF = 0;
  IS_EOB = 0;   // End of file AND buffer
  IS_FAIL = 0;
  IS_LENGTH = 0;
  BAM_READS_BEGIN = 0;
  BAM_BLOCK_CURSOR = 0;
  
  begin_block_offset = block_begin; 
  begin_read_offset = begin_offset;
  end_block_offset = block_end;
  end_read_offset = end_offset;
  
  BAM_BLOCK_CURSOR = block_begin; 
  // Rcout << "block begin: " << BAM_BLOCK_CURSOR << '\n';
  // Rcout << begin_block_offset << " " << begin_read_offset
    // << ", " << end_block_offset << " " << end_read_offset << '\n';
  
  IN = in_stream;
}

// Destructor
BAMReader_Multi::~BAMReader_Multi() {
  if(buffer.size() > 0) {
    for(unsigned int i = 0; i < buffer.size(); i++) {
      buffer.at(i).clear_buffer();
    }
  }
}

void BAMReader_Multi::SetInputHandle(std::istream *in_stream) {
	IN = in_stream;
  //  get length of file:
  if(in_stream != &std::cin) {
    IN->seekg (0, std::ios_base::end);
    IS_LENGTH = IN->tellg();
    IN->seekg (BAM_BLOCK_CURSOR, std::ios_base::beg);    
  }
}

// OK.
unsigned int BAMReader_Multi::readBamHeader(
    std::vector<uint64_t> &block_begins, 
    std::vector<unsigned int> &read_offsets,
    unsigned int n_workers) {
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
  if(buffer_pos == comp_buffer_count) {
    BAM_READS_BEGIN = tellg();
    BAM_BLOCK_CURSOR = BAM_READS_BEGIN;
  }
  
  return(ProfileBAM(block_begins, read_offsets, n_workers));
}

void BAMReader_Multi::fillChrs(std::vector<chr_entry> &chrs_dest) {
  for(unsigned int i = 0; i < chrs.size(); i++) {
    chrs_dest.push_back(chrs.at(i));
  }
}

int BAMReader_Multi::getBGZFstarts(std::vector<uint64_t> & BGZF_begins) {
  BGZF_begins.clear();
  
  IN->clear();
  IN->seekg (BAM_READS_BEGIN, std::ios_base::beg);
  
  unsigned int bgzf_size = 0;
  
  Progress p(IS_LENGTH, true);
  size_t cursor = BAM_READS_BEGIN;
  p.increment(cursor);
  while(!IN->eof() && bgzf_size != 10) {
    BGZF_begins.push_back(IN->tellg());
    
    stream_uint16 u16;
    char GzipCheck[16];
    
    IN->read(GzipCheck, 16);
    if(strncmp(bamGzipHead, GzipCheck, 16) != 0) {
      Rcout << "This does not seem to be a legit BAM file\n";
      IN->clear();
      IN->seekg (BAM_READS_BEGIN, std::ios_base::beg);
      return(-1);
    }
    IN->read(u16.c, 2);
    bgzf_size = u16.u + 1 - 2  - 16;
    
    IN->ignore(bgzf_size);
    p.increment((size_t)IN->tellg() - cursor);
    cursor = IN->tellg();
  }
  IN->clear();
  IN->seekg (BAM_READS_BEGIN, std::ios_base::beg);

  if(bgzf_size == 10) {
    return(0);
  } else {
    return(-1);
  }
}

unsigned int BAMReader_Multi::ProfileBAM(
    std::vector<uint64_t> &block_begins, 
    std::vector<unsigned int> &read_offsets, 
    unsigned int target_n_threads) {
      
  if(BAM_READS_BEGIN == 0) return(0);
  IN->seekg (BAM_READS_BEGIN, std::ios_base::beg);
  
  std::vector<uint64_t> temp_begins;
  std::vector<unsigned int> temp_last_read_offsets;
  stream_uint32 u32;
  unsigned int bytes_read;

  // scan file to obtain a list of bgzf offsets
  int ret = getBGZFstarts(temp_begins);
  if(ret != 0) {
    return(0);
  }
  // assign n blocks to check if they are self-contained bgzf (i.e. they start and end at read boundary)
  unsigned int divisor = (temp_begins.size()/ target_n_threads);
  unsigned int i = 0;
  Rcout << "Checking whether bgzf blocks contain whole number of reads\n";
  while(i < temp_begins.size() && block_begins.size() < target_n_threads) {
    block_begins.push_back(temp_begins.at(i));
    read_offsets.push_back(0);
    i+=divisor;
  }
  bool is_self_contained = true;
  for(unsigned int j = 0; j < block_begins.size(); j++) {
    buffer_chunk * temp_buffer = new buffer_chunk;
    
    IN->seekg (block_begins.at(j), std::ios_base::beg);
    
    temp_buffer->read_from_file(IN);
    temp_buffer->decompress();
    
    while(temp_buffer->GetPos() < temp_buffer->GetMaxBufferDecompressed()) {
      bytes_read = temp_buffer->read(u32.c, 4);
      if(bytes_read < 4) {
        is_self_contained = false;
        delete temp_buffer; break;
      }
      bytes_read = temp_buffer->ignore(u32.u);
      if(bytes_read < u32.u) {
        is_self_contained = false;
        delete temp_buffer; break;
      }
    }
    if(!is_self_contained) {
      delete temp_buffer; break;
    }
    delete temp_buffer;
  }
  if(is_self_contained) {
  // if so, exit with these n blocks (where n == target_n_threads)
    // push EOF
    block_begins.push_back(temp_begins.at(temp_begins.size() - 1));
    read_offsets.push_back(0);
    return(temp_begins.size());
  } else {
    Rcout << "BAM reads appear to be split across BGZF blocks, requiring full indexing...\n";
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    temp_begins.clear(); block_begins.clear(); read_offsets.clear();    
  }
 
  uint64_t new_begin = 0;
  unsigned int last_read_offset = 0;
  unsigned int new_offset = 0;
    
  bool end_of_buffer = false;
  bool break_at_read_head = false;
  bool break_at_read_body = false;
  unsigned int head_offset = 0;
  uint64_t block_begin;
    
  IN->seekg (BAM_READS_BEGIN, std::ios_base::beg);
  while(!eof()) {
    buffer_chunk * temp_buffer = new buffer_chunk;
    block_begin = tellg();
    temp_buffer->read_from_file(IN);
    temp_buffer->decompress();

    if(!temp_buffer->is_eof_block()) {
      // Deal with previous block and output block start:
      new_begin = block_begin; 
      if(break_at_read_body) {
        bytes_read = temp_buffer->ignore(new_offset); // assume this read is completed
      } else if(break_at_read_head) {
        bytes_read = temp_buffer->read(u32.c + head_offset, 4 - head_offset);
        new_offset = u32.u - head_offset + 4;
        bytes_read = temp_buffer->ignore(u32.u); // assume this read is completed
      } else {
        new_offset = 0;
      }
      temp_begins.push_back(new_begin);
      // first_read_offsets.push_back(new_offset);
      Rcout << new_begin << '\t' << new_offset << '\n';
      
      head_offset = 0;
      end_of_buffer = false;
      break_at_read_head = false;
      break_at_read_body = false;
      while(!end_of_buffer) {
        last_read_offset = temp_buffer->GetPos();
        bytes_read = temp_buffer->read(u32.c, 4);
        // Rcout << bytes_read << '\t';
        if(bytes_read < 4) {
          end_of_buffer = true;
          break_at_read_head = true;
          head_offset = bytes_read;
          break;
        }
        bytes_read = temp_buffer->ignore(u32.u);
        // Rcout << bytes_read << '\n';
        if(bytes_read < u32.u) {
          new_offset = u32.u - bytes_read;
          end_of_buffer = true;
          break_at_read_body = true;
          break;
        }
        if(temp_buffer->is_at_end()) {
          new_offset = 0;
          end_of_buffer = true;
        }
      }
    } else {
      IS_EOF = 1; break;
    }
    delete temp_buffer;
    
    if(!break_at_read_head && !break_at_read_body) last_read_offset = 0;
    temp_last_read_offsets.push_back(last_read_offset);
  }
  // Rcout << "EOF block recorded: " << block_begin << "\n";
    
  // reset
  IN->seekg (BAM_READS_BEGIN, std::ios_base::beg);
  IS_EOF = 0;
  
  // divide cake into n_threads:
  divisor = (temp_begins.size()/ target_n_threads);
  i = 0;
  while(i < temp_begins.size() && block_begins.size() < target_n_threads) {
    block_begins.push_back(temp_begins.at(i));
    read_offsets.push_back(temp_last_read_offsets.at(i));
    i+=divisor;
  }
  // Return position of EOF block:
  block_begins.push_back(block_begin);
  read_offsets.push_back(0);

  return(temp_begins.size());
}

int BAMReader_Multi::read_from_file(unsigned int n_blocks) {
  unsigned int i = 0;
  if(IS_EOF == 1) return(i);
  
  IN->clear();    // In case another thread has hit the EOF bit
  IN->seekg (BAM_BLOCK_CURSOR, std::ios_base::beg);
  
  // Rcout << "Reading blocks starting from " << comp_buffer_count << " at file position "
    // << BAM_BLOCK_CURSOR << '\n';
  
  buffer.resize(comp_buffer_count + n_blocks);
  while(i < n_blocks) {
    if(IS_EOF == 1) return(i);
    if(BAM_BLOCK_CURSOR == end_block_offset && end_read_offset == 0) {
      IS_EOF = 1;
      return(i);
    }
    buffer.at(comp_buffer_count).SetBGZFPos(BAM_BLOCK_CURSOR);
    int ret = buffer.at(comp_buffer_count).read_from_file(IN);
    if(ret != 0) {
      if(ret == 1) {
        IS_EOF = 1;
      } else {
        Rcout << "Error reading file, error code: " << ret << ", bgzf pos = " 
          << buffer.at(comp_buffer_count).GetBGZFPos() << '\n';
      }
      return(i);
    } else {
      i++;
    }
    // Set begin cursor if BAM_BLOCK_CURSOR == begin_block_offset
    if(BAM_BLOCK_CURSOR == begin_block_offset) {
      buffer.at(comp_buffer_count).SetPos(begin_read_offset);
    }
    if(BAM_BLOCK_CURSOR == end_block_offset && end_block_offset > 0) {
      buffer.at(comp_buffer_count).SetEndPos(end_read_offset);
      IS_EOF = 1;   // Set virtual EOF
      i--;
    }

    comp_buffer_count++; 
    if(buffer.at(comp_buffer_count - 1).GetMaxBuffer() == 10) {
      IS_EOF = 1;
      i--;
    }
    BAM_BLOCK_CURSOR = IN->tellg();
  }
  BAM_BLOCK_CURSOR = IN->tellg();
  return(i);
}

int BAMReader_Multi::decompress(unsigned int n_blocks) {
  if(IS_EOB == 1) return(0);
  unsigned int end_blocks = min(comp_buffer_count, n_blocks + buffer_count);

  // Rcout << "Decompressing blocks starting from " << buffer_count << '\n';
  
  for(unsigned int i = buffer_count; i < end_blocks; i++) {
    if(!buffer.at(i).is_decompressed()) buffer.at(i).decompress();
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
      if(auto_load_data) {
        read_from_file(n_bgzf);
        decompress(n_bgzf);        
      } else {
        return(cursor);
      }
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
      if(auto_load_data) {
        read_from_file(n_bgzf);
        decompress(n_bgzf);        
      } else {
        return(cursor);
      }
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
