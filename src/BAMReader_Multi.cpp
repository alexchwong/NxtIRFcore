#include <stddef.h>
#include "BAMReader_Multi.h"
#include <stdexcept>

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
  // if(buffer) free(buffer);
  // if(decompressed_buffer) free(decompressed_buffer);
  buffer = NULL;
  decompressed_buffer = NULL;
  decompressed = false;
  bgzf_pos = 0;
  max_buffer = 0;
  max_decompressed = 0;
  pos = 0;
  end_pos = 65536;
}

int buffer_chunk::clear_buffer() {
  if(buffer) free(buffer);
  if(decompressed_buffer) free(decompressed_buffer);
  buffer = NULL;
  decompressed_buffer = NULL;
  decompressed = false;
  bgzf_pos = 0;
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

unsigned int buffer_chunk::peek(char * dest, unsigned int len) {  
  if(!is_decompressed()) return(0);
  if(pos >= max_decompressed) return(0);
  
  unsigned int bytes_to_read = min(max_decompressed, pos + len) - pos;
  memcpy(dest, decompressed_buffer + pos, bytes_to_read);
  // pos += bytes_to_read;
  return(bytes_to_read);
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
    bool verbose,
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
  
  return(ProfileBAM(block_begins, read_offsets, verbose, n_workers));
}

void BAMReader_Multi::fillChrs(std::vector<chr_entry> &chrs_dest) {
  for(unsigned int i = 0; i < chrs.size(); i++) {
    chrs_dest.push_back(chrs.at(i));
  }
}

int BAMReader_Multi::getBGZFstarts(std::vector<uint64_t> & BGZF_begins, bool verbose) {
  BGZF_begins.clear();
  
  IN->clear();
  IN->seekg (BAM_READS_BEGIN, std::ios_base::beg);
  
  // Min 100Mb file to bother showing any progress
  bool verbose_2 = verbose && (IS_LENGTH > 100000000);
  
  unsigned int bgzf_size = 0;
  unsigned int bgzf_check_threshold = 10000;    // Only check Gzip block every 10k runs
  Progress p(IS_LENGTH, verbose_2);
  p.increment(BAM_READS_BEGIN);
  uint64_t last_reported_pos = BAM_READS_BEGIN;
  while(!IN->eof() && bgzf_size != 10) {
    BGZF_begins.push_back(IN->tellg());
    // Rcout << "BGZF pos " << IN->tellg() << '\t';
    stream_uint16 u16;
    
    if(BGZF_begins.size() % bgzf_check_threshold == 1) {
      char * GzipCheck = (char*)malloc(16);
      IN->read(GzipCheck, 16);
      if(strncmp(bamGzipHead, GzipCheck, 16) != 0) {
        Rcout << "This does not seem to be a legit BAM file\n";
        IN->clear();
        IN->seekg (BAM_READS_BEGIN, std::ios_base::beg);
        return(-1);
      }
      free(GzipCheck);
      p.increment((uint64_t)IN->tellg() - last_reported_pos);
      last_reported_pos = IN->tellg();
    } else {
      // IN->ignore(16);
      IN->seekg (16, std::ios_base::cur);
    }

    IN->read(u16.c, 2);
    bgzf_size = u16.u + 1 - 2  - 16;
    // Rcout << " bgzf_size " << bgzf_size << '\n';
    // IN->ignore(bgzf_size);
    IN->seekg (bgzf_size, std::ios_base::cur);
  }
  p.increment(IS_LENGTH - last_reported_pos);
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
    bool verbose,
    unsigned int target_n_threads) {
      
  if(BAM_READS_BEGIN == 0) return(0);
  IN->seekg (BAM_READS_BEGIN, std::ios_base::beg);
  
  std::vector<uint64_t> temp_begins;
  std::vector<unsigned int> temp_last_read_offsets;
  stream_uint32 u32;
  // unsigned int bytes_read;

  // scan file to obtain a list of bgzf offsets
  int ret = getBGZFstarts(temp_begins, verbose);
  if(ret != 0) {
    return(0);
  }
  // assign n blocks to check if they are self-contained bgzf (i.e. they start and end at read boundary)
  unsigned int divisor = (temp_begins.size()/ target_n_threads);
  unsigned int i = 0;
  // Rcout << "temp_begins.size() == " << temp_begins.size() <<
    // " target_n_threads " << target_n_threads <<
    // " divisor " << divisor << '\n';
  // if(verbose) Rcout << "Checking whether bgzf blocks contain whole number of reads\n";
  while(i < temp_begins.size() && block_begins.size() < target_n_threads) {
    block_begins.push_back(temp_begins.at(i));
    // Rcout << "block " << temp_begins.at(i) << '\t';
    read_offsets.push_back(0);
    i+=divisor;
    // Rcout << '\n';
  }
  bool is_self_contained = true;
  
  if(target_n_threads == 1) {
    // No need to profile
    block_begins.push_back(temp_begins.at(temp_begins.size() - 1));
    read_offsets.push_back(0);
    return(temp_begins.size());    
  }
  
  for(unsigned int j = 0; j < block_begins.size(); j++) {
    buffer_chunk * temp_buffer = new buffer_chunk;
    
    IN->seekg (block_begins.at(j), std::ios_base::beg);
    // Rcout << "BGZF begin " << block_begins.at(j) << '\t';
    temp_buffer->read_from_file(IN);
    temp_buffer->decompress();
    
    if(!temp_buffer->is_eof_block()) {
      while(temp_buffer->GetMaxBufferDecompressed() - temp_buffer->GetPos() > 36) {
        temp_buffer->read(u32.c, 4);
        if(u32.u <= temp_buffer->GetMaxBufferDecompressed() - temp_buffer->GetPos()) {
          temp_buffer->ignore(u32.u);
        } else {
          break;
        }
      }
      if(!temp_buffer->is_at_end()) {
        is_self_contained = false;
      }
    }

    delete temp_buffer;
    if(!is_self_contained) break;
  }
  if(is_self_contained) {
    if(verbose) Rcout << "BAM is self contained\n";
  // if so, exit with these n blocks (where n == target_n_threads)
    // push EOF
    block_begins.push_back(temp_begins.at(temp_begins.size() - 1));
    read_offsets.push_back(0);
    return(temp_begins.size());
  } else {
    if(verbose) Rcout << "BAM reads appear to be split across BGZF blocks, requiring full indexing...\n";

    block_begins.resize(0); read_offsets.resize(0);    
  }

  // setwd("d:/Alex/Vignette/Cpp Optim/")
  // IRFinder("Aligned.out.bam", "test", file.path(tempdir(), "Reference"), verbose = TRUE)
  
  if(IN->eof()) IN->clear();
  SetAutoLoad(false);
  IN->seekg (BAM_READS_BEGIN, std::ios_base::beg);
  BAM_BLOCK_CURSOR = BAM_READS_BEGIN;
  
  Progress p(temp_begins.size(), verbose);
  
  while(read_from_file(100) > 0) {
    p.increment(comp_buffer_count - buffer_count);
    decompress(true);
    
    if(buffer.at(buffer_pos).GetPos() > 0) {
      // Rcout << "Trying to go to the beginning of next buffer\n";
      if(!GotoNextRead(false)) break; // goes to the first full read of the next line     
    }
    // Rcout << "BGZF# " << buffer_pos << '\t';
    // Rcout << "Block read offset " << buffer.at(buffer_pos).GetPos() << '\n';
    temp_last_read_offsets.push_back(buffer.at(buffer_pos).GetPos()); // Gets current position of current buffer
    
    while(isReadable()) {
      while(1) {
        // strict step onto next read
        if(!GotoNextRead(true)) break;
      }

      // If ends at full read, ignore() automatically goes to pos zero of next buffer
      if(buffer.at(buffer_pos).GetPos() == 0) {   
        if(isReadable()) {
          temp_last_read_offsets.push_back(0); 
          // Rcout << "BGZF# " << buffer_pos << '\t';
          // Rcout << "Block read offset 0";   // First read is at zero of new bgzf
        } else if(buffer.at(buffer_pos).is_eof_block()) {
          temp_last_read_offsets.push_back(0);
          break;
        }
      } else {
        // Goto next read non-strictly (i.e. read across bgzf boundaries)
        // if this is not possible, then read some more buffers
        if(buffer_pos == comp_buffer_count - 1) break;
        if(!GotoNextRead(false)) break; 
        temp_last_read_offsets.push_back(buffer.at(buffer_pos).GetPos()); // record next position
        // Rcout << "BGZF# " << buffer_pos << '\t';
        // Rcout << "Block read offset " << buffer.at(buffer_pos).GetPos() << '\n';
      }
    }
  }
  if(verbose) Rcout << "Extended profiling finished\n";
  if(temp_last_read_offsets.size() == temp_begins.size() - 1) {
    // if(verbose) Rcout << "Pushing EOF block position\n";
    temp_last_read_offsets.push_back(0); 
  } else if(temp_last_read_offsets.size() != temp_begins.size()) {
    Rcout << "BGZF block counts mismatch between BGZF positions and first read positions\n";
    Rcout << "temp_begins.size() " << temp_begins.size() << '\t'
      << "temp_last_read_offsets.size() " << temp_last_read_offsets.size() << '\n';
    return(0);
  }
  
  IN->clear();
  IN->seekg (BAM_READS_BEGIN, std::ios_base::beg);
  IS_EOF = 0;
  
  // divide cake into n_threads:
  divisor = 1 + (temp_begins.size()/ target_n_threads);
  i = 0;
  while(i < temp_begins.size() && block_begins.size() < target_n_threads) {
    block_begins.push_back(temp_begins.at(i));
    read_offsets.push_back(temp_last_read_offsets.at(i));
    i+=divisor;
  }
  // Return position of EOF block:
  block_begins.push_back(temp_begins.at(temp_begins.size() - 1));
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
      if(i == 0) IS_EOB = 1;
      return(i);
    }
    buffer.at(comp_buffer_count).SetBGZFPos(BAM_BLOCK_CURSOR);
    int ret = buffer.at(comp_buffer_count).read_from_file(IN);
    if(ret != 0) {
      if(ret == 1) {
        IS_EOF = 1;
        if(i == 0) IS_EOB = 1;
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

int BAMReader_Multi::decompress(bool allow_openmp) {
  if(IS_EOB == 1) return(0);
  if(allow_openmp) {
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(unsigned int i = buffer_count; i < comp_buffer_count; i++) {
      if(!buffer.at(i).is_decompressed()) {
        buffer.at(i).decompress();
      }
    }
    
    for(unsigned int i = buffer_count; i < comp_buffer_count; i++) {
      if(buffer.at(i).is_decompressed()) buffer_count+=1;
    }
  } else {
    for(unsigned int i = buffer_count; i < comp_buffer_count; i++) {
      if(!buffer.at(i).is_decompressed()) {
        buffer.at(i).decompress();
        buffer_count += 1;
      }
    }
  }
  
  if(buffer_count != comp_buffer_count) {
    Rcout << "Some buffers were not decompressed";
  }

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
        decompress();
      } else {
        return(cursor);
      }
    } // reading will always start with reading buffer if current is empty
    if(!buffer.at(buffer_pos).is_eof_block()) {
      cursor += buffer.at(buffer_pos).read(dest + cursor, len - cursor);
    }
    if(buffer.at(buffer_pos).is_at_end()) {
      buffer.at(buffer_pos).clear_buffer(); // destroy current buffer
      buffer_pos++; // increment
      if(buffer_pos == comp_buffer_count) return(cursor);
      if(IS_EOF == 1 && buffer_pos == comp_buffer_count - 1 && 
          buffer.at(buffer_pos).GetRemainingBytes() == 0) {
        IS_EOB = 1; // Rcout << "EOB reached\n";
        return(cursor);
      }
      // if(buffer_pos < comp_buffer_count && !buffer.at(buffer_pos).is_decompressed()) {
        // decompress();
      // }
    } // reading will always end with end of buffer being cleared

  }
  return(cursor);
}

unsigned int BAMReader_Multi::peek(char * dest, unsigned int len) {  
  // Read from current buffer
  unsigned int cursor = 0;
  unsigned int temp_buffer_pos = buffer_pos;
  if(IS_EOB == 1) return(cursor);
  while(cursor < len) {
    if(temp_buffer_pos == comp_buffer_count && IS_EOF != 1) {
      if(auto_load_data) {
        read_from_file(n_bgzf);
        decompress();
      } else {
        return(cursor);
      }
    } // reading will always start with reading buffer if current is empty
    if(!buffer.at(temp_buffer_pos).is_eof_block()) {
      cursor += buffer.at(temp_buffer_pos).peek(dest + cursor, len - cursor);
    }
    if(buffer.at(temp_buffer_pos).is_at_end()) {
      // buffer.at(buffer_pos).clear_buffer(); // destroy current buffer
      temp_buffer_pos++; // increment
      if(temp_buffer_pos == comp_buffer_count) return(cursor);
      if(IS_EOF == 1 && temp_buffer_pos == comp_buffer_count - 1 &&
          buffer.at(temp_buffer_pos).GetRemainingBytes() == 0) {
        // IS_EOB = 1; // Rcout << "EOB reached\n";
        return(cursor);
      }
      // if(buffer_pos < comp_buffer_count && !buffer.at(buffer_pos).is_decompressed()) {
        // decompress();
      // }
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
        decompress();        
      } else {
        return(cursor);
      }
    } // reading will always start with reading buffer if current is empty
    if(!buffer.at(buffer_pos).is_eof_block()) {
      cursor += buffer.at(buffer_pos).ignore(len - cursor);
    }
    if(buffer.at(buffer_pos).is_at_end()) {
      buffer.at(buffer_pos).clear_buffer(); // destroy current buffer
      buffer_pos++; // increment
      if(buffer_pos == comp_buffer_count) return(cursor);
      if(IS_EOF == 1 && buffer_pos == comp_buffer_count - 1 && 
          buffer.at(buffer_pos).GetRemainingBytes() == 0) {
        IS_EOB = 1; // Rcout << "EOB reached\n";
        return(cursor);
      }
      // if(buffer_pos < comp_buffer_count && !buffer.at(buffer_pos).is_decompressed()) {
        // decompress();
      // }
    } // reading will always end with end of buffer being cleared
  }
  return(cursor);
}

bool BAMReader_Multi::isReadable() {
  // If not last decompressed buffer:
  if(buffer_pos < buffer_count - 2) return(true);
  if(buffer_pos >= buffer_count) return(false);
  // If next buffer is not EOF buffer:
  // Else, buffer is last available
  
  stream_uint32 u32;
  
  if(buffer_pos == buffer_count - 2) {
    if(buffer.at(buffer_pos).GetRemainingBytes() +
      buffer.at(buffer_pos + 1).GetRemainingBytes() < 4) return(false);
    
    peek(u32.c, 4);
    
    if(buffer.at(buffer_pos).GetRemainingBytes() +
      buffer.at(buffer_pos + 1).GetRemainingBytes() < 4 + u32.u) return(false);
  } else {
    if(buffer.at(buffer_pos).GetRemainingBytes() < 4) return(false);
    
    peek(u32.c, 4);
    
    if(buffer.at(buffer_pos).GetRemainingBytes() < 4 + u32.u) return(false);
  }

  return(true);
}

bool BAMReader_Multi::isReadableStrict() {
  // Strict means is the next read readable without going to the next buffer?
  if(buffer_pos >= buffer_count) return(false);
  // If next buffer is not EOF buffer:
  // Else, buffer is last availabl
  if(buffer.at(buffer_pos).GetRemainingBytes() < 4) return(false);
  
  stream_uint32 u32;
  buffer.at(buffer_pos).peek(u32.c, 4);
  if(buffer.at(buffer_pos).GetRemainingBytes() < 4 + u32.u) return(false);
  
  return(true);
}

bool BAMReader_Multi::GotoNextRead(bool strict) {
  if(strict) {
    if(buffer_pos >= buffer_count) return(false);
    if(buffer.at(buffer_pos).GetRemainingBytes() < 4) return(false);
    
    stream_uint32 u32;
    buffer.at(buffer_pos).peek(u32.c, 4);
    if(buffer.at(buffer_pos).GetRemainingBytes() <= 4 + u32.u) return(false);
    
    ignore(4 + u32.u);
    return(true);
  } else {
    if(buffer_pos >= buffer_count) return(false);
    if(buffer_pos == buffer_count - 2 && buffer.at(buffer_count - 1).is_eof_block()) {
      // cautious
      if(buffer.at(buffer_pos).GetRemainingBytes() < 4) return(false);

      stream_uint32 u32;
      buffer.at(buffer_pos).peek(u32.c, 4);
      if(buffer.at(buffer_pos).GetRemainingBytes() <= 4 + u32.u) return(false);
      ignore(4 + u32.u);
      return(true);
    } else if(buffer_pos < buffer_count - 1) {
      // Doable:
      stream_uint32 u32;
      if(read(u32.c, 4) < 4) Rcout << "Unable to read 4 bytes when supposed to\n";
      if(ignore(u32.u) < u32.u) Rcout << "Unable to read rest of read when supposed to\n";
      return(true);
    } else {
      if(buffer.at(buffer_pos).GetRemainingBytes() < 4) return(false);

      stream_uint32 u32;
      buffer.at(buffer_pos).peek(u32.c, 4);
      if(buffer.at(buffer_pos).GetRemainingBytes() <= 4 + u32.u) return(false);

      ignore(4 + u32.u);
      return(true);
    }
  }
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
