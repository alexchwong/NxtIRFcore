#include "covFile.h"

const char covFile::bamEOF[covFile::bamEOFlength+1] =
		"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00";
const char covFile::bamGzipHead[covFile::bamGzipHeadLength+1] = 
		"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00";
const char covBuffer::bamGzipHead[covBuffer::bamGzipHeadLength+1] = 
		"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00";

// Constructor
covFile::covFile() {
    bufferPos = 0;
    bufferMax = 0;
    index_begin = 0;
    body_begin = 0;
    
    mode = "";
    chr_index = NULL;
}

covBuffer::covBuffer() {
    bufferPos = 0;
    bufferMax = 65536 - 18 - 8;
    l_file_buffer = 0;
    file_bufferPos = 0;
    
    file_buffer = NULL;
}

// Destructor
covFile::~covFile() {
  if(chr_index != NULL) {
    free(chr_index);
  }
}

covBuffer::~covBuffer() {
  if(file_buffer != NULL) {
    free(file_buffer);
  }
}

// covBuffer functions

bool covBuffer::BufferIsFull(unsigned int threshold) {
  // Returns true if writing another threshold bits will lead to overflow (default threshold = 8)
  return(bufferPos + threshold > bufferMax);
}

int covBuffer::write(char * src, unsigned int len) {
  
  unsigned int remaining_bytes = 0;
  unsigned int src_pos = 0;
  // Read next block if buffer empty or if pos is at end of buffer
  if(bufferPos == bufferMax) {
    WriteBuffer();        
  }
  
  if (len <= bufferMax - bufferPos) {
    memcpy(&buffer[bufferPos], &src[src_pos], len);
    bufferPos += len;
    return(Z_OK);
  } else {
    // If len will overflow buffer, then fill to bufferMax then flush
    remaining_bytes = len - (bufferMax - bufferPos);
    memcpy(&buffer[bufferPos], &src[src_pos], bufferMax - bufferPos);
    src_pos += bufferMax - bufferPos;
    bufferPos = bufferMax;
    WriteBuffer();
    
    while(remaining_bytes > bufferMax) {        
      remaining_bytes -= bufferMax;
      memcpy(&buffer[0], &src[src_pos], bufferMax);
      src_pos += bufferMax;
      bufferPos = bufferMax;
      WriteBuffer();
    }
    
    memcpy(&buffer[0], &src[src_pos], remaining_bytes);
    bufferPos += remaining_bytes;
    src_pos += remaining_bytes;
  }
  return(0);
}

int covBuffer::WriteBuffer() {
  
  stream_uint16 u16;
  stream_uint32 u32;
  uint32_t crc;
  z_stream zs;
  
  zs.zalloc = NULL; zs.zfree = NULL;
  zs.msg = NULL;
  zs.next_in  = (Bytef*)buffer;
  zs.avail_in = bufferPos;
  zs.next_out = (Bytef*)compressed_buffer;
  zs.avail_out = 65536 - 18 - 8;
  int ret = deflateInit2(&zs, 6, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY); // -15 to disable zlib header/footer
  
  ret = deflate(&zs, Z_FINISH);
  ret = deflateEnd(&zs);
  
  int block_len = zs.total_out + 18 + 8;
  
  // Write to file_buffer
  char *data_tmp;  
  file_buffer = (char *)realloc((data_tmp = file_buffer), l_file_buffer += block_len);
  
  memcpy(&file_buffer[file_bufferPos], bamGzipHead, bamGzipHeadLength);
  file_bufferPos += bamGzipHeadLength;
  
  u16.u = block_len - 1;
  memcpy(&file_buffer[file_bufferPos], u16.c, 2);
  file_bufferPos += 2;
  
  memcpy(&file_buffer[file_bufferPos], compressed_buffer, zs.total_out);
  file_bufferPos += zs.total_out;
  
  crc = crc32(crc32(0L, NULL, 0L), (Bytef*)buffer, bufferPos);
  u32.u = crc;
  
  memcpy(&file_buffer[file_bufferPos], u32.c, 4);
  file_bufferPos += 4;
  
  u32.u = bufferPos;
  memcpy(&file_buffer[file_bufferPos], u32.c, 4);
  file_bufferPos += 4;
  
  bufferPos = 0;
  bufferMax = 65536 - 18 - 8;
  return(ret);
}

void covFile::SetOutputHandle(std::ostream *out_stream) {
  OUT = out_stream;

  mode = "write";
  bufferMax = 65536 - 18 - 8;
  bufferPos = 0;
  out_cur_seqID = 0;
  
  chr_index = (char*)malloc(1200);
  chr_index_alloc = 1200;
  chr_index_pos = 0;
  chr_coord = 0;
}

void covFile::SetInputHandle(std::istream *in_stream) {
    IS_EOF = 0;
    IS_FAIL = 0;
    IS_LENGTH = 0;

    IN = in_stream;
  
  // Identify EOF
    IN->seekg (0, std::ios_base::end);
    IS_LENGTH = IN->tellg();
    
    // Check EOF bit
    IN->seekg (-bamEOFlength, std::ios_base::end);
    
    char check_eof_buffer[covFile::bamEOFlength+1];
    IN->read(check_eof_buffer, bamEOFlength);
         
    if(strncmp(check_eof_buffer, bamEOF, bamEOFlength) == 0) {
        EOF_POS = IS_LENGTH - bamEOFlength;
        mode = "read";

    } else {
        // Rcout << "EOF bit not detected\n";
        EOF_POS = 0;
				IS_EOF = 1;
				IS_FAIL = 1;				
        mode = "";
    }
    IN->seekg (0, std::ios_base::beg);
}

int covFile::ReadBuffer() {
    // read compressed buffer
    if((size_t)IN->tellg() >= EOF_POS) {
        IS_EOF = 1;
        return(Z_STREAM_END);
    } else if(fail()) {
        return(Z_STREAM_ERROR);
    }

    stream_uint16 u16;
		int ret = 0;
		
    char GzipCheck[bamGzipHeadLength];
    IN->read(GzipCheck, bamGzipHeadLength);

     if(strncmp(bamGzipHead, GzipCheck, bamGzipHeadLength) != 0) {
        std::ostringstream oss;
        oss << "Exception during BAM decompression - BGZF header corrupt: (at " << IN->tellg() << " bytes) ";
        // throw(std::runtime_error(oss.str()));
				return(Z_BUF_ERROR);
    }

    IN->read(u16.c, 2);
    IN->read(compressed_buffer, u16.u + 1 - 2  - bamGzipHeadLength);

    bufferMax = 65536;
    z_stream zs;
    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.msg = NULL;
    zs.next_in = (Bytef*)compressed_buffer;
    zs.avail_in = u16.u + 1 - 2  - bamGzipHeadLength;
    zs.next_out = (Bytef*)buffer;
    zs.avail_out = bufferMax;

    stream_uint32 u32;
    memcpy(u32.c, &compressed_buffer[u16.u + 1 - 2 - bamGzipHeadLength - 8],4);

    ret = inflateInit2(&zs, -15);
    if(ret != Z_OK) {
        std::ostringstream oss;
        oss << "Exception during BAM decompression - inflateInit2() fail: (" << ret << ") ";
        // throw(std::runtime_error(oss.str()));
				return(ret);
    }
    ret = inflate(&zs, Z_FINISH);
    if(ret != Z_OK && ret != Z_STREAM_END) {
        std::ostringstream oss;
        oss << "Exception during BAM decompression - inflate() fail: (" << ret << ") ";
        // throw(std::runtime_error(oss.str()));
				return(ret);
    }
    ret = inflateEnd(&zs);
    
    bufferMax -= zs.avail_out;
    
    // check CRC
    uint32_t crc = crc32(crc32(0L, NULL, 0L), (Bytef*)buffer, bufferMax);
    // CRC check:
    if(u32.u != crc) {
        std::ostringstream oss;
        oss << "CRC fail during BAM decompression: (at " << IN->tellg() << " bytes) ";
        // throw(std::runtime_error(oss.str()));
				return(ret);
    }
    bufferPos = 0;
    
    return(ret);
}

int covFile::read(char * dest, unsigned int len) {
    
    unsigned int remaining_bytes = 0;
    unsigned int dest_pos = 0;
		int ret = 0;
    // Read next block if buffer empty or if pos is at end of buffer
    if(bufferMax == 0 || bufferPos == bufferMax) {
        ret = ReadBuffer();
				if(ret != Z_OK) return(ret);
    }
    
    if (len <= bufferMax - bufferPos) {
        memcpy(&dest[0], &buffer[bufferPos], len);
        bufferPos += len;
        return(Z_OK);
    } else {
        memcpy(&dest[dest_pos], &buffer[bufferPos], bufferMax - bufferPos);
        remaining_bytes = len - (bufferMax - bufferPos);
        dest_pos += bufferMax - bufferPos;
        bufferMax = 0;
        bufferPos = 0;
        ret = ReadBuffer();
				if(ret != Z_OK) return(ret);

        while(remaining_bytes > bufferMax) {
            memcpy(&dest[dest_pos], &buffer[0], bufferMax);
            remaining_bytes -= bufferMax;
            dest_pos += bufferMax;
            bufferMax = 0;
            bufferPos = 0;
						ret = ReadBuffer();
						if(ret != Z_OK) return(ret);
        }
        
        memcpy(&dest[dest_pos], &buffer[bufferPos], remaining_bytes);
        bufferPos += remaining_bytes;
        dest_pos += remaining_bytes;
    }
    return(Z_OK);
}



int covFile::write(char * src, unsigned int len) {
    
    unsigned int remaining_bytes = 0;
    unsigned int src_pos = 0;
    // Read next block if buffer empty or if pos is at end of buffer
    if(bufferPos == bufferMax) {
        WriteBuffer();        
    }
    
    if (len < bufferMax - bufferPos) {
        memcpy(&buffer[bufferPos], &src[src_pos], len);
        bufferPos += len;
        return(Z_OK);
    } else {
        remaining_bytes = len - (bufferMax - bufferPos);
        memcpy(&buffer[bufferPos], &src[src_pos], bufferMax - bufferPos);
        src_pos += bufferMax - bufferPos;
        bufferPos = bufferMax;
        WriteBuffer();

      while(remaining_bytes > bufferMax) {                 
          remaining_bytes -= bufferMax;
          memcpy(&buffer[0], &src[src_pos], bufferMax);
          src_pos += bufferMax;
          bufferPos = bufferMax;
          WriteBuffer();
      }
        
        memcpy(&buffer[bufferPos], &src[src_pos], remaining_bytes);
        bufferPos += remaining_bytes;
        src_pos += remaining_bytes;
    }
    return(0);
}


int covFile::WriteBuffer() {
    
    stream_uint16 u16;
    stream_uint32 u32;
    uint32_t crc;
    z_stream zs;

    zs.zalloc = NULL; zs.zfree = NULL;
    zs.msg = NULL;
    zs.next_in  = (Bytef*)buffer;
    zs.avail_in = bufferPos;
    zs.next_out = (Bytef*)compressed_buffer;
    zs.avail_out = 65536 - 18 - 8;
    int ret = deflateInit2(&zs, 6, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY); // -15 to disable zlib header/footer

    ret = deflate(&zs, Z_FINISH);
    ret = deflateEnd(&zs);
    
    int block_len = zs.total_out + 18 + 8;
    
    OUT->write(bamGzipHead, bamGzipHeadLength);
    u16.u = block_len - 1;
    OUT->write(u16.c, 2);
    
    OUT->write(compressed_buffer, zs.total_out);
    
    crc = crc32(crc32(0L, NULL, 0L), (Bytef*)buffer, bufferPos);
    u32.u = crc;
    OUT->write(u32.c, 4);
    u32.u = bufferPos;
    OUT->write(u32.c, 4);
    
    bufferPos = 0;
    bufferMax = 65536 - 18 - 8;
    return(ret);
}

int covFile::ignore(unsigned int len) {
    // Essentially copy read() but without memcpy etc.
    unsigned int remaining_bytes = len;

    if(bufferMax == 0 || bufferPos == bufferMax) {
        ReadBuffer();        
    }
    
    if (len <= bufferMax - bufferPos) {
        bufferPos += len;
        return(Z_OK);
    } else {
        remaining_bytes = len - (bufferMax - bufferPos);

        bufferMax = 0;
        bufferPos = 0;
        ReadBuffer();

        while(remaining_bytes > bufferMax) {
            remaining_bytes -= bufferMax;
            bufferMax = 0;
            bufferPos = 0;
            ReadBuffer();
        }
        
        bufferPos += remaining_bytes;
        // Note bufferPos can equal bufferMax. IN->tellg() will return the position of the next bgzf block
    }
    return(0);
}

bool covFile::eof() {
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

bool covFile::fail() {
    if(IS_FAIL == 1) {
        return (true);
    } else {
        if(IN->fail()) {
            IS_FAIL = 1;
            return (true);
        } else {
            return (false);
        }
    }
}

int covFile::ReadHeader() {
    IN->seekg (0, std::ios_base::beg);    
    chr_names.clear();
    chr_lens.clear();
    bufferPos = 0;
    bufferMax = 0;    
    
    char cov_header[4];
		int ret = read(cov_header,4);
		if(ret != Z_OK) {
        Rcout << "File is not BGZF compressed; unlikely to be COV file\n";
				return(ret);
		}
    std::string s_cov_header = "COV\x01";
    if(strncmp(cov_header, s_cov_header.c_str(), 4) != 0) {
        Rcout << "COV file has incorrect header!\n";
				return(-1);
    }
    
    stream_uint32 n_ref;
    std::string chrName;

    read(n_ref.c, 4);
    for(unsigned int i = 0; i < n_ref.u; i++) {
        stream_uint32 l_name;
        read(l_name.c, 4);
        
        char * c_name = new char[l_name.u];
        read(c_name, l_name.u);
        chrName = string(c_name, l_name.u - 1);
        chr_names.push_back(chrName);

        stream_uint32 l_ref;
        read(l_ref.c, 4);
        chr_lens.push_back(l_ref.u);
        
        delete[] c_name;
    }
    
    index_begin = IN->tellg();      // should be the start point of bgzf block containing index
    bufferPos = 0;
    bufferMax = 0;    
    // Rcout << "index_begin: " << index_begin << '\n';
    
    // keep profiling to identify where body_begin is
    stream_uint32 u32;
    for(unsigned int j = 0; j < 3; j++) {
        for(unsigned int i = 0; i < chr_names.size(); i++) {
            read(u32.c, 4);
            ignore(u32.u);            
             // Rcout << "refID: " << i << ", strand: " << j << ", index bytes: " << u32.u << '\n';
        }
    }
    body_begin = IN->tellg();
    // Rcout << "body_begin: " << body_begin << '\n';
    
    return(n_ref.u);
}

int covFile::FetchPos(const std::string seqname, const uint32_t start, const int strand,
    uint64_t * file_offset, uint32_t * block_start) {

    // Inputs seqname, and desired start position of query
    // Alters file_offset to point to the compressed offset position of bgzf block to start reading
    // Alters block_start to display the actual coordinate start of the first entry of block
    // Returns:
        // 0 = success
        // -1 = fail
    if(strand < 0 || strand > 2) {
        return -1;
    }        
       
    if(index_begin == 0) {
        ReadHeader();
        if(index_begin == 0) {
            return -1;
        }
    }
    
    int ref_index;
    auto it_chr = std::find(chr_names.begin(), chr_names.end(), seqname);
    if(it_chr == chr_names.end()) {
        return -1;
    } else {
        ref_index = distance(chr_names.begin(), it_chr);
        ref_index += strand * chr_names.size();
    }

    int i = 0;
    stream_uint32 u32;
    IN->seekg(index_begin, std::ios_base::beg);    
    bufferPos = 0;
    bufferMax = 0;    

//    ignore(4);     // ignore index size
    while(i < ref_index) {
        read(u32.c, 4);
        ignore(u32.u);
        i += 1;
    }
    
    // Now we are at the start of the ref_index
    stream_uint32 chr_block_size;
    stream_uint32 cur_block_start;
    stream_uint64 cur_offset;
    uint32_t prev_block_start;
    uint64_t prev_offset = 0;

    uint32_t block_counter = 0;
    prev_block_start = 0;
    
    read(chr_block_size.c, 4);     // use this to know when to stop
    while(block_counter < chr_block_size.u) {
        read(cur_block_start.c, 4);
        read(cur_offset.c, 8);
        block_counter += 12;
        if(cur_block_start.u > start) { // first entry should always equal zero so should not be immediately called
            break;
        } else {
            prev_block_start = cur_block_start.u;  
            prev_offset = cur_offset.u;
        }
    }

    // Rcout << "file offset: " << prev_offset + body_begin << ", block_start: " << prev_block_start << "\n";
    *file_offset = prev_offset + body_begin;
    *block_start = prev_block_start;
    return 0;
}

int covFile::FetchRLE(const std::string seqname, const uint32_t start, const uint32_t end, const int strand,
  std::vector<int> * values, std::vector<unsigned int> * lengths) {
 
  stream_int32 i32;
  stream_uint32 u32;
  
//  std::vector<int> values;
//  std::vector<unsigned int> lengths;
  
  uint64_t file_offset = 0;
  uint32_t block_start = 0;

  
  auto it_chr = std::find(chr_names.begin(), chr_names.end(), seqname);
  if(it_chr == chr_names.end()) {
      return -1;
  } else {
      int ref_index = distance(chr_names.begin(), it_chr);
      if(end > chr_lens[ref_index]) {
          return -1;
      }
  }
  
  if(FetchPos(seqname, start, strand, &file_offset, &block_start) != 0) {
      return -1;
  }
  
  IN->seekg(file_offset, std::ios_base::beg);
  bufferPos = 0;
  bufferMax = 0;    
  
  uint32_t prev_start = block_start;
  
  do {
      read(i32.c, 4);
      read(u32.c, 4);
      prev_start += u32.u;
  } while(prev_start < start);
  if(prev_start > start) {
      if(prev_start >= end) {
          (*values).push_back(i32.i);
          (*lengths).push_back((unsigned int)(end - start));
      } else {
          (*values).push_back(i32.i);
          (*lengths).push_back((unsigned int)(prev_start - start));         
      }
  }
  // main loop
  while(prev_start < end) {
      read(i32.c, 4);
      read(u32.c, 4);
      if(prev_start + u32.u >= end) {
          (*values).push_back(i32.i);
          (*lengths).push_back((unsigned int)(end - prev_start));
          break;
      } else {
          (*values).push_back(i32.i);
          (*lengths).push_back(u32.u);       
      }
      prev_start += u32.u;
  };
  return(0);
}

// Write operations

int covFile::FlushBody() {
  // Final operation after everything has been written
  
  stream_uint32 u32;
  // Flush final chromosome
  u32.u = chr_index_pos;
  write(u32.c, 4);
  // Now write the rest of the index buffer
  if(chr_index_pos > 0) {
    write(chr_index, chr_index_pos);
  }
  WriteBuffer();
  
  body.WriteBuffer();

  OUT->write(body.file_buffer, body.file_bufferPos);
  OUT->write(covFile::bamEOF, covFile::bamEOFlength);
  
  OUT->flush();
  return 0;
}

int covFile::WriteHeader(std::vector<std::string> s_chr, std::vector<int32_t> u_lens) {

  std::map< std::string, unsigned int > chrmap;
  
  for(unsigned int i = 0; i < s_chr.size(); i++) {
      chrmap.insert({s_chr[i], (unsigned int)u_lens[i]});
  }

  char zero = '\0';
  char wh_buffer[1000];
  
  std::string header = "COV\x01";
  strncpy(wh_buffer, header.c_str(), 4);
  write(wh_buffer,4);
  
  stream_uint32 u32;
  stream_uint64 u64;

  u32.u = s_chr.size();
  write(u32.c ,4);

  for (auto chr = chrmap.begin(); chr != chrmap.end(); chr++) {
    u32.u = chr->first.length() + 1;
    write(u32.c ,4);
    strncpy(wh_buffer, chr->first.c_str(), chr->first.length());
    write(wh_buffer, chr->first.length());
    write(&zero, 1);
    
    u32.u = chr->second;
    write(u32.c ,4);
    
    chr_names.push_back(chr->first);
    chr_lens.push_back(chr->second);
  }
  
  WriteBuffer();       // flush output
  
  // Initialize first chromosome entry here
  u32.u = chr_coord;                              // should equal zero
  memcpy(&chr_index[chr_index_pos], u32.c, 4);
  chr_index_pos += 4;
  u64.u = body.file_bufferPos;                    // should equal zero
  memcpy(&chr_index[chr_index_pos], u64.c, 8);
  chr_index_pos += 8;
  
  return 0;
}

int covFile::WriteEntry(unsigned int seqID, int value, unsigned int length) {
//  Assume sorted output. seqID is only there to signify if chromosome has changed
  stream_int32 i32;
  stream_uint32 u32;
  stream_uint64 u64;
  
  if(seqID >= chr_names.size() * 3) {
    return -1;
  } else if(seqID < out_cur_seqID) {
    return -1;      // Does not support non-sorted writing
  }
  
  // Assume non contiguous seq_ID... have to write empty entries for all seq_ID in between
  while(out_cur_seqID < seqID) {
    // Finalise current chromosome if exists
    
    // Write size of chr_index
    u32.u = chr_index_pos;
    write(u32.c, 4);
    // Now write the rest of the index buffer
    if(chr_index_pos > 0) {
      write(chr_index, chr_index_pos);
    }
    body.WriteBuffer();
    
    out_cur_seqID += 1;
    
    // Initialize current chromosome here
    
    free(chr_index);
    chr_index = (char*)malloc(1200);
    chr_index_alloc = 1200;
    chr_index_pos = 0;
    chr_coord = 0;
// Write first chr index entry:
    u32.u = chr_coord;
    memcpy(&chr_index[chr_index_pos], u32.c, 4);
    chr_index_pos += 4;
    u64.u = body.file_bufferPos;
    memcpy(&chr_index[chr_index_pos], u64.c, 8);
    chr_index_pos += 8;    
  }
  
  // Now seqID == out_cur_seqID
  // body should not be full if new chromosome entry. This ensures
  if(body.BufferIsFull(8)) {
    // i.e. if writing an extra 8 bytes will overflow buffer
    body.WriteBuffer();     // write bgzf block
    
    // If chr index buffer too full then reallocate more memory
    // Checking equality suffices as each memory allocation is in a block divisible by 12
    if(chr_index_pos == chr_index_alloc) {
      char * data_tmp;
      chr_index = (char*)realloc((data_tmp = chr_index), chr_index_alloc += 1200);
    }
    
    // Write file offset of new block into the index - 12 bytes per record
    u32.u = chr_coord;      // The begin coordinate of the next bgzf block
    memcpy(&chr_index[chr_index_pos], u32.c, 4);
    chr_index_pos += 4;
    u64.u = body.file_bufferPos;        // the file offset of the next bgzf block
    memcpy(&chr_index[chr_index_pos], u64.c, 8);
    chr_index_pos += 8;
  }

  i32.i = value;
  body.write(i32.c, 4);
  
  u32.u = length;
  body.write(u32.c, 4);
  
  chr_coord += length;

  return(0);
}


