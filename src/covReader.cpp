/* covReader.cpp Reads COV files (BAM coverage)

Copyright (C) 2021 Alex Chit Hei Wong

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.  */

#include "covReader.h"

// Constructor
covReader::covReader() {
  bufferPos = 0;
  bufferMax = 0;
  index_begin = 0;
  body_begin = 0;

  compressed_buffer = (char*)malloc(65536);
  buffer = (char*)malloc(65536);
}

// Destructor
covReader::~covReader() {
  free(buffer);
  free(compressed_buffer);
}

void covReader::SetInputHandle(std::istream *in_stream) {
  IS_EOF = 0;
  IS_FAIL = 0;
  IS_LENGTH = 0;

  IN = in_stream;

  // Identify EOF
  IN->seekg (0, std::ios_base::end);
  IS_LENGTH = IN->tellg();
  
  // Check EOF bit
  IN->seekg (-bamEOFlength, std::ios_base::end);
  
  char check_eof_buffer[bamEOFlength+1];
  IN->read(check_eof_buffer, bamEOFlength);
       
  if(strncmp(check_eof_buffer, bamEOF, bamEOFlength) == 0) {
    EOF_POS = IS_LENGTH - bamEOFlength;
  } else {
    // Rcout << "EOF bit not detected\n";
    EOF_POS = 0;
    IS_EOF = 1;
    IS_FAIL = 1;				
  }
  IN->seekg (0, std::ios_base::beg);
}

int covReader::ReadBuffer() {
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
    Rcout << "Exception during BAM decompression - BGZF header corrupt: (at " 
      << IN->tellg() << " bytes) ";
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
    Rcout << "Exception during BAM decompression - inflateInit2() fail: (" 
      << ret << ") ";
    return(ret);
  }
  ret = inflate(&zs, Z_FINISH);
  if(ret != Z_OK && ret != Z_STREAM_END) {
    Rcout << "Exception during BAM decompression - inflate() fail: (" 
      << ret << ") ";
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
    return(ret);
  }
  bufferPos = 0;
  
  return(ret);
}

int covReader::read(char * dest, unsigned int len) {
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

int covReader::ignore(unsigned int len) {
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

bool covReader::eof() {
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

bool covReader::fail() {
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

int covReader::ReadHeader() {
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
  
  // should be the start point of bgzf block containing index
  index_begin = IN->tellg();      
  
  bufferPos = 0;
  bufferMax = 0;    

  // keep profiling to identify where body_begin is
  stream_uint32 u32;
  for(unsigned int j = 0; j < 3; j++) {
    for(unsigned int i = 0; i < chr_names.size(); i++) {
      read(u32.c, 4);
      ignore(u32.u);            
    }
  }
  body_begin = IN->tellg();

  return(n_ref.u);
}

int covReader::GetChrs(std::vector<chr_entry> &chrs) {
  if(chr_names.size() > 0) {
    for(unsigned int i = 0; i < chr_names.size(); i++) {
      chrs.push_back(chr_entry(i, chr_names.at(i), chr_lens.at(i)));
    }
  }
  return(0);
}

int covReader::FetchPos(
  const std::string seqname, const uint32_t start, const int strand,
  uint64_t * file_offset, uint32_t * block_start
) {

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
    if(cur_block_start.u > start) { 
      // first entry should always equal zero so should not be immediately called
      break;
    } else {
      prev_block_start = cur_block_start.u;  
      prev_offset = cur_offset.u;
    }
  }

  *file_offset = prev_offset + body_begin;
  *block_start = prev_block_start;
  return 0;
}

int covReader::FetchRLE(
    const std::string seqname, 
    const uint32_t start, const uint32_t end, 
    const int strand,
    std::vector<int> * values, 
    std::vector<unsigned int> * lengths
) {
 
  stream_int32 i32;
  stream_uint32 u32;
  
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
