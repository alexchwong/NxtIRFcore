/* GZWriter.cpp Writes Gzipped files

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

#include "GZWriter.h"

void GZWriter::SetOutputHandle(std::ostream *out_stream) {
  OUT = out_stream;
}

// Writes a line to a gzipped file, including '\n'
int GZWriter::writeline(const std::string& s_src) {
  unsigned int s_size = s_src.size() + 1;
  char * line = new char[s_size];
  memcpy(line, s_src.data(), s_size - 1);
  line[s_size - 1] = '\n';
  int ret = writebuffer(line, s_size);
  delete[] line;
  return(ret);
}

// Writes a string to a gzipped file, excluding '\n'
int GZWriter::writestring(const std::string& s_src) {
  unsigned int s_size = s_src.size();
  char * line = new char[s_size];
  memcpy(line, s_src.data(), s_size);
  int ret = writebuffer(line, s_size);
  delete[] line;
  return(ret);
}

// Writes from given char* src buffer of given length len. Used internally in IRFinder
int GZWriter::writebuffer(const char * src, unsigned int len) {
  unsigned int bytesremaining = len;
  unsigned int srcpos = 0;
  int ret;
  if(bufferPos >= CHUNK_gz) {
    ret = flush(0);
    if(ret != Z_OK) return(ret);
  }  
  while (bytesremaining + bufferPos > CHUNK_gz) {
    memcpy(&buffer[bufferPos], &src[srcpos], CHUNK_gz - bufferPos);
    srcpos += CHUNK_gz - bufferPos;
    bytesremaining -= CHUNK_gz - bufferPos;
    bufferPos = CHUNK_gz;
    ret = flush(0);
    if(ret != Z_OK) return(ret);
  }
  memcpy(&buffer[bufferPos], &src[srcpos], bytesremaining);
  bufferPos += bytesremaining;
  bytesremaining = 0;
  if(bufferPos >= CHUNK_gz) {
    ret = flush(0);
    if(ret != Z_OK) return(ret);
  }  
  return(0);
}

// Writes from memory to file via ostream
// Returns Z_OK if success, and error message otherwise
int GZWriter::flush(bool final) {
  if(bufferPos > 0) {
    int ret;
    unsigned int have;
    z_stream strm;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    
    ret = deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 31, 8, Z_DEFAULT_STRATEGY);
    if (ret != Z_OK) {
      Rcout << "Exception during zlib initialization: (" << ret << ") "  << strm.msg;
			return(ret);
    }

    strm.avail_in = bufferPos;
    strm.next_in = (Bytef*)buffer;
    strm.avail_out = CHUNK_gz;
    strm.next_out = (Bytef*)compressed_buffer;
    
    ret = deflate(&strm, Z_FINISH);  

    if (ret != Z_OK && ret != Z_STREAM_END) {
        Rcout << "Exception during zlib deflate: (" << ret << ") " << strm.msg;
				return(ret);
    }

    have = strm.total_out;
  
    OUT->write(compressed_buffer, have);
    if(final) {   
      OUT->flush();
    }
    
    deflateEnd(&strm);
    bufferPos=0;
  }
  return(Z_OK);
}