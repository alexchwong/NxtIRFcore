#include "GZWriter.h"
#include <stdexcept>

GZWriter::GZWriter() {
  bufferPos = 0;
}
void GZWriter::SetOutputHandle(std::ostream *out_stream) {
  OUT = out_stream;
}

int GZWriter::writeline(const std::string& s_src) {
  unsigned int s_size = s_src.size() + 1;
  char * line = new char[s_size];
  memcpy(line, s_src.data(), s_size - 1);
  line[s_size - 1] = '\n';
  writebuffer(line, s_size);
  delete[] line;
  return(0);
}

int GZWriter::writestring(const std::string& s_src) {
  unsigned int s_size = s_src.size();
  char * line = new char[s_size];
  memcpy(line, s_src.data(), s_size);
  writebuffer(line, s_size);
  delete[] line;
  return(0);
}

int GZWriter::writebuffer(const char * src, unsigned int len) {
  unsigned int bytesremaining = len;
  unsigned int srcpos = 0;
  if(bufferPos >= CHUNK_gz) {
    flush(0);
  }  
  while (bytesremaining + bufferPos > CHUNK_gz) {
    memcpy(&buffer[bufferPos], &src[srcpos], CHUNK_gz - bufferPos);
    srcpos += CHUNK_gz - bufferPos;
    bytesremaining -= CHUNK_gz - bufferPos;
    bufferPos = CHUNK_gz;
    flush(0);
  }
  memcpy(&buffer[bufferPos], &src[srcpos], bytesremaining);
  bufferPos += bytesremaining;
  bytesremaining = 0;
  if(bufferPos >= CHUNK_gz) {
    flush(0);
  }  
  return(0);
}

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
      std::ostringstream oss;
      oss << "Exception during zlib initialization: (" << ret << ") "  << strm.msg;
      throw(std::runtime_error(oss.str()));
    }

    strm.avail_in = bufferPos;
    strm.next_in = (Bytef*)buffer;
    strm.avail_out = CHUNK_gz;
    strm.next_out = (Bytef*)compressed_buffer;
    
    ret = deflate(&strm, Z_FINISH);  

    if (ret != Z_OK && ret != Z_STREAM_END) {
        std::ostringstream oss;
        oss << "Exception during zlib deflate: (" << ret << ") " << strm.msg;
        throw(std::runtime_error(oss.str()));
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