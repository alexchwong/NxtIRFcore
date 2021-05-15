#include "GZReader.h"
  
GZReader::GZReader() {
  bufferLen = 0;
  bufferPos = 0;
  buffer = NULL;
}

GZReader::~GZReader() {
  if(buffer != NULL) {
    // delete[] buffer;
    free(buffer);
  }
}

int GZReader::getline(std::string & s_myLine, const char delim) {
  
  int ret = 0;
  unsigned long i = bufferPos;
  while(ret != 1) {
    if(i == bufferLen) {
      ret = GetBuffer();
    }
    while(i < bufferLen) {
      if(buffer[i] == delim) {
        break;
      }
      i++;
    }
    if(i < bufferLen || ret == 1) {
      s_myLine.clear();
      if(i > bufferPos) {
        char * str_ptr = new char[i - bufferPos + 1];
        memcpy(str_ptr, &buffer[bufferPos], i - bufferPos);
        str_ptr[i - bufferPos] = '\0';
        s_myLine.assign(str_ptr, i - bufferPos + 1);
        delete[] str_ptr;
      }
      bufferPos = i + 1;
			return(ret);
    }
  }
	return(ret);
}

int GZReader::GetBuffer() {
  // gets a chunk of data and appends to buffer
  unsigned char *data = NULL;
  int data_alloc = 0;
  int curpos = 0;
  
  int err;
  int bytes_read;
  unsigned char *data_tmp;
  
  data = (unsigned char *)realloc((data_tmp = data), data_alloc += CHUNK_gz - 1);
  bytes_read = gzread (gz_in, data + curpos, CHUNK_gz - 1);
  curpos += bytes_read;
  
  if (bytes_read < CHUNK_gz - 1) {
    if (gzeof (gz_in)) {
      data = (unsigned char *)realloc((data_tmp = data), data_alloc -= (CHUNK_gz - 1) - bytes_read );
    } else {
      const char * error_string;
      error_string = gzerror (gz_in, & err);
      if (err) {
        // std::ostringstream oss;
        Rcout << "Exception during zlib decompression: (" << err << ") " << error_string;
        // throw(std::runtime_error(oss.str()));
        free(data);
				return(err);
      }
    }
  }

  // if(bufferLen > 0) {
    // reallocate buffer with bufferLen + curpos
    char *buffer_tmp;
    buffer = (char*)realloc(buffer_tmp = buffer, bufferLen + curpos);
    memcpy(&buffer[bufferLen], data, curpos);
  // } else { 
    // buffer = (char*)realloc(buffer_tmp = buffer, curpos);
    // memcpy(buffer, data, curpos);
  // }
  bufferLen += curpos;
  free(data);
  if (gzeof (gz_in)) {
    return(1);
  } else {
    return(0);
  }
}

int GZReader::LoadGZ(std::string s_filename, bool asStream, bool lazy) {
  gz_in = gzopen(s_filename.c_str(), "r");
  
  if(lazy == false) {
    unsigned char *data = NULL;
    int data_alloc = 0;
    int curpos = 0;
    
    while(true) {
      int err;
      int bytes_read;
      unsigned char *data_tmp;
      
      data = (unsigned char *)realloc((data_tmp = data), data_alloc += CHUNK_gz - 1);
      bytes_read = gzread (gz_in, data + curpos, CHUNK_gz - 1);
      curpos += bytes_read;
      
      if (bytes_read < CHUNK_gz - 1) {
        if (gzeof (gz_in)) {
          data = (unsigned char *)realloc((data_tmp = data), data_alloc -= (CHUNK_gz - 1) - bytes_read );
          break;
        }
        else {
          const char * error_string;
          error_string = gzerror (gz_in, & err);
          if (err) {
            // std::ostringstream oss;
            Rcout << "Exception during zlib decompression: (" << err << ") " << error_string;
            // throw(std::runtime_error(oss.str()));
            free(data);
						return(err);
          }
        }
      }
    }
    if(asStream) {
      iss.str((char*)data);
    } else {
      char *buffer_tmp;
      buffer = (char*)realloc(buffer_tmp = buffer, curpos);
      memcpy(buffer, data, curpos);
      bufferLen = curpos;
    }
    gzclose(gz_in);
    free(data);
  } else {
    // lazy load just opens the file
  }
	return(0);
}

void GZReader::read(char * dest, const unsigned long len) {
  memcpy(dest, &buffer[bufferPos], len);
  bufferPos += len;
}
void GZReader::ignore(const unsigned long len) {
  bufferPos += len;
}

bool GZReader::eof() {
  return(gzeof(gz_in) && bufferPos == bufferLen);
}

int GZReader::closeGZ() {
	gzclose(gz_in);
	return(0);
}