/* GZReader.cpp Reads Gzipped files

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

#include "GZReader.h"
  
GZReader::GZReader() {
  bufferLen = 0;
  bufferPos = 0;
  buffer = NULL;
  
  loaded = false; lazy = false; streamed = false;
}

GZReader::~GZReader() {
  if(buffer != NULL) {
    free(buffer);
  }
}

// Only allowed for lazy = true
// Returns -1 if error; 0 if success; 1 if success and EOF
int GZReader::getline(std::string & s_myLine, const char delim) {
  if(!lazy || !loaded || streamed) return(-1);
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

// gets a single chunk of data and appends to buffer. Only for lazy = TRUE
int GZReader::GetBuffer() {
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
        Rcout << "Exception during zlib decompression: (" << err << ") " << error_string;
        free(data);
				return(err);
      }
    }
  }

  char *buffer_tmp;
  buffer = (char*)realloc(buffer_tmp = buffer, bufferLen + curpos);
  memcpy(&buffer[bufferLen], data, curpos);

  bufferLen += curpos;
  free(data);
  if (gzeof (gz_in)) {
    return(1);
  } else {
    return(0);
  }
}

// Loads a file
//   Options:
//   - lazy = FALSE: opens file as well as reads entire file into memory
//   - asStream = TRUE: copies read data into istringstream object
int GZReader::LoadGZ(std::string s_filename, bool asStream, bool lazymode) {
  gz_in = gzopen(s_filename.c_str(), "r");
  
  if(lazymode == false) {
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
            Rcout << "Exception during zlib decompression: (" << err << ") " << error_string;
            free(data);
						return(err);
          }
        }
      }
    }
    if(asStream) {
      iss.str((char*)data);
      loaded = true; streamed = true; lazy = false;
    } else {
      char *buffer_tmp;
      buffer = (char*)realloc(buffer_tmp = buffer, curpos);
      memcpy(buffer, data, curpos);
      bufferLen = curpos;
      loaded = true; streamed = false; lazy = false;
    }
    gzclose(gz_in);
    free(data);
  } else {
    loaded = true; streamed = false; lazy = true;
  }
	return(0);
}

// Operations to read gzip as if it was ifstream object. Not used in IRFinder

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

// Only required for lazy mode
int GZReader::closeGZ() {
	gzclose(gz_in);
	return(0);
}