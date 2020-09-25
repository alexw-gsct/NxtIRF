#include "GZReader.h"
#include <stdexcept>

#include "RcppArmadillo.h"

  using namespace Rcpp;
  
GZReader::GZReader() {
  bufferLen = 0;
  bufferPos = 0;
}

GZReader::~GZReader() {
  delete[] buffer;
}

void GZReader::LoadGZ(std::string s_filename, bool asStream) {
  gzFile gz_in;
  gz_in = gzopen(s_filename.c_str(), "r");
  
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
    Rcout << "Bytes read " << bytes_read << ", curpos " << curpos << "\n";
    
    if (bytes_read < CHUNK_gz - 1) {
      if (gzeof (gz_in)) {
        data = (unsigned char *)realloc((data_tmp = data), data_alloc -= (CHUNK_gz - 1) - bytes_read );
        break;
      }
      else {
        const char * error_string;
        error_string = gzerror (gz_in, & err);
        if (err) {
          std::ostringstream oss;
          oss << "Exception during zlib decompression: (" << err << ") " << error_string;
          throw(std::runtime_error(oss.str()));
        }
      }
    }
  }
  if(asStream) {
    iss.str((char*)data);
    gzclose(gz_in);
  } else {
    buffer = new char[curpos];
    memcpy(buffer, data, curpos);
  }
  free(data);
}

void GZReader::read(char * dest, const size_t len) {
  memcpy(dest, &buffer[bufferPos], len);
  bufferPos += len;
}
void GZReader::ignore(const size_t len) {
  bufferPos += len;
}
