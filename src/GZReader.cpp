#include "GZReader.h"
#include <stdexcept>

#include "RcppArmadillo.h"

  using namespace Rcpp;
  
GZReader::GZReader() {
}

void GZReader::LoadGZ(std::string s_filename) {
  gzFile gz_in;
  gz_in = gzopen(s_filename.c_str(), "r");
  
  unsigned char *data = NULL;
  int data_alloc = 0;
  int curpos = 0;
  
  while(true) {
    int err;
    int bytes_read;
    unsigned char buffer[CHUNK_gz];
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
          std::ostringstream oss;
          oss << "Exception during zlib decompression: (" << err << ") " << error_string;
          throw(std::runtime_error(oss.str()));
        }
      }
    }
  }
  
  iss.str((char*)data);
  gzclose(gz_in);
  free(data);
}

