

// GZ File writer

#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include "includedefine.h"

#define CHUNK_gz 65536

class GZWriter {
private:
  ostream * OUT;
//  gzFile gz_out;
  
  char compressed_buffer[CHUNK_gz];
  
  char buffer[CHUNK_gz];
  unsigned int bufferPos;
  
public:
  GZWriter();
  void SetOutputHandle(std::ostream *out_stream);
//  int Open(const std::string s_file);
  int writebuffer(const char * src, unsigned int len);
  int writeline(const std::string s_src);
  int flush(bool final = false);
};
