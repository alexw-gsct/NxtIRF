

// GZ File writer

#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include "includedefine.h"

#define CHUNK_gz 262144

class GZWriter {
private:
  ostream * OUT;

  char compressed_buffer[CHUNK_gz];
  
  char buffer[CHUNK_gz];
  unsigned int bufferPos;
  
public:
  GZWriter();
  void SetOutputHandle(std::ostream *out_stream);
  int writebuffer(char * src, unsigned int len);
  int writeline(std::string s_src);
  int flush(int final = 0);
};
