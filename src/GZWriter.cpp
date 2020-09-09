#include "GZWriter.h"
  
#include "RcppArmadillo.h"

GZWriter::GZWriter() {
  bufferPos = 0;
//  strm = new z_stream;
}

void GZWriter::SetOutputHandle(std::ostream *out_stream) {
  OUT = out_stream;
}

int GZWriter::writeline(std::string s_src) {
  unsigned int s_size = s_src.size();
  char line[s_size + 1];
  strcpy(line, s_src.c_str());
  line[s_size] = '\n';
  writebuffer(line, s_size + 1);
  return(0);
}

int GZWriter::writebuffer(char * src, unsigned int len) {
  if(len + bufferPos > CHUNK_gz) {
    flush(0);
  }  
  memcpy(&buffer[bufferPos], src, len);
  bufferPos += len;
  return(0);
}

int GZWriter::flush(int final) {
  int ret;
  unsigned int have;
  z_stream strm;
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  
  deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 31, 8, Z_DEFAULT_STRATEGY);

  strm.avail_in = bufferPos;
  strm.next_in = (Bytef*)buffer;
  strm.avail_out = CHUNK_gz;
  strm.next_out = (Bytef*)compressed_buffer;
  
  if(final == 0) {
    ret = deflate(&strm, Z_FINISH);    /* no bad return value */
  } else {
    ret = deflate(&strm, Z_FINISH);    /* no bad return value */
  }
  have = CHUNK_gz - strm.avail_out;

  OUT->write(compressed_buffer, have);
  OUT->flush();
  deflateEnd(&strm);
  bufferPos=0;
  return(Z_OK);
}
