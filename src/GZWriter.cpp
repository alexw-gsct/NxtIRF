#include "GZWriter.h"
#include <stdexcept>

#include "RcppArmadillo.h"

GZWriter::GZWriter() {
  bufferPos = 0;
}
void GZWriter::SetOutputHandle(std::ostream *out_stream) {
  OUT = out_stream;
}
/*
int GZWriter::Open(const std::string s_file) {
  gz_out = gzopen(s_file.c_str(), "w");
}
*/
int GZWriter::writeline(const std::string s_src) {
  unsigned int s_size = s_src.size();
  char line[s_size + 1];
  strcpy(line, s_src.data());
  line[s_size] = '\n';
  writebuffer(line, s_size + 1);
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
      oss << "Exception during zlib initialization: (" << ret << ") ";
      throw(std::runtime_error(oss.str()));
    }
    
  
    strm.avail_in = bufferPos;
    strm.next_in = (Bytef*)buffer;
    strm.avail_out = CHUNK_gz;
    strm.next_out = (Bytef*)compressed_buffer;
    
//    if(final == 0) {
      ret = deflate(&strm, Z_FINISH);  
//    } else {
//      ret = deflate(&strm, Z_FINISH);
//    }
  if (ret != Z_OK) {
    std::ostringstream oss;
    oss << "Exception during zlib deflate: (" << ret << ") ";
    throw(std::runtime_error(oss.str()));
  }

    have = strm.total_out;
  
    OUT->write(compressed_buffer, have);
    OUT->flush();
    deflateEnd(&strm);
    bufferPos=0;
  }
  return(Z_OK);
}