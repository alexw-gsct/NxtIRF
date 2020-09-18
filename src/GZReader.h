// GZ File reader

#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <zconf.h>
#include <cstdlib>
#include "includedefine.h"

#define CHUNK_gz 1024

class GZReader {
private:

public:
  GZReader();
  ~GZReader();
  void LoadGZ(std::string s_filename, bool asStream = false);
  void read(char * dest, const size_t len);
  void ignore(const size_t len);
  std::istringstream iss;
  char * buffer;
  size_t bufferLen;
  size_t bufferPos;
};
