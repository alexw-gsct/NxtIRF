// GZ File reader

#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <zconf.h>
#include <cstdlib>
#include "includedefine.h"

#define CHUNK_gz 262144

class GZReader {
private:

public:
  GZReader();
  void LoadGZ(std::string s_filename);

  std::istringstream iss;
};
