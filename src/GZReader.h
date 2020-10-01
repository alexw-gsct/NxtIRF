// GZ File reader

#include <stdexcept>
#include "includedefine.h"

#define CHUNK_gz 262144

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
