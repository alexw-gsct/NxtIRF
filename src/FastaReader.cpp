#include "FastaReader.h"

void FastaReader::SetInputHandle(std::istream *in_stream) {
  IN = in_stream;
  FirstSeq = true;
}

bool FastaReader::ReadSeq() {
  std::string myLine;
  std::string sequence_raw;
  std::string line;
  std::string subline;
  std::string subline2;
  
  sequence.clear();
  if(FirstSeq) {
    std::getline(*IN, myLine, '>');
    FirstSeq = false;
  }
  std::getline(*IN, seqname, '\n');
  std::getline(*IN, sequence_raw, '>');
  
  std::stringstream sn(sequence_raw);
  while(std::getline(sn, line, '\n')){
    // remove /r
    std::stringstream sl(line);
    while(std::getline(sl, subline, '\r')){
    // append subline to sequence
    std::stringstream sl2(subline);
      while(std::getline(sl2, subline2, ' ')){
        sequence.append(subline2);
      }
    }
  }
}