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
  std::getline(*IN, myLine, '\n');
  std::istringstream iss;
  iss.str(myLine);
  std::getline(iss, seqname, '\t');  // chromosome name is tab-terminated
  
  std::getline(*IN, sequence_raw, '>');

  sequence_raw.erase( std::remove(sequence_raw.begin(), sequence_raw.end(), ' '), sequence_raw.end() );
  sequence_raw.erase( std::remove(sequence_raw.begin(), sequence_raw.end(), '\r'), sequence_raw.end() );
  sequence_raw.erase( std::remove(sequence_raw.begin(), sequence_raw.end(), '\n'), sequence_raw.end() );
  sequence.append(sequence_raw);
  return(true);
}