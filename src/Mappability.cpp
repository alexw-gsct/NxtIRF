#include "Mappability.h"
#include <cassert>

char c_complement(char n)
{   
  switch(n)
  {   
  case 'A':
    return 'T';
  case 'T':
    return 'A';
  case 'G':
    return 'C';
  case 'C':
    return 'G';
  }   
  assert(false);
  return ' ';
}

void reverseit(char arr[])
{
  int len= strlen(arr) - 1;
  for(int i=0; i<=len/2; i++)
  {
    char temp=arr[i];
    arr[i]=arr[len-i];
    arr[len-i]=temp;
  }
}

std::string reverse_complement(std::string sequence) {
  char * buffer = new char[sequence.length() + 1];
  strcpy(buffer, sequence.c_str());
  reverseit(buffer);
  string nucs = string(buffer);
  delete[] buffer;
  transform(
    begin(nucs),
    end(nucs),
    begin(nucs),
    c_complement);
  return(nucs);
}

std::string GenerateReadError(char * input_read, unsigned int read_len, unsigned int error_pos,
                              unsigned int direction, unsigned int error_seed) {
  
  char * new_read = new char[read_len + 1];
  new_read[read_len] = '\0';
  if(direction == 0) {
    memcpy(&new_read[0], input_read, read_len);  
  } else {
    for(unsigned int i = 0; i < read_len; i++) {
      switch(input_read[i])
      {   
      case 'A':
        new_read[read_len - i - 1] = 'T'; break;
      case 'T':
        new_read[read_len - i - 1] = 'A'; break;
      case 'G':
        new_read[read_len - i - 1] = 'C'; break;
      case 'C':
        new_read[read_len - i - 1] = 'G'; break;
      default :
        new_read[read_len - i - 1] = input_read[i];
      }         
    }
  }
  /*
  if(direction == 1) {
    string s_new_read = reverse_complement(string(new_read));
    memcpy(new_read, s_new_read.c_str(), s_new_read.length());
  }
  */
  char error_nuc;
  switch(error_seed % 2) {
  case 0:
    switch(new_read[error_pos - 1]) {
    case 'A':
      error_nuc = 'G'; break;
    case 'C':
      error_nuc = 'A'; break;
    case 'G':
      error_nuc = 'T'; break;
    case 'T':
      error_nuc = 'C'; break;
    }
  case 1:
    switch(new_read[error_pos - 1]) {
    case 'A':
      error_nuc = 'T'; break;
    case 'C':
      error_nuc = 'G'; break;
    case 'G':
      error_nuc = 'C'; break;
    case 'T':
      error_nuc = 'A'; break;
    }    
  case 2:
    switch(new_read[error_pos - 1]) {
    case 'A':
      error_nuc = 'C'; break;
    case 'C':
      error_nuc = 'T'; break;
    case 'G':
      error_nuc = 'A'; break;
    case 'T':
      error_nuc = 'G'; break;
    }    
  }
  memcpy(&new_read[error_pos - 1], &error_nuc, 1);
  
  string return_str = string(new_read);
  delete[] new_read;
  return(return_str);
}

bool checkDNA(char * input_read, unsigned int read_len) {
  for(unsigned int i = 0; i < read_len; i++) {
    if(input_read[i]!='A' && input_read[i]!='T' && input_read[i]!='G' && input_read[i]!='C') {
      return false;
    }
  }
  return true;
}

/*
bool checkDNA(const std::string& strand) {
  unsigned short size = strand.size();
  for(unsigned short i=0; i<size; ++i) {
    if(strand[i]!='A' && strand[i]!='T' && strand[i]!='G' && strand[i]!='C') {
      return false;
    }
  }
  return true;
}
*/

// [[Rcpp::export]]
int IRF_GenerateMappabilityReads(std::string genome_file, std::string out_fa,
	int read_len, int read_stride, int error_pos) {
  
  std::ifstream inGenome;
  inGenome.open(genome_file, std::ifstream::in);
  
  std::ofstream outFA;
  outFA.open(out_fa, std::ios::binary);
  // GZWriter outGZ;
  // outGZ.SetOutputHandle(&outFA);
    
  unsigned int direction = 0;
  char * read = new char[read_len + 1];
  unsigned int seed = 0;
  
  string chr;
  string sequence;

  FastaReader inFA;
  inFA.SetInputHandle(&inGenome);
  
  while(!inGenome.eof() && !inGenome.fail()) {

    inFA.ReadSeq();
    sequence = inFA.sequence;
    chr = inFA.seqname;
    char * buffer = new char[sequence.length() + 1];
    std::strcpy (buffer, sequence.c_str());
    
    for(unsigned int bufferPos = 1; (bufferPos < sequence.length() - read_len - 1); bufferPos += read_stride) {
      memcpy(read, &buffer[bufferPos - 1], read_len);
      if(checkDNA(read, read_len)) {
        std::string write_name;
        write_name = (direction == 0 ? ">RF!" : ">RR!");
        write_name.append(chr);
        write_name.append("!");
        write_name.append(std::to_string(bufferPos));

        // outGZ.writeline(write_name);
        outFA << write_name << '\n';
        std::string write_seq = GenerateReadError(read, read_len, error_pos, direction, seed) ;
        // outGZ.writeline(write_seq);
        outFA << write_seq << '\n';
        
        seed += 1;
        direction = (direction == 0 ? 1 : 0);
      }
      if((seed % 100000 == 0) & (seed > 0)) {
        Rcout << "Processed " << bufferPos << " coord of chrom:" << chr << '\n';
      }
    }
    delete[] buffer;
  }
  delete[] read;
  
  inGenome.close();
  // outGZ.flush(true);
  outFA.flush();
  outFA.close();
  return(0);
}

// [[Rcpp::export]]
int IRF_GenerateMappabilityRegions(std::string bam_file, std::string output_file, int threshold){
  std::string s_inBAM = bam_file;
  std::string s_outFile = output_file;
  
  FragmentsMap oFragMap;
  
  BAM2blocks BB;
  
  BB.registerCallbackChrMappingChange( std::bind(&FragmentsMap::ChrMapUpdate, &oFragMap, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&FragmentsMap::ProcessBlocks, &oFragMap, std::placeholders::_1) );
  
  BAMReader inbam;
  std::ifstream inbam_stream;
  inbam_stream.open(s_inBAM, std::ifstream::binary);
  inbam.SetInputHandle(&inbam_stream);
  
  BB.openFile(&inbam);
  
  std::string BBreport;
  BB.processAll(BBreport);
  
  std::ofstream outFragsMap;
  outFragsMap.open(s_outFile, std::ifstream::out);
  oFragMap.WriteOutput(&outFragsMap, BB.chr_names, BB.chr_lens, threshold);
  outFragsMap.flush(); outFragsMap.close();
  
  return(0);
}
