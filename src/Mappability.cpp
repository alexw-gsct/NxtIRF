#include "RcppArmadillo.h"
#include "Mappability.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

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
  memcpy(&new_read[0], input_read, read_len);
  
  if(direction == 1) {
    string s_new_read = reverse_complement(string(new_read));
    memcpy(new_read, s_new_read.c_str(), s_new_read.length());
  }
  
  char error_nuc;
  switch(error_seed % 2) {
  case 0:
    switch(new_read[error_pos - 1]) {
    case 'A':
      error_nuc = 'G';
    case 'C':
      error_nuc = 'A';
    case 'G':
      error_nuc = 'T';
    case 'T':
      error_nuc = 'C';
    }
  case 1:
    switch(new_read[error_pos - 1]) {
    case 'A':
      error_nuc = 'T';
    case 'C':
      error_nuc = 'G';
    case 'G':
      error_nuc = 'C';
    case 'T':
      error_nuc = 'A';
    }    
  case 2:
    switch(new_read[error_pos - 1]) {
    case 'A':
      error_nuc = 'C';
    case 'C':
      error_nuc = 'T';
    case 'G':
      error_nuc = 'A';
    case 'T':
      error_nuc = 'G';
    }    
  }
  memcpy(&new_read[error_pos - 1], &error_nuc, 1);
  
  string return_str = string(new_read);
  delete[] new_read;
  return(string(new_read));
}

bool checkDNA(const std::string& strand) {
  unsigned short size = strand.size();
  for(unsigned short i=0; i<size; ++i) {
    if(strand[i]!='A' && strand[i]!='T' && strand[i]!='G' && strand[i]!='C') {
      return false;
    }
  }
  return true;
}

// [[Rcpp::export]]
int IRF_SupplyMappaRegionReads(std::string genome_file, std::string region_file, std::string out_fa, 
                               int read_len, int read_stride, int error_pos) {
  
  std::ifstream inGenome;
  inGenome.open(genome_file, std::ifstream::in);
  FastaReader inFA;
  inFA.SetInputHandle(&inGenome);
  
  std::ofstream outFA;
  outFA.open(out_fa, std::ofstream::binary);
  GZWriter outGZ;
  outGZ.SetOutputHandle(&outFA);
  
  std::ifstream inRegions;
  inRegions.open(region_file, std::ifstream::in);
  
  string myLine;
  string myField;
  std::map< std::string, std::vector< std::pair<unsigned int,unsigned int> > > region_list;
  while(!inRegions.eof() && !inRegions.fail()) {
    int start;
    int end;
    string s_chr;
    
    std::getline(inRegions, myLine, '\n');
    std::stringstream region_line(myLine);
    
    std::getline(region_line, s_chr, '\t');
    std::getline(region_line, myField, '\t');
    start = stol(myField);
    std::getline(region_line, myField, '\t');
    end = stol(myField);
    if(start == end) break;
    region_list[s_chr].push_back(std::make_pair(start, end));
  }
  inRegions.close();
  
  unsigned int direction = 0;
  char * read = new char[read_len + 1];
  unsigned int seed = 0;
  
  string chr;
  string sequence;
  
  
  while(!inGenome.eof() && !inGenome.fail()) {
    
    inFA.ReadSeq();
    sequence = inFA.sequence;
    chr = inFA.seqname;

    std::map< std::string, std::vector< std::pair<unsigned int,unsigned int> > >::iterator it_chr;
    it_chr = region_list.find(chr);

    
    if(it_chr != region_list.end()) {

      char * buffer = new char[sequence.length() + 1];
      std::strcpy (buffer, sequence.c_str());
      
      // Iterate through each element in vector
      for(auto it_region = it_chr->second.begin(); it_region!=it_chr->second.end(); it_region++) {

        for(unsigned int bufferPos = it_region->first; (bufferPos < it_region->second - read_len - 1); bufferPos += read_stride) {
          memcpy(read, &buffer[bufferPos - 1], read_len);
          
          if(checkDNA(string(read))) {
            std::string write_name;
            write_name = (direction == 0 ? ">RF!" : ">RR!");
            write_name.append(chr);
            write_name.append("!");
            write_name.append(std::to_string(bufferPos));
            outGZ.writeline(write_name);
            //        read_names.push_back(write_name);
            std::string write_seq = GenerateReadError(read, read_len, error_pos, direction, seed) ;
            outGZ.writeline(write_seq);
            //        read_seqs.push_back(write_seq);
            seed += 1;
            direction = (direction == 0 ? 1 : 0);
          }
           
          if((seed % 1000000 == 0) & (seed > 0)) {
            Rcout << "Processed " << bufferPos << " coord of chrom:" << chr << '\n';
          }
        }
       
      }
      delete[] buffer;
    }
  }
  delete[] read;
  
  inGenome.close();
  outGZ.flush(true);// outGZ.close();
  outFA.flush();
  outFA.close();
  
  return(0);
}


// [[Rcpp::export]]
int IRF_SupplyMappaReads(std::string genome_file, std::string out_fa, int read_len, int read_stride, int error_pos) {
  
  std::ifstream inGenome;
  inGenome.open(genome_file, std::ifstream::in);
  
  std::ofstream outFA;
  outFA.open(out_fa, std::ofstream::binary);
  GZWriter outGZ;
  outGZ.SetOutputHandle(&outFA);
    
  unsigned int direction = 0;
  char * read = new char[read_len + 1];
  unsigned int seed = 0;
  
  string myLine;
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
      if(checkDNA(string(read))) {
        std::string write_name;
        write_name = (direction == 0 ? ">RF!" : ">RR!");
        write_name.append(chr);
        write_name.append("!");
        write_name.append(std::to_string(bufferPos));
        outGZ.writeline(write_name);
//        read_names.push_back(write_name);
        std::string write_seq = GenerateReadError(read, read_len, error_pos, direction, seed) ;
        outGZ.writeline(write_seq);
//        read_seqs.push_back(write_seq);
        seed += 1;
        direction = (direction == 0 ? 1 : 0);
      }
      if((seed % 1000000 == 0) & (seed > 0)) {
        Rcout << "Processed " << bufferPos << " coord of chrom:" << chr << '\n';
      }
    }
    delete[] buffer;
  }
  delete[] read;
  
  inGenome.close();
  outGZ.flush(true); //outGZ.close();
  outFA.flush();
  outFA.close();
  return(0);
}

// [[Rcpp::export]]
int IRF_genmap(std::string bam_file, std::string output_path){
  std::string s_inBAM = bam_file;
  std::string outputDir = output_path;
  
  FragmentsMap oFragMap;
  
  BAM2blocks BB;
  
  BB.registerCallbackChrMappingChange( std::bind(&FragmentsMap::ChrMapUpdate, &oFragMap, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&FragmentsMap::ProcessBlocks, &oFragMap, std::placeholders::_1) );
  
  BAMReader inbam;
  std::ifstream inbam_stream;
  inbam_stream.open(s_inBAM, std::ifstream::binary);
  inbam.SetInputHandle(&inbam_stream);
  
  BB.openFile(&inbam); // This file needs to be a decompressed BAM. (setup via fifo / or expect already decompressed via stdin).
  
  std::string BBreport;
  BB.processAll(BBreport);
  
  std::ofstream outFragsMap;
  outFragsMap.open(outputDir + "/Mappability.txt", std::ifstream::out);
  oFragMap.WriteOutput(&outFragsMap);
  outFragsMap.flush(); outFragsMap.close();
  
  return(0);
}
