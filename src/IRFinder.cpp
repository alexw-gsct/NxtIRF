// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

#include "ReadBlockProcessor.h"
#include "ReadBlockProcessor_CoverageBlocks.h"
#include "ReadBlockProcessor_OutputBAM.h"
#include "BAM2blocks.h"
#include "includedefine.h"
#include "GZReader.h"

#include <zlib.h>

union stream_uint32 {
  char c[4];
  uint32_t u;
};
union stream_int32 {
  char c[4];
  int32_t i;
};


using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List IRF_RLE_From_Cov(std::string s_in, int strand) {
  stream_uint32 u32;
  stream_int32 i32;
  char buffer[1000];
  std::string chrName;

  GZReader gz_in;
  gz_in.LoadGZ(s_in);
  gz_in.ignore(4);
  
  std::vector<std::string> chr_names;
  std::vector<int> chr_lens;
  
  gz_in.read(i32.c ,4);
  int n_chr = i32.i;
  
  for (unsigned int i = 0; i < n_chr; i++) {
    gz_in.read(i32.c ,4);
    gz_in.read(buffer , i32.i);
    chrName = string(buffer, i32.i-1);
    chr_names.push_back(chrName);
    gz_in.read(i32.c ,4);
    chr_lens.push_back(i32.i);
  }

  List RLEList;

  for(int j = 0; j < 3; j++) {
    for(int i = 0; i < n_chr; i++) {
      gz_in.read(u32.c ,4);
      unsigned int block_size = u32.u;
      if(j == strand) {
        
        std::vector<int> values;
        std::vector<unsigned int> lengths;
        unsigned int pos = 0;
        while( pos < block_size ) {
          gz_in.read(i32.c ,4);
          values.push_back(i32.i);
          gz_in.read(u32.c ,4);
          lengths.push_back(u32.u);
          pos += 8;
        }
        List RLE = List::create(
          _["values"] = values,
          _["lengths"] = lengths 
        );
        RLEList.push_back(RLE, chr_names[i]);
        
      } else {
        gz_in.ignore(block_size);
      }
    }
  }
  return(RLEList);
}

// [[Rcpp::export]]
int IRF_gunzip(std::string s_in, std::string s_out) {
  
  GZReader gz_in;
  gz_in.LoadGZ(s_in, true);
  
  std::ofstream out;
  out.open(s_out, std::ifstream::out);
  std::string myLine;
  
  while(!gz_in.iss.eof()) {
    getline(gz_in.iss, myLine, '\n');
    out << myLine << "\n";
  }
  out.flush(); out.close();
  
  return(0);
}

// [[Rcpp::export]]
int IRF_main(std::string bam_file, std::string reference_file, std::string output_file){
  
  std::string s_bam = bam_file;
  
  std::string s_output = output_file;
  
  std::string s_ref = reference_file;
  
  // Read single reference file:
  Rcout << "Running IRFinder on " << s_bam << "\nReference: " << reference_file << "\n";
  Rcout << "Output file: " << output_file << "\n\n";

  Rcout << "Reading reference file\n";
  
  GZReader gz_in;
  gz_in.LoadGZ(reference_file, true);

  std::string myLine;
  std::string myBuffer;
  
    getline(gz_in.iss, myLine, '>');    // discard first >
    getline(gz_in.iss, myLine, '\n');   // ignore file names for now
    getline(gz_in.iss, myBuffer, '>');  // this is the data block for ref-cover.bed

  CoverageBlocksIRFinder oCoverageBlocks;
  std::istringstream inCoverageBlocks;
  inCoverageBlocks.str(myBuffer);
  oCoverageBlocks.loadRef(inCoverageBlocks);

    getline(gz_in.iss, myLine, '\n');
    getline(gz_in.iss, myBuffer, '>');

  SpansPoint oSpansPoint;
  oSpansPoint.setSpanLength(5,4);
  std::istringstream inSpansPoint;
  inSpansPoint.str(myBuffer);
  oSpansPoint.loadRef(inSpansPoint);

    getline(gz_in.iss, myLine, '\n');
    getline(gz_in.iss, myBuffer, '>');
  
  FragmentsInROI oFragmentsInROI;
  FragmentsInChr oFragmentsInChr;

  std::istringstream inFragmentsInROI;
  inFragmentsInROI.str(myBuffer);
  oFragmentsInROI.loadRef(inFragmentsInROI);

    getline(gz_in.iss, myLine, '\n');
    getline(gz_in.iss, myBuffer, '>');

  JunctionCount oJuncCount;
  std::istringstream inJuncCount;
  inJuncCount.str(myBuffer);
  oJuncCount.loadRef(inJuncCount);
  
  FragmentsMap oFragMap;
  

  BAM2blocks BB;
  
  BB.registerCallbackChrMappingChange( std::bind(&JunctionCount::ChrMapUpdate, &oJuncCount, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&JunctionCount::ProcessBlocks, &oJuncCount, std::placeholders::_1) );
  
  BB.registerCallbackChrMappingChange( std::bind(&FragmentsInChr::ChrMapUpdate, &oFragmentsInChr, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&FragmentsInChr::ProcessBlocks, &oFragmentsInChr, std::placeholders::_1) );
  
  BB.registerCallbackChrMappingChange( std::bind(&SpansPoint::ChrMapUpdate, &oSpansPoint, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&SpansPoint::ProcessBlocks, &oSpansPoint, std::placeholders::_1) );
  
  BB.registerCallbackChrMappingChange( std::bind(&FragmentsInROI::ChrMapUpdate, &oFragmentsInROI, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&FragmentsInROI::ProcessBlocks, &oFragmentsInROI, std::placeholders::_1) );

  BB.registerCallbackChrMappingChange( std::bind(&CoverageBlocks::ChrMapUpdate, &oCoverageBlocks, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&CoverageBlocks::ProcessBlocks, &oCoverageBlocks, std::placeholders::_1) );

  BB.registerCallbackChrMappingChange( std::bind(&FragmentsMap::ChrMapUpdate, &oFragMap, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&FragmentsMap::ProcessBlocks, &oFragMap, std::placeholders::_1) );
  
  Rcout << "Processing BAM file\n";
    
  BAMReader inbam;
  std::ifstream inbam_stream;
  inbam_stream.open(s_bam, std::ifstream::binary);
  inbam.SetInputHandle(&inbam_stream);
  
  BB.openFile(&inbam); // This file needs to be a decompressed BAM. (setup via fifo / or expect already decompressed via stdin).
  
  std::ostringstream outBAMsummary;
  BB.processAll(&outBAMsummary);

// Write output to file:  
  Rcout << "Writing output file\n";

  std::ofstream out;
  out.open(s_output + ".txt.gz", std::ios::binary);

// GZ compression:
  GZWriter outGZ;
  outGZ.SetOutputHandle(&out);

// Write stats here:

  outGZ.writeline("BAM Report\tValue");
  myLine = outBAMsummary.str();
  outGZ.writebuffer(myLine.data(), myLine.size());
  outGZ.writeline("");

  std::ostringstream outDirectional;
  int directionality = oJuncCount.Directional(&outDirectional);
  outGZ.writeline("Directionality\tValue");
  myLine = outDirectional.str();
  outGZ.writebuffer(myLine.data(), myLine.size());
  outGZ.writeline("");
//  Rcout << ">Directionality\n" << directionality << "\n\n";  
    
  std::ostringstream outFragmentsInROI;
  oFragmentsInROI.WriteOutput(&outFragmentsInROI);
  outGZ.writeline("ROIname\ttotal_hits\tpositive_strand_hits\tnegative_strand_hits");
  myLine = outFragmentsInROI.str();
  outGZ.writebuffer(myLine.data(), myLine.size());
  outGZ.writeline("");
//  Rcout << "ROIname\ttotal_hits\tpositive_strand_hits\tnegative_strand_hits\n" << outFragmentsInROI.str() << "\n";
  
  std::ostringstream outJuncCount;
  oJuncCount.WriteOutput(&outJuncCount);
  Rcout << "JuncCount written to ostringstream\n";
  outGZ.writeline("JC_seqname\tstart\tend\tstrand\ttotal\tpos\tneg");
  Rcout << "JuncCount header written\n";
  myLine = outJuncCount.str();
  outGZ.writebuffer(myLine.data(), myLine.size());
  outGZ.writeline("");
  Rcout << "JuncCount body written\n";
  
  outGZ.writeline("");
  /*
  myLine = outJuncCount.str();
  outGZ.writebuffer(myLine.data(), myLine.size());
  */
//  Rcout << "JC_seqname\tstart\tend\tstrand\ttotal\tpos\tneg\n" << outJuncCount.str() << "\n";
  
  std::ostringstream outSpansPoint;
  oSpansPoint.WriteOutput(&outSpansPoint);
  outGZ.writeline("SP_seqname\tpos\ttotal\tpos\tneg");
  myLine = outSpansPoint.str();
  outGZ.writebuffer(myLine.data(), myLine.size());
  outGZ.writeline("");
  
//  out << "SP_seqname\tpos\ttotal\tpos\tneg\n" << outSpansPoint.str() << "\n";
  
  std::ostringstream outFragmentsInChr;
  oFragmentsInChr.WriteOutput(&outFragmentsInChr);
  outGZ.writeline("ChrCoverage_seqname\ttotal\tpos\tneg");
  myLine = outFragmentsInChr.str();
  outGZ.writebuffer(myLine.data(), myLine.size());
  outGZ.writeline("");
//  out << "ChrCoverage_seqname\ttotal\tpos\tneg\n" << outFragmentsInChr.str() << "\n";
  
  std::ostringstream outCoverageBlocks_nondir;
  oCoverageBlocks.WriteOutput(&outCoverageBlocks_nondir, oJuncCount, oSpansPoint);
  myLine = outCoverageBlocks_nondir.str();
  outGZ.writebuffer(myLine.data(), myLine.size());
  outGZ.writeline("");
  //  out << "NonDir_" << outCoverageBlocks_nondir.str() << "\n";
  
  if (directionality != 0) {
    std::ostringstream outCoverageBlocks_dir;
    oCoverageBlocks.WriteOutput(&outCoverageBlocks_dir, oJuncCount, oSpansPoint, directionality); // Directional.
    myLine = outCoverageBlocks_dir.str();
    outGZ.writebuffer(myLine.data(), myLine.size());
    outGZ.writeline("");
//    out << "Directional_" << outCoverageBlocks_nondir.str() << "\n";
  }
  outGZ.flush(true); // outGZ.close();
  out.flush(); out.close();
  
  // Write Coverage Binary file:
  
  std::ofstream COVout;
  COVout.open(output_file + ".cov.gz", std::ofstream::binary);
   
  GZWriter outGZ2 = GZWriter(9);
  outGZ2.SetOutputHandle(&COVout);
  
  oFragMap.WriteBinary(&outGZ2, BB.chr_names, BB.chr_lens);
  outGZ2.flush(true); COVout.flush(); COVout.close();
  
  // if (s_inROI != "NULL") {
  /*  Output computed statistics from data structures.
    oFragmentsInROI -- this tells us if the data was directional or not -- if we need to know for other output modules. */
    // std::ofstream outFragmentsInROI;
    // outFragmentsInROI.open(s_output + "/IRFinder-ROI.txt", std::ifstream::out);
    // oFragmentsInROI.WriteOutput(&outFragmentsInROI);
    // outFragmentsInROI.flush(); outFragmentsInROI.close();
  // }
  
  // std::ofstream outJuncCount;
  // outJuncCount.open(s_output + "/IRFinder-JuncCount.txt", std::ifstream::out);
  // oJuncCount.WriteOutput(&outJuncCount);
  // outJuncCount.flush(); outJuncCount.close();
  
  // int directionality = oJuncCount.Directional();
  // cout << "RNA-Seq directionality -1/0/+1:\t" << directionality << "\n";
  
  
  // std::ofstream outSpansPoint;
  // outSpansPoint.open(s_output + "/IRFinder-SpansPoint.txt", std::ifstream::out);
  // oSpansPoint.WriteOutput(&outSpansPoint);
  // outSpansPoint.flush(); outSpansPoint.close();
  
  // std::ofstream outFragmentsInChr;
  // outFragmentsInChr.open(s_output + "/IRFinder-ChrCoverage.txt", std::ifstream::out);
  // oFragmentsInChr.WriteOutput(&outFragmentsInChr);
  // outFragmentsInChr.flush(); outFragmentsInChr.close();
  
  // std::ofstream outCoverageBlocks;
  // outCoverageBlocks.open(s_output + "/IRFinder-IR-nondir.txt", std::ifstream::out);
  // oCoverageBlocks.WriteOutput(&outCoverageBlocks, oJuncCount, oSpansPoint);
  // outCoverageBlocks.flush(); outCoverageBlocks.close();
  
  // if (directionality != 0) {
    // outCoverageBlocks.open(s_output + "/IRFinder-IR-dir.txt", std::ifstream::out);
    // oCoverageBlocks.WriteOutput(&outCoverageBlocks, oJuncCount, oSpansPoint, directionality); // Directional.
    // outCoverageBlocks.flush(); outCoverageBlocks.close();
  // }
  
  return(0);
}