// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

#include "ReadBlockProcessor.h"
#include "ReadBlockProcessor_CoverageBlocks.h"
#include "ReadBlockProcessor_OutputBAM.h"
#include "BAM2blocks.h"
#include "includedefine.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
int IRF_main(std::string bam_file, std::string reference_file, std::string output_file){
  
  std::string s_bam = bam_file;
  
  std::string s_output = output_file;
  
  std::string s_ref = reference_file;
  
  // std::string s_inCoverageBlocks = reference_path + "/ref-cover.bed";
  // std::string s_inSJ = reference_path + "/ref-sj.ref";
  // std::string s_inSpansPoint = reference_path + "/ref-read-continues.ref";
  // std::string s_inROI = reference_path + "/ref-ROI.bed";

    // Read single reference file:
  
  std::ifstream inRef;
  inRef.open(s_ref, std::ifstream::in);

    std::string myLine;
    std::string myBuffer;
    
    getline(inRef, myLine, '>');    // discard first >
    getline(inRef, myLine, '\n');   // ignore file names for now
    getline(inRef, myBuffer, '>');  // this is the data block for ref-cover.bed

  CoverageBlocksIRFinder oCoverageBlocks;
  std::istringstream inCoverageBlocks;
  inCoverageBlocks.str(myBuffer);
  oCoverageBlocks.loadRef(inCoverageBlocks);

    getline(inRef, myLine, '\n');
    getline(inRef, myBuffer, '>');

  SpansPoint oSpansPoint;
  oSpansPoint.setSpanLength(5,4);
  std::istringstream inSpansPoint;
  inSpansPoint.str(myBuffer);
  oSpansPoint.loadRef(inSpansPoint);

    getline(inRef, myLine, '\n');
    getline(inRef, myBuffer, '>');
  
  FragmentsInROI oFragmentsInROI;
  FragmentsInChr oFragmentsInChr;

    std::istringstream inFragmentsInROI;
    inFragmentsInROI.str(myBuffer);
    oFragmentsInROI.loadRef(inFragmentsInROI);

    getline(inRef, myLine, '\n');
    getline(inRef, myBuffer, '>');

  JunctionCount oJuncCount;
  std::istringstream inJuncCount;
  inJuncCount.str(myBuffer);
  oJuncCount.loadRef(inJuncCount);
  // inJuncCount.close();
  inRef.close();  

  // FragmentsInROI oFragmentsInROI;
  // FragmentsInChr oFragmentsInChr;
  
  // JunctionCount oJuncCount;
  // std::ifstream inJuncCount;
  // inJuncCount.open(s_inSJ, std::ifstream::in);
  // oJuncCount.loadRef(inJuncCount);
  // inJuncCount.close();
  
  // SpansPoint oSpansPoint;
  // oSpansPoint.setSpanLength(5,4);
  // std::ifstream inSpansPoint;
  // inSpansPoint.open(s_inSpansPoint, std::ifstream::in);
  // oSpansPoint.loadRef(inSpansPoint);
  // inSpansPoint.close();
  
  // CoverageBlocksIRFinder oCoverageBlocks;
  // std::ifstream inCoverageBlocks;
  // inCoverageBlocks.open(s_inCoverageBlocks, std::ifstream::in);
  // oCoverageBlocks.loadRef(inCoverageBlocks);
  // inCoverageBlocks.close();
  
  BAM2blocks BB;
  
  BB.registerCallbackChrMappingChange( std::bind(&JunctionCount::ChrMapUpdate, &oJuncCount, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&JunctionCount::ProcessBlocks, &oJuncCount, std::placeholders::_1) );
  
  BB.registerCallbackChrMappingChange( std::bind(&FragmentsInChr::ChrMapUpdate, &oFragmentsInChr, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&FragmentsInChr::ProcessBlocks, &oFragmentsInChr, std::placeholders::_1) );
  
  BB.registerCallbackChrMappingChange( std::bind(&SpansPoint::ChrMapUpdate, &oSpansPoint, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&SpansPoint::ProcessBlocks, &oSpansPoint, std::placeholders::_1) );
  
  // if (s_inROI != "NULL") {	
    // std::ifstream inFragmentsInROI;
    // inFragmentsInROI.open(s_inROI, std::ifstream::in);
    // oFragmentsInROI.loadRef(inFragmentsInROI);
    // inFragmentsInROI.close();
    
    BB.registerCallbackChrMappingChange( std::bind(&FragmentsInROI::ChrMapUpdate, &oFragmentsInROI, std::placeholders::_1) );
    BB.registerCallbackProcessBlocks( std::bind(&FragmentsInROI::ProcessBlocks, &oFragmentsInROI, std::placeholders::_1) );
  // }
  
  BB.registerCallbackChrMappingChange( std::bind(&CoverageBlocks::ChrMapUpdate, &oCoverageBlocks, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&CoverageBlocks::ProcessBlocks, &oCoverageBlocks, std::placeholders::_1) );
  
  BAMReader inbam;
  std::ifstream inbam_stream;
  inbam_stream.open(s_bam, std::ifstream::binary);
  inbam.SetInputHandle(&inbam_stream);
  
  BB.openFile(&inbam); // This file needs to be a decompressed BAM. (setup via fifo / or expect already decompressed via stdin).
  
  BB.processAll();

// Write output to file:  
  std::ofstream out;
  out.open(s_output, std::ifstream::out);

  std::ostringstream outFragmentsInROI;
  oFragmentsInROI.WriteOutput(&outFragmentsInROI);
  out << "ROIname\ttotal_hits\tpositive_strand_hits\tnegative_strand_hits\n" << outFragmentsInROI.str() << "\n";
  
  std::ostringstream outJuncCount;
  oJuncCount.WriteOutput(&outJuncCount);
  out << "JC_seqname\tstart\tend\tstrand\ttotal\tpos\tneg\n" << outJuncCount.str() << "\n";
  
  int directionality = oJuncCount.Directional();
  out << ">Directionality\n" << directionality << "\n\n";  
  
  std::ostringstream outSpansPoint;
  oSpansPoint.WriteOutput(&outSpansPoint);
  out << "SP_seqname\tpos\ttotal\tpos\tneg\n" << outSpansPoint.str() << "\n";
  
  std::ostringstream outFragmentsInChr;
  oFragmentsInChr.WriteOutput(&outFragmentsInChr);
  out << "ChrCoverage_seqname\ttotal\tpos\tneg\n" << outFragmentsInChr.str() << "\n";
  
  std::ostringstream outCoverageBlocks_nondir;
  oCoverageBlocks.WriteOutput(&outCoverageBlocks_nondir, oJuncCount, oSpansPoint);
  out << "NonDir_" << outCoverageBlocks_nondir.str() << "\n";
  
  if (directionality != 0) {
    std::ostringstream outCoverageBlocks_dir;
    oCoverageBlocks.WriteOutput(&outCoverageBlocks_dir, oJuncCount, oSpansPoint, directionality); // Directional.
    out << "Directional_" << outCoverageBlocks_nondir.str() << "\n";
  }
  
  out.close();
  
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