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

  BB.processAll();

  std::ofstream outFragsMap;
  outFragsMap.open(outputDir + "/Mappability.txt", std::ifstream::out);
  oFragMap.WriteOutput(&outFragsMap);
  outFragsMap.flush(); outFragsMap.close();

  return(0);
}

// [[Rcpp::export]]
int IRF_main(std::string bam_file, std::string reference_path, std::string output_path){
  
  std::string s_inBAM = bam_file;
  
  std::string outputDir = output_path;
  
  std::string s_inCoverageBlocks = reference_path + "/ref-cover.bed";
  std::string s_inSJ = reference_path + "/ref-sj.ref";
  std::string s_inSpansPoint = reference_path + "/ref-read-continues.ref";
  std::string s_inROI = reference_path + "/ref-ROI.bed";
  
  FragmentsInROI oFragmentsInROI;
  FragmentsInChr oFragmentsInChr;
  
  JunctionCount oJuncCount;
  std::ifstream inJuncCount;
  inJuncCount.open(s_inSJ, std::ifstream::in);
  oJuncCount.loadRef(inJuncCount);
  inJuncCount.close();
  
  SpansPoint oSpansPoint;
  oSpansPoint.setSpanLength(5,4);
  std::ifstream inSpansPoint;
  inSpansPoint.open(s_inSpansPoint, std::ifstream::in);
  oSpansPoint.loadRef(inSpansPoint);
  inSpansPoint.close();
  
  CoverageBlocksIRFinder oCoverageBlocks;
  std::ifstream inCoverageBlocks;
  inCoverageBlocks.open(s_inCoverageBlocks, std::ifstream::in);
  oCoverageBlocks.loadRef(inCoverageBlocks);
  inCoverageBlocks.close();
  
  BAM2blocks BB;
  
  BB.registerCallbackChrMappingChange( std::bind(&JunctionCount::ChrMapUpdate, &oJuncCount, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&JunctionCount::ProcessBlocks, &oJuncCount, std::placeholders::_1) );
  
  BB.registerCallbackChrMappingChange( std::bind(&FragmentsInChr::ChrMapUpdate, &oFragmentsInChr, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&FragmentsInChr::ProcessBlocks, &oFragmentsInChr, std::placeholders::_1) );
  
  BB.registerCallbackChrMappingChange( std::bind(&SpansPoint::ChrMapUpdate, &oSpansPoint, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&SpansPoint::ProcessBlocks, &oSpansPoint, std::placeholders::_1) );
  
  if (s_inROI != "NULL") {	
    std::ifstream inFragmentsInROI;
    inFragmentsInROI.open(s_inROI, std::ifstream::in);
    oFragmentsInROI.loadRef(inFragmentsInROI);
    inFragmentsInROI.close();
    
    BB.registerCallbackChrMappingChange( std::bind(&FragmentsInROI::ChrMapUpdate, &oFragmentsInROI, std::placeholders::_1) );
    BB.registerCallbackProcessBlocks( std::bind(&FragmentsInROI::ProcessBlocks, &oFragmentsInROI, std::placeholders::_1) );
  }
  
  BB.registerCallbackChrMappingChange( std::bind(&CoverageBlocks::ChrMapUpdate, &oCoverageBlocks, std::placeholders::_1) );
  BB.registerCallbackProcessBlocks( std::bind(&CoverageBlocks::ProcessBlocks, &oCoverageBlocks, std::placeholders::_1) );
  
  BAMReader inbam;
  std::ifstream inbam_stream;
  inbam_stream.open(s_inBAM, std::ifstream::binary);
  inbam.SetInputHandle(&inbam_stream);
  
  BB.openFile(&inbam); // This file needs to be a decompressed BAM. (setup via fifo / or expect already decompressed via stdin).
  
  BB.processAll();
  
  if (s_inROI != "NULL") {
    // Output computed statistics from data structures.
    //oFragmentsInROI -- this tells us if the data was directional or not -- if we need to know for other output modules.
    std::ofstream outFragmentsInROI;
    outFragmentsInROI.open(outputDir + "/IRFinder-ROI.txt", std::ifstream::out);
    oFragmentsInROI.WriteOutput(&outFragmentsInROI);
    outFragmentsInROI.flush(); outFragmentsInROI.close();
  }
  
  std::ofstream outJuncCount;
  outJuncCount.open(outputDir + "/IRFinder-JuncCount.txt", std::ifstream::out);
  oJuncCount.WriteOutput(&outJuncCount);
  outJuncCount.flush(); outJuncCount.close();
  
  int directionality = oJuncCount.Directional();
  cout << "RNA-Seq directionality -1/0/+1:\t" << directionality << "\n";
  
  
  std::ofstream outSpansPoint;
  outSpansPoint.open(outputDir + "/IRFinder-SpansPoint.txt", std::ifstream::out);
  oSpansPoint.WriteOutput(&outSpansPoint);
  outSpansPoint.flush(); outSpansPoint.close();
  
  std::ofstream outFragmentsInChr;
  outFragmentsInChr.open(outputDir + "/IRFinder-ChrCoverage.txt", std::ifstream::out);
  oFragmentsInChr.WriteOutput(&outFragmentsInChr);
  outFragmentsInChr.flush(); outFragmentsInChr.close();
  
  std::ofstream outCoverageBlocks;
  outCoverageBlocks.open(outputDir + "/IRFinder-IR-nondir.txt", std::ifstream::out);
  oCoverageBlocks.WriteOutput(&outCoverageBlocks, oJuncCount, oSpansPoint);
  outCoverageBlocks.flush(); outCoverageBlocks.close();
  
  if (directionality != 0) {
    outCoverageBlocks.open(outputDir + "/IRFinder-IR-dir.txt", std::ifstream::out);
    oCoverageBlocks.WriteOutput(&outCoverageBlocks, oJuncCount, oSpansPoint, directionality); // Directional.
    outCoverageBlocks.flush(); outCoverageBlocks.close();
  }
  
  return(1);
}