#include "ReadBlockProcessor.h"
#include "ReadBlockProcessor_CoverageBlocks.h"
#include "BAM2blocks.h"
#include "GZReader.h"
#include "Mappability.h"

#include "includedefine.h"

#ifndef GALAXY

// [[Rcpp::export]]
List IRF_RLE_From_Cov(std::string s_in, std::string seqname, int start, int end, int strand) {
// Returns an RLE covering the region described above
// s_in: The coverage file
// strand: 0 = +, 1 = -, 2 = *
  
  List NULL_RLE = List::create(
    _["values"] = 0,
    _["lengths"] = 0 
  );
  
  if(start > end){
    return(NULL_RLE);
  }

  std::ifstream inCov_stream;
  inCov_stream.open(s_in, std::ifstream::binary);
  
  covFile inCov;
  inCov.SetInputHandle(&inCov_stream);
  
  inCov.ReadHeader();
  
  // Find corresponding seqname
  int ref_index;
  auto it_chr = std::find(inCov.chr_names.begin(), inCov.chr_names.end(), seqname);
  if(it_chr == inCov.chr_names.end()) {
    return(NULL_RLE);
  } else {
    ref_index = distance(inCov.chr_names.begin(), it_chr);
  }
  // end = 0 implies fetch whole chromosome
  int eff_end = 0;
  if(end == 0) {
    eff_end = inCov.chr_lens.at(ref_index);
  } else {
    eff_end = end;
  }
  
  std::vector<int> values;
  std::vector<unsigned int> lengths;
  // Push first value
  values.push_back(0);
  lengths.push_back((unsigned int)start);
  
  inCov.FetchRLE(seqname, (uint32_t)start, (uint32_t)eff_end, strand, &values, &lengths);

  inCov_stream.close();
  // Push last value
  if((uint32_t)eff_end < inCov.chr_lens.at(ref_index)) {
    values.push_back(0);
    lengths.push_back(inCov.chr_lens.at(ref_index) - eff_end);
  }
    
  List RLE = List::create(
    _["values"] = values,
    _["lengths"] = lengths 
  );

  return(RLE);
}

// [[Rcpp::export]]
List IRF_RLEList_From_Cov(std::string s_in, int strand) {
  // Returns an RLEList
  // s_in: The coverage file
  // strand: 0 = +, 1 = -, 2 = *
  
  List NULL_RLE = List::create(
    _["values"] = 0,
    _["lengths"] = 0 
  );
  
  List RLEList;
  
  std::ifstream inCov_stream;
  inCov_stream.open(s_in, std::ifstream::binary);
  
  covFile inCov;
  inCov.SetInputHandle(&inCov_stream);
  
  inCov.ReadHeader();
  
  for (unsigned int i = 0; i < inCov.chr_names.size(); i++) {
    uint32_t eff_end = inCov.chr_lens.at(i);
    uint32_t start = 0;
    
    std::vector<int> values;
    std::vector<unsigned int> lengths;

    inCov.FetchRLE(inCov.chr_names.at(i), (uint32_t)start, (uint32_t)eff_end, strand, &values, &lengths);
    
    List RLE = List::create(
      _["values"] = values,
      _["lengths"] = lengths 
    );
    RLEList.push_back(RLE, inCov.chr_names.at(i));
  }

  inCov_stream.close();
  
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

#else
	// galaxy
#endif

// [[Rcpp::export]]
int IRF_main(std::string bam_file, std::string reference_file, std::string output_file){
  
  std::string s_bam = bam_file;
  
  std::string s_output = output_file;
  
  std::string s_ref = reference_file;
  
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
  BB.processAll(myLine);

// Write output to file:  
  Rcout << "Writing output file\n";

  std::ofstream out;
  out.open(s_output + ".txt.gz", std::ios::binary);

// GZ compression:
  GZWriter outGZ;
  outGZ.SetOutputHandle(&out);

// Write stats here:

  outGZ.writeline("BAM_report\tValue");
  outGZ.writestring(myLine);
  outGZ.writeline("");

  int directionality = oJuncCount.Directional(myLine);
  outGZ.writeline("Directionality\tValue");
  outGZ.writestring(myLine);
  outGZ.writeline("");

  oFragmentsInROI.WriteOutput(myLine);
  outGZ.writeline("ROIname\ttotal_hits\tpositive_strand_hits\tnegative_strand_hits");
  outGZ.writestring(myLine);
  outGZ.writeline("");
  
  oJuncCount.WriteOutput(myLine);
  outGZ.writeline("JC_seqname\tstart\tend\tstrand\ttotal\tpos\tneg");
  outGZ.writestring(myLine);
  outGZ.writeline("");

  oSpansPoint.WriteOutput(myLine);
  outGZ.writeline("SP_seqname\tpos\ttotal\tpos\tneg");
  outGZ.writestring(myLine);
  outGZ.writeline("");
  
  oFragmentsInChr.WriteOutput(myLine);
  outGZ.writeline("ChrCoverage_seqname\ttotal\tpos\tneg");
  outGZ.writestring(myLine);
  outGZ.writeline("");
  
  oCoverageBlocks.WriteOutput(myLine, oJuncCount, oSpansPoint);
  outGZ.writestring(myLine);
  outGZ.writeline("");
  
  if (directionality != 0) {
    std::ostringstream outCoverageBlocks_dir;
    oCoverageBlocks.WriteOutput(myLine, oJuncCount, oSpansPoint, directionality); // Directional.
    outGZ.writestring(myLine);
    outGZ.writeline("");
  }
  outGZ.flush(true);
  out.flush(); out.close();
  
  // Write Coverage Binary file:
  
  std::ofstream ofCOV;
  ofCOV.open(output_file + ".cov", std::ofstream::binary);
   
  covFile outCOV;
  outCOV.SetOutputHandle(&ofCOV);
  
  oFragMap.WriteBinary(&outCOV, BB.chr_names, BB.chr_lens);
  ofCOV.close();
  
  return(0);
}

#ifdef GALAXY
// Galaxy main
int main(int argc, char * argv[]) {
	// Usage:
    // irfinder_galaxy main sample.bam IRFinder.ref.gz OutputHeader
    // irfinder_galaxy gen_map_reads genome.fa reads_to_map.fa 70 10
    // irfinder_galaxy process_mappability_bam mappedreads.bam mappability.bed

    if(strcmp(argv[1],"gen_map_reads")) {
        std::string s_genome = argv[2];
        std::string s_output = argv[3];
        int read_len = atoi(argv[4]);
        int read_stride = atoi(argv[5]);
        int read_error = atoi(argv[4]) / 2;
        IRF_GenerateMappabilityReads(s_genome, s_output, read_len, read_stride, read_error);
        exit(0);
    } else if(strcmp(argv[1],"process_mappability_bam")) {
        std::string s_bam = argv[2];
        std::string s_output = argv[3];
		int threshold = atoi(argv[4]);
        IRF_GenerateMappabilityRegions(s_bam, s_output, threshold);
        exit(0);
    } else if(strcmp(argv[1],"main")) {
        std::string s_bam = argv[2];
        std::string s_ref = argv[3];		
        std::string s_output = argv[4];		
		IRF_main(s_bam, s_ref, s_output);
	}
}
	
#endif