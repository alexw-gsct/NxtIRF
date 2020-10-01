#include "ReadBlockProcessor.h" // No need to include GZWriter.h as this is included in ReadBlockProcessor.h
#include "BAM2blocks.h"
#include "FastaReader.h"
#include "includedefine.h"

char c_complement(char n);
void reverseit(char arr[]);
std::string reverse_complement(std::string sequence);

std::string GenerateReadError(char * input_read, unsigned int read_len, unsigned int error_pos,
                              unsigned int direction, unsigned int error_seed);

bool checkDNA(const std::string& strand);

int IRF_GenerateMappabilityReads(std::string genome_file, std::string out_fa,
	int read_len, int read_stride, int error_pos);

int IRF_GenerateMappabilityRegions(std::string bam_file, std::string output_file, int threshold);