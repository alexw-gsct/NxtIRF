#include "ReadBlockProcessor.h"
#include "includedefine.h"
// using namespace std;

//chrName_junc_count holds the data structure -- ChrName(string) -> Junc Start/End -> count.
//chrID_junc_count holds the ChrID -> ...
//  where the ChrID is the ChrID relating to the appropriate ChrName, as understood by the currently processed BAM file.
void JunctionCount::ChrMapUpdate(const std::vector<std::string> &chrmap) {
	chrID_junc_count.resize(0);
	chrID_juncLeft_count.resize(0);
	chrID_juncRight_count.resize(0);
	// Below could be done with an iterator - i is not used except for element access of the single collection.
	for (unsigned int i = 0; i < chrmap.size(); i++) {
		chrID_junc_count.push_back( &chrName_junc_count[chrmap.at(i)] );
		chrID_juncLeft_count.push_back( &chrName_juncLeft_count[chrmap.at(i)] );
		chrID_juncRight_count.push_back( &chrName_juncRight_count[chrmap.at(i)] );
	}
	
	
//		std::map<string, std::map<unsigned int,unsigned int[2]>> chrName_juncLeft_count;
//		std::vector<std::map<unsigned int,unsigned int[2]>*> chrID_juncLeft_count;

//		std::map<string, std::map<unsigned int,unsigned int[2]>> chrName_juncRight_count;
//		std::vector<std::map<unsigned int,unsigned int[2]>*> chrID_juncRight_count;

}

void JunctionCount::loadRef(std::istringstream &IN) {
	// ChrName, Start, End, direction(+/-/.).
	std::string myLine;
	std::string myField;
	myLine.reserve(1000);
	myField.reserve(100);
	unsigned int start;
	unsigned int end;
	string s_chr;
	s_chr.reserve(30);
	string direction;

	while(!IN.eof() && !IN.fail()) {
		getline(IN, myLine, '\n');

		if (IN.eof() || IN.fail()) {
			if (myLine.length() == 0) {
				// This line is empty - just a blank line at the end of the file.
				// Checking at this stage allows correct handling of files both with and without a trailing \n after the last record.
				break;
			}else{
				// Error line in input, ignore.
				break;
			}
		}

		std::istringstream lineStream;
		lineStream.str(myLine);
		
		getline(lineStream, s_chr, '\t');
		getline(lineStream, myField, '\t');
		start = stol(myField);
		getline(lineStream, myField, '\t');
		end = stol(myField);
		getline(lineStream, direction, '\t');
		
		if (direction == "-")  {
			chrName_junc_count[s_chr][make_pair(start,end)][2] = 1;
		}else if (direction == "+") {
			chrName_junc_count[s_chr][make_pair(start,end)][2] = 2;
		}
	}
}

void JunctionCount::ProcessBlocks(const FragmentBlocks &blocks) {
	for (int index = 0; index < blocks.readCount; index ++) {
		//Walk each *pair* of blocks. ie: ignore a read that is just a single block.
		for (unsigned int j = 1; j < blocks.rLens[index].size(); j++) {
			if ((blocks.rLens[index][j-1] >= 5) && (blocks.rLens[index][j] >= 5)) {
				(*chrID_junc_count[blocks.chr_id])[
					make_pair(
						blocks.readStart[index] + blocks.rStarts[index][j-1] + blocks.rLens[index][j-1],
						blocks.readStart[index] + blocks.rStarts[index][j])
					][blocks.direction]++;
				(*chrID_juncLeft_count[blocks.chr_id])[
						blocks.readStart[index] + blocks.rStarts[index][j-1] + blocks.rLens[index][j-1]
					][blocks.direction]++;
				(*chrID_juncRight_count[blocks.chr_id])[
						blocks.readStart[index] + blocks.rStarts[index][j]
					][blocks.direction]++;
			}
		}
	}
}

int JunctionCount::WriteOutput(std::string& output, std::string& QC) const {
  std::ostringstream oss; std::ostringstream oss_qc; 
	int junc_anno = 0;
	int junc_unanno = 0;
	for (auto itChr=chrName_junc_count.begin(); itChr!=chrName_junc_count.end(); itChr++) {
		string chr = itChr->first;
		for (auto itJuncs=itChr->second.begin(); itJuncs!=itChr->second.end(); ++itJuncs) {
			if((itJuncs->second)[2] != 0) {
				junc_anno += ((itJuncs->second)[1] + (itJuncs->second)[0]);
			} else {
				junc_unanno += ((itJuncs->second)[1] + (itJuncs->second)[0]);
			}
			oss << chr << "\t" << itJuncs->first.first << "\t" << itJuncs->first.second
				<< "\t" << ( (itJuncs->second)[2] == 1 ? "-" : (itJuncs->second)[2] == 2 ? "+" : "." )
				<< "\t" << ((itJuncs->second)[1] + (itJuncs->second)[0])
				<< "\t" << (itJuncs->second)[1]
				<< "\t" << (itJuncs->second)[0] << "\n";
		}
	}
	oss_qc << "Annotated Junctions" << "\t" << junc_anno << "\n"
		<< "Unannotated Junctions" << "\t" << junc_unanno << "\n";
	
  output = oss.str();
	QC.append(oss_qc.str());
	return 0;
}

int JunctionCount::Directional(std::string& output) const {
	unsigned int dir_same = 0;
	unsigned int dir_diff = 0;

	unsigned int dir_evidence = 0;
	unsigned int nondir_evidence = 0;
	unsigned int dir_evidence_known = 0;
	unsigned int nondir_evidence_known = 0;

    std::ostringstream oss;    

	for (auto itChr=chrName_junc_count.begin(); itChr!=chrName_junc_count.end(); itChr++) {
		for (auto itJuncs=itChr->second.begin(); itJuncs!=itChr->second.end(); ++itJuncs) {
			if (((itJuncs->second)[1] + (itJuncs->second)[0]) > 8) {
				if ((itJuncs->second)[0] > (itJuncs->second)[1] * 4) {
					dir_evidence++;
					if ((itJuncs->second)[2] == 1) { //Ref is "-"
						dir_same++;
					}else if ((itJuncs->second)[2] == 2) {
						dir_diff++;
					}
				}else if ((itJuncs->second)[1] > (itJuncs->second)[0] * 4) {
					dir_evidence++;
					if ((itJuncs->second)[2] == 2) { //Ref is "+"
						dir_same++;
					}else if ((itJuncs->second)[2] == 1) {
						dir_diff++;
					}				
				}else{
					nondir_evidence++;
					if ((itJuncs->second)[2] > 0) {
						nondir_evidence_known++;
					}
				}
			}
		}
	}
	dir_evidence_known = dir_same + dir_diff;
	oss << "Dir evidence\t" << dir_evidence << "\n";
	oss << "Nondir evidence\t" << nondir_evidence << "\n";
	oss << "Dir evidence known junctions\t" << dir_evidence_known << "\n";
	oss << "Nondir evidence known junctions\t" << nondir_evidence_known << "\n";
	oss << "Dir matches ref\t" << dir_same << "\n";
	oss << "Dir opposed to ref\t" << dir_diff << "\n";
	oss << "Dir score all (0-10000)\t" << ((long long)dir_evidence * 10000 / (dir_evidence + nondir_evidence + 1)) << "\n"; //+1 to prevent divide by zero errors.
	long dir_score_known = ((long long)dir_evidence_known * 10000 / (dir_evidence_known + nondir_evidence_known + 1));
	oss << "Dir score known junctions (0-10000)\t" << dir_score_known << "\n";

	if ((dir_same > dir_diff * 100) && (dir_score_known >= 9000)) {
        oss << "Overall Directionality\t" << 1 << '\n';
        output = oss.str();
		return 1;
	}else if ((dir_diff > dir_same * 100) && (dir_score_known >= 9000)) {
        oss << "Overall Directionality\t" << -1 << '\n';
        output = oss.str();
		return -1;
	}else{
        oss << "Overall Directionality\t" << 0 << '\n';
        output = oss.str();
		return 0;
	}
}

unsigned int JunctionCount::lookup(std::string ChrName, unsigned int left, unsigned int right, bool direction) const {
	try {
		return chrName_junc_count.at(ChrName).at(make_pair(left, right))[direction];
	}catch (const std::out_of_range& e) {
	}
	return 0;
}
unsigned int JunctionCount::lookup(std::string ChrName, unsigned int left, unsigned int right) const {
	try {
		return chrName_junc_count.at(ChrName).at(make_pair(left, right))[0] + chrName_junc_count.at(ChrName).at(make_pair(left, right))[1];
	}catch (const std::out_of_range& e) {
	}
	return 0;
}
unsigned int JunctionCount::lookupLeft(std::string ChrName, unsigned int left, bool direction) const {
	try {
		return chrName_juncLeft_count.at(ChrName).at(left)[direction];
	}catch (const std::out_of_range& e) {
	}
	return 0;
}
unsigned int JunctionCount::lookupLeft(std::string ChrName, unsigned int left) const {
	try {
		return chrName_juncLeft_count.at(ChrName).at(left)[0] + chrName_juncLeft_count.at(ChrName).at(left)[1];
	}catch (const std::out_of_range& e) {
	}
	return 0;
}
unsigned int JunctionCount::lookupRight(std::string ChrName, unsigned int right, bool direction) const {
	try {
		return chrName_juncRight_count.at(ChrName).at(right)[direction];
	}catch (const std::out_of_range& e) {
	}
	return 0;
}
unsigned int JunctionCount::lookupRight(std::string ChrName, unsigned int right) const {
	try {
		return chrName_juncRight_count.at(ChrName).at(right)[0] + chrName_juncRight_count.at(ChrName).at(right)[1];
	}catch (const std::out_of_range& e) {
	}
	return 0;
}





int SpansPoint::WriteOutput(std::string& output, std::string& QC) const {
  std::ostringstream oss; std::ostringstream oss_qc; 
	int spans_reads = 0;	
	for (auto itChrPos=chrName_pos.begin(); itChrPos!=chrName_pos.end(); itChrPos++) {
		string chr = itChrPos->first;

		auto itCountPos=chrName_count[1].at(chr).begin();
		auto itCountNeg=chrName_count[0].at(chr).begin();

		
		for (auto itPosition=itChrPos->second.begin(); itPosition!=itChrPos->second.end(); ++itPosition) {
			spans_reads += (*itCountPos + *itCountNeg);
			oss << chr << "\t" << *itPosition << "\t" << (*itCountPos + *itCountNeg) << "\t" << *itCountPos << "\t" << *itCountNeg << "\n";
			++itCountPos;
			++itCountNeg;
		}
	}
	oss_qc << "Spans Reads\t" << spans_reads << "\n";
    output = oss.str();
		QC.append(oss_qc.str());
	return 0;
}


unsigned int SpansPoint::lookup(std::string chrName, unsigned int pos, bool direction) const {
	auto it_pos = std::lower_bound(chrName_pos.at(chrName).begin(), chrName_pos.at(chrName).end(), pos);
	if (it_pos == chrName_pos.at(chrName).end() || *it_pos != pos) {
		// throw not-found/out-of-bounds exception?
		throw std::out_of_range("Pos not found - SpansPoint::lookup");
		return 0;
	}else{
		// Then use that offset into the other vectors.
		return chrName_count[direction].at(chrName).at(it_pos - chrName_pos.at(chrName).begin());
	}
}

unsigned int SpansPoint::lookup(std::string chrName, unsigned int pos) const {
	//	std::map<string, std::vector<int>> chrName_pos;
	auto it_pos = std::lower_bound(chrName_pos.at(chrName).begin(), chrName_pos.at(chrName).end(), pos);
	if (it_pos == chrName_pos.at(chrName).end() || *it_pos != pos) {
		// throw not-found/out-of-bounds exception?
		throw std::out_of_range("Pos not found - SpansPoint::lookup");
		return 0;
	}else{
		// Then use that offset into the other vectors.
		return (
			chrName_count[0].at(chrName).at(it_pos - chrName_pos.at(chrName).begin())
			+ chrName_count[1].at(chrName).at(it_pos - chrName_pos.at(chrName).begin())
			);
	}
}


void SpansPoint::setSpanLength(unsigned int overhang_left, unsigned int overhang_right) {
	overhangLeft = overhang_left;
	overhangRight = overhang_right;
	overhangTotal = overhang_right + overhang_left;
}

void SpansPoint::ProcessBlocks(const FragmentBlocks &blocks) {
	std::vector<unsigned int>::iterator it_position;

	//Walk each read within the fragment (1 or 2).
	for (int index = 0; index < blocks.readCount; index ++) {
		//Walk each block within each read.
		for (unsigned int j = 0; j < blocks.rLens[index].size(); j++) {
			if ( blocks.rLens[index][j] > overhangTotal ) {
				//Block is long enough it may sufficiently overhang a point of interest.
				it_position = std::upper_bound(
						(*chrID_pos.at(blocks.chr_id)).begin(),
						(*chrID_pos.at(blocks.chr_id)).end(),
						blocks.readStart[index] + blocks.rStarts[index][j] + overhangLeft - 1
				);  // -1 --- as the test is > rather than >=.
				while (it_position != (*chrID_pos.at(blocks.chr_id)).end() && *it_position < (blocks.readStart[index] + blocks.rStarts[index][j] + blocks.rLens[index][j] )) {
					//increment corresponding counter.
					(*chrID_count[blocks.direction].at(blocks.chr_id)).at(it_position - (*chrID_pos.at(blocks.chr_id)).begin())++;
					it_position++;
				}
			}
		}
	}
}


void SpansPoint::loadRef(std::istringstream &IN) {
	// TODO: will we ever want to store some additional info -- eg: String name of each position? Not right now.
	std::string myLine;
	std::string myField;
	myLine.reserve(1000);
	myField.reserve(100);
	int pos;
	string s_chr;
	s_chr.reserve(30);
	string direction;

	while ( !IN.eof() && !IN.fail() ) {
		getline(IN, myLine, '\n');

		if (IN.eof() || IN.fail()) {
			if (myLine.length() == 0) {
				// This line is empty - just a blank line at the end of the file.
				// Checking at this stage allows correct handling of files both with and without a trailing \n after the last record.
				break;
			}else{
				// Error line in input, ignore.
				break;
			}
		}

		std::istringstream lineStream;
		lineStream.str(myLine);
		
		getline(lineStream, s_chr, '\t');
		getline(lineStream, myField, '\t');
		pos = stol(myField);
		
		getline(lineStream, direction, '\t');
		
		chrName_pos[s_chr].push_back(pos);
	}
	
	for (std::map<string, std::vector<unsigned int>>::iterator it_chr=chrName_pos.begin(); it_chr!=chrName_pos.end(); it_chr++) {	
		std::sort( it_chr->second.begin(), it_chr->second.end() );
		// We now have chrName_pos sorted by position.
		chrName_count[0][it_chr->first].resize(it_chr->second.size(), 0);
		chrName_count[1][it_chr->first].resize(it_chr->second.size(), 0);
		// Just created a vector of the same size as the position one to store the counter.
		// initialise the count vector to the same length, with zero start.
	}
}

void SpansPoint::ChrMapUpdate(const std::vector<std::string> &chrmap) {
	chrID_pos.resize(0);
	chrID_count[0].resize(0);
	chrID_count[1].resize(0);
	for (unsigned int i = 0; i < chrmap.size(); i++) {
		chrID_pos.push_back( &chrName_pos[chrmap.at(i)] );
		chrID_count[0].push_back( &chrName_count[0][chrmap.at(i)] );
		chrID_count[1].push_back( &chrName_count[1][chrmap.at(i)] );
	}
}


void FragmentsInROI::ChrMapUpdate(const std::vector<std::string> &chrmap) {
	chrID_ROI.resize(0);
	chrID_count[0].resize(0);
	chrID_count[1].resize(0);
	for (unsigned int i = 0; i < chrmap.size(); i++) {
		chrID_ROI.push_back( &chrName_ROI[chrmap.at(i)] );
		chrID_count[0].push_back( &chrName_count[0][chrmap.at(i)] );
		chrID_count[1].push_back( &chrName_count[1][chrmap.at(i)] );
	}
}

int FragmentsInROI::WriteOutput(std::string& output, std::string& QC) const {
	std::ostringstream oss; std::ostringstream oss_QC;
	int count_Intergenic = 0; int count_rRNA = 0; int count_NonPolyA = 0;
	for (std::map<string, unsigned long>::const_iterator itID=RegionID_counter[1].begin(); 
			itID!=RegionID_counter[1].end(); ++itID) {
		std::istringstream iss;
		iss.str(itID->first);
		string ROI_type;
		getline(iss, ROI_type, '/');
		if(ROI_type.compare(0, 10, "Intergenic") == 0) {
			count_Intergenic += (itID->second + RegionID_counter[0].at(itID->first));
		} else if(ROI_type.compare(0, 4, "rRNA") == 0) {
			count_rRNA += (itID->second + RegionID_counter[0].at(itID->first));
		} else if(ROI_type.compare(0, 8, "NonPolyA") == 0) {
			count_NonPolyA += (itID->second + RegionID_counter[0].at(itID->first));	
		}
		oss << itID->first << "\t" << (itID->second + RegionID_counter[0].at(itID->first)) << "\t" 
			<< itID->second << "\t" << RegionID_counter[0].at(itID->first) << "\n";
		//Outputs tab separated: ROIname, total hits, positive-strand hits, negative-strand hits.
	}
    output = oss.str();
		
		oss_QC << "Intergenic Reads\t" << count_Intergenic << "\n"
			<< "rRNA Reads\t" << count_rRNA << "\n"
			<< "NonPolyA Reads\t" << count_NonPolyA << "\n";
		QC.append(oss_QC.str());
	return 0;
}

void FragmentsInROI::loadRef(std::istringstream &IN) {

	std::string myLine;
	std::string myField;
	myLine.reserve(1000);
	myField.reserve(100);
	int start;
	int end;
	string s_chr;
	s_chr.reserve(30);
	string s_name;
	s_name.reserve(200);

	while ( !IN.eof() && !IN.fail() ) {
		// Input ref:  chr - start - end - name - dir(?). (name\tdir could be considered a single variable)
		getline(IN, myLine, '\n');

		if (IN.eof() || IN.fail()) {
			if (myLine.length() == 0) {
				// This line is empty - just a blank line at the end of the file.
				// Checking at this stage allows correct handling of files both with and without a trailing \n after the last record.
				break;
			}else{
				// Error line in input, ignore.
				break;
			}
		}

		std::istringstream lineStream;
		lineStream.str(myLine);
		
		getline(lineStream, s_chr, '\t');
		getline(lineStream, myField, '\t');
		start = stol(myField);
		getline(lineStream, myField, '\t');
		end = stol(myField);

		getline(lineStream, s_name, '\t');

		chrName_ROI[s_chr].push_back(std::make_pair(end, start));
		chrName_count[0][s_chr].push_back(&RegionID_counter[0][s_name]);
		chrName_count[1][s_chr].push_back(&RegionID_counter[1][s_name]);

	}
}

void FragmentsInROI::ProcessBlocks(const FragmentBlocks &blocks) {
	std::vector<std::pair<unsigned int,unsigned int>>::iterator it_ROI;

	unsigned int frag_start = blocks.readStart[0];
	unsigned int frag_end = blocks.readEnd[0];
	if (blocks.readCount > 1 && blocks.readEnd[1] > frag_end) {
		frag_end = blocks.readEnd[1];
	}
	
	// Frag start, Frag end.
	// See if this is fully inside one of the ref-regions.

	it_ROI = std::lower_bound(
			(*chrID_ROI.at(blocks.chr_id)).begin(),
			(*chrID_ROI.at(blocks.chr_id)).end(),
			std::make_pair(frag_end, frag_end)
	);
	
	if (it_ROI != (*chrID_ROI.at(blocks.chr_id)).end() ) {
		if (frag_start >= it_ROI->second && frag_end <= it_ROI->first) {
			(*(*chrID_count[blocks.direction].at(blocks.chr_id)).at(it_ROI - (*chrID_ROI.at(blocks.chr_id)).begin()))++;			
		}
	}
}



void FragmentsInChr::ProcessBlocks(const FragmentBlocks &blocks) {
	(*chrID_count.at(blocks.chr_id))[blocks.direction]++;
}

void FragmentsInChr::ChrMapUpdate(const std::vector<string> &chrmap) {
	chrID_count.resize(0);
	for (unsigned int i = 0; i < chrmap.size(); i++) {
		chrName_count[chrmap.at(i)].resize(2); // This data structure isn't auto initializing - unlike all the other structures. Or maybe just a vector can't access via [] until a position exists? But a map is fine. Makes sense.
		chrID_count.push_back( &chrName_count[chrmap.at(i)] );
	}
}

int FragmentsInChr::WriteOutput(std::string& output, std::string& QC) const {
	std::ostringstream oss; std::ostringstream oss_QC;
	int count_Mito = 0; int count_ERCC = 0;
    for (auto itChr=chrName_count.begin(); itChr!=chrName_count.end(); itChr++) {
		string chr = itChr->first;
		if(chr.compare(0, 1, "M") == 0 || chr.compare(0, 2, "MT") == 0) {
			count_Mito += ((itChr->second)[1] + (itChr->second)[0]);
		} else if(chr.compare(0, 4, "ERCC") == 0) {
			count_ERCC += ((itChr->second)[1] + (itChr->second)[0]);
		}		
		oss << itChr->first << "\t"
			<< ((itChr->second)[1] + (itChr->second)[0]) << "\t"
			<< (itChr->second)[1] << "\t"
			<< (itChr->second)[0] << "\n";
	}
    output = oss.str();
		oss_QC << "Mitochondrial Reads\t" << count_Mito << "\n"
			<< "ERCC Reads\t" << count_ERCC << "\n";

		QC.append(oss_QC.str());
	return 0;
}

void FragmentsMap::ChrMapUpdate(const std::vector<string> &chrmap) {
  chrID_count[0].resize(0);
  chrID_count[1].resize(0);
  chrID_count[2].resize(0);
  for (unsigned int i = 0; i < chrmap.size(); i++) {
    chrName_count[0][chrmap.at(i)].insert({0,0}); // Insert dummy pair
    chrID_count[0].push_back( &chrName_count[0][chrmap.at(i)] );
    chrName_count[1][chrmap.at(i)].insert({0,0}); // Insert dummy pair
    chrID_count[1].push_back( &chrName_count[1][chrmap.at(i)] );
    chrName_count[2][chrmap.at(i)].insert({0,0}); // Insert dummy pair
    chrID_count[2].push_back( &chrName_count[2][chrmap.at(i)] );
  }
}

void FragmentsMap::ProcessBlocks(const FragmentBlocks &blocks) {
  std::map<unsigned int, int>::iterator it_position;
  // Contains chr, then pos, then incremental coverage (+ indicates increased coverage compared to prior loci)
  
  //Walk each read within the fragment (1 or 2).
  for (int index = 0; index < blocks.readCount; index ++) {
    //Walk each block within each read.
    for (unsigned int j = 0; j < blocks.rLens[index].size(); j++) {
      // Stranded 
      it_position = (*chrID_count[blocks.direction].at(blocks.chr_id)).find(blocks.readStart[index] + blocks.rStarts[index][j]);
      if (it_position == (*chrID_count[blocks.direction].at(blocks.chr_id)).end()) {
        (*chrID_count[blocks.direction].at(blocks.chr_id)).insert({ blocks.readStart[index] + blocks.rStarts[index][j], 1});
      } else {
        it_position->second += 1;
        if(it_position->second == 0) {
          (*chrID_count[blocks.direction].at(blocks.chr_id)).erase(it_position);
        }
      }
      it_position = (*chrID_count[blocks.direction].at(blocks.chr_id)).find(blocks.readStart[index] + blocks.rStarts[index][j] + blocks.rLens[index][j]);
      if (it_position == (*chrID_count[blocks.direction].at(blocks.chr_id)).end()) {
        (*chrID_count[blocks.direction].at(blocks.chr_id)).insert({ blocks.readStart[index] + blocks.rStarts[index][j] + blocks.rLens[index][j], -1});
      } else {
        it_position->second -= 1;
        if(it_position->second == 0) {
          (*chrID_count[blocks.direction].at(blocks.chr_id)).erase(it_position);
        }
      }
      // Unstranded 
      it_position = (*chrID_count[2].at(blocks.chr_id)).find(blocks.readStart[index] + blocks.rStarts[index][j]);
      if (it_position == (*chrID_count[2].at(blocks.chr_id)).end()) {
        (*chrID_count[2].at(blocks.chr_id)).insert({ blocks.readStart[index] + blocks.rStarts[index][j], 1});
      } else {
        it_position->second += 1;
        if(it_position->second == 0) {
          (*chrID_count[2].at(blocks.chr_id)).erase(it_position);
        }
      }
      it_position = (*chrID_count[2].at(blocks.chr_id)).find(blocks.readStart[index] + blocks.rStarts[index][j] + blocks.rLens[index][j]);
      if (it_position == (*chrID_count[2].at(blocks.chr_id)).end()) {
        (*chrID_count[2].at(blocks.chr_id)).insert({ blocks.readStart[index] + blocks.rStarts[index][j] + blocks.rLens[index][j], -1});
      } else {
        it_position->second -= 1;
        if(it_position->second == 0) {
          (*chrID_count[2].at(blocks.chr_id)).erase(it_position);
        }
      }
    }
  }
}

int FragmentsMap::WriteBinary(covFile *os, const std::vector<std::string> chr_names, const std::vector<int32_t> chr_lens) const {
  // Write COV file as binary

  // Issue is map constructs auto-sort
  // Need to put chrs and lengths into a map structure
  std::map< std::string, int32_t > chrmap;
  
  // Arrange chromosomes in same order as arranged by mapping chrs
  for(unsigned int i = 0; i < chr_names.size(); i++) {
      chrmap.insert({chr_names[i], chr_lens[i]});
  }
  // Re-push into alphabetical ordered chromosomes
  std::vector<std::string> sort_chr_names;
  std::vector<int32_t> sort_chr_lens;
  for (auto chr = chrmap.begin(); chr != chrmap.end(); chr++) {
    sort_chr_names.push_back(chr->first);
    sort_chr_lens.push_back(chr->second);
  }
  os->WriteHeader(sort_chr_names, sort_chr_lens);

  unsigned int refID = 0;
  for(unsigned int j = 0; j < 3; j++) {
    for (auto itChr=chrName_count[j].begin(); itChr!=chrName_count[j].end(); itChr++) {
      unsigned int coordpos = 0;
      unsigned int coorddepth = 0;
      bool writefirst = true;
            
      for(auto it_pos = itChr->second.begin(); it_pos != itChr->second.end(); it_pos++) {
        if(writefirst) {
          writefirst = false;
          if(it_pos->first == 0) {
            // Write coverage only
            coorddepth += it_pos->second;
          } else {
            coorddepth = 0;

            os->WriteEntry(refID, coorddepth, it_pos->first);
            
            coorddepth += it_pos->second;
            coordpos = it_pos->first;
          }
        } else {
          // coorddepth should already be recorded previously
          os->WriteEntry(refID, coorddepth, it_pos->first - coordpos);
          
          coorddepth += it_pos->second;
          
          coordpos = it_pos->first;
        }
      }
      // Write last entry for remainder of chromosome length
      os->WriteEntry(refID, coorddepth, chrmap[itChr->first] - coordpos);
      
      refID += 1; 
    }
  }
  os->FlushBody();
  return(0);
}

int FragmentsMap::WriteOutput(std::ostream *os, 
    const std::vector<std::string> chr_names, const std::vector<int32_t> chr_lens, 
    int threshold) const {
    // This is called on mappability
  // Issue is map constructs auto-sort
  // Need to put chrs and lengths into a map structure
  std::map< std::string, int32_t > chrmap;
  
  // Arrange chromosomes in same order as arranged by mapping chrs
  for(unsigned int i = 0; i < chr_names.size(); i++) {
      chrmap.insert({chr_names[i], chr_lens[i]});
  }
  // Re-push into alphabetical ordered chromosomes
  std::vector<std::string> sort_chr_names;
  std::vector<int32_t> sort_chr_lens;
  for (auto chr = chrmap.begin(); chr != chrmap.end(); chr++) {
    sort_chr_names.push_back(chr->first);
    sort_chr_lens.push_back(chr->second);
  }    
    
  for (auto itChr=chrName_count[2].begin(); itChr!=chrName_count[2].end(); itChr++) {
    int coverage = 0;
    bool covered = false;
    
    if (itChr->second.begin()->first == 0 && itChr->second.begin()->second > threshold) {
      covered = true;
    } else {
		// Write first coordinate
	  *os << itChr->first << "\t0\t";
	}
    for(auto it_pos = itChr->second.begin(); it_pos != itChr->second.end(); it_pos++) {
      coverage += it_pos->second;
      if(coverage > threshold) {
        if(covered) {
          // do nothing
        } else {
          *os << it_pos->first << '\n';
          covered = true;
        }
      } else {
        if(covered) {
          *os << itChr->first << "\t"
              << it_pos->first << "\t";
          covered = false;
        } else {
          
          // do nothing
        }
      }
    }
    // Write last entry
    if(!covered) {
      *os << chrmap[itChr->first] << "\n";    
    }
  }
  return 0;
}


JunctionCount::~JunctionCount() {
    chrName_junc_count.clear();
    chrName_juncLeft_count.clear();
    chrName_juncRight_count.clear();
}

SpansPoint::~SpansPoint() {
    chrName_pos.clear();
    chrName_count[0].clear();
    chrName_count[1].clear();
}

FragmentsInChr::~FragmentsInChr() {
    chrName_count.clear();
}

FragmentsMap::~FragmentsMap() {
  chrName_count[0].clear();
  chrName_count[1].clear();
  chrName_count[2].clear();
}


FragmentsInROI::~FragmentsInROI() {
    RegionID_counter[0].clear();
    RegionID_counter[1].clear();
    chrName_ROI.clear();
    chrName_count[0].clear();
    chrName_count[1].clear();
    chrName_ROI_text.clear();
}