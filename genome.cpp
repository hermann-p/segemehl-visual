#include "genome.h"
#include "readcontainer.h"
#include "utils.h"

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <fstream>
#include <cstdlib>

#include <sstream>
#include <cstdio>
#include <iostream>

chr_num_t Genome::getChrNum ( const std::string name ) const {
  auto lookup = chrNums.find(name);
  assume(lookup != chrNums.end(), "Chromosome " + name + " not found");
  return lookup->second;
}


std::string Genome::getChrName ( const chr_num_t num ) const {
  assume(num < chrNames.size() && num >= 0, "Index out of range");
  return chrNames[num];
}


void Genome::createChromosome ( const std::string name, const uint32_t length ) {
  chr_num_t N = chrNames.size();
  chrNames.push_back(name);
  chrLens.push_back(length);
  chrNums.emplace( name, N );
  readPos.push_back(new std::map<chr_pos_t, ReadContainer*>());
}


ReadContainer* Genome::createRead( const chr_num_t chr, const chr_pos_t pos, const short len ) {
  ReadContainer* rc = getReadAt(chr, pos);
  if (rc == nullptr) {                      // no read on chr:pos exists yet
    rc = new ReadContainer(chr, pos, len);
    readPos[chr]->emplace(pos, rc);          // register read with genome navigation map
  }
  return rc;
}


ReadContainer* Genome::getReadAt ( const chr_num_t chr, const chr_pos_t pos ) const {
  auto chromosome = readPos[chr];
  auto lookup = chromosome->find(pos);
  if (lookup == chromosome->end()) {
    return nullptr;
  }
  return lookup->second;
}


ReadContainer* Genome::getReadAt ( const std::string chr, const chr_pos_t pos ) const {
  return getReadAt( getChrNum(chr), pos );
}


void Genome::read ( std::ifstream& input ) {
  log("Genome: reading file...");
  unsigned int L(1);
  for (std::string line; getline(input, line);) {
    if (line[0] == '@') {
      assume(parseHeaderLine(line), "Could not parse header in line " + std::to_string(L), false);
    }
    else {
      assume(parseDataLine(line), "Could not parse data in line " + std::to_string(L), false);
    }
    L++;
  }
}


void Genome::read ( std::string& fileName ) {
  std::ifstream input(fileName);
  read(input);
}


std::ifstream& operator>> ( std::ifstream& input, Genome& target ) {
  target.read(input);
  return input;
}


bool Genome::parseHeaderLine ( const std::string& line ) {
  try {
    switch (line[1]) {
      case 'H': // Token @HD: General file info
	log(line);
	break;
      case 'S': { // Token @SQ: Reference dictionary entry
	  auto tokens = strsplit(line, "\t", false);
	  std::string name = strsplit(tokens[1], ":")[1];
	  unsigned short len = atoi( strsplit(tokens[2], ":")[2].c_str() );
	  createChromosome(name, len);
	}
	break;
      case 'P': // Token @PG: program used
	assume(strsplit(line, "\t")[1] == "ID:segemehl", "This program is designed for segemehl output only!", false);
      default:  // Ignore other tokens
	break;
    }
  } catch (const std::exception e) { // keep going on errors, but warn
    return false;
  }
  return true;
}


bool Genome::parseDataLine ( const std::string& line ) {
  try {
    auto tokens = strsplit(line, "\t");
    // int flag = atoi( tokens[FLAG].c_str() );
    chr_num_t chr = getChrNum(tokens[RNAME]);
    chr_pos_t pos5 = atoi( tokens[POS].c_str() );
    chr_pos_t pos3 = pos5 + atoi( tokens[TLEN].c_str() );
    bool hasNext = false;
    chr_num_t nextChr = 0;
    chr_pos_t nextPos = 0;
    bool hasPrev = false;
    chr_num_t prevChr = 0;
    chr_pos_t prevPos = 0;
    int start = 0;
    int end = 0;

    auto iter = tokens.begin();
    for (int i(0); i < QUAL; i++) iter++;
    for (;iter != tokens.end(); iter++) {
      if ((*iter)[0] == 'X') {
	switch((*iter)[1]) {
	  case 'X': // start of current split
	    start = atoi( (*iter).substr(5).c_str() );
	    break;
	  case 'Y': // enf of current split
	    end = atoi( (*iter).substr(5).c_str() );
	    break;
	  case 'Q': // number of current split
	    break;
	  
	  case 'P': // refseq of prev split
	    hasPrev = true;
	    prevChr = (chr_num_t)getChrNum( (*iter).substr(5) );
	    break;
	  case 'U': // 3' of prev
	    prevPos = (chr_pos_t)atoi( (*iter).substr(5).c_str() );
	    break;
	    
	  case 'C': // refseq of next split
	    hasNext = true;
	    nextChr = (chr_num_t)getChrNum( (*iter).substr(5) );
	    break;
	    
	  case 'V': // 5' of next split
	    nextPos = (chr_pos_t)atoi( (*iter).substr(5).c_str() );
	    break;
	  default:
	    break;
	}
      }
    }
    
    if (pos3 == pos5 && end - start > 0) {     // Try to fix length (information often missing or wrong in .sam)
      pos3 = pos5 + (end-start);
    }
    if (pos3 == pos5) {                        // Try to calc length from cigar-string (slow but should always work)
      pos3 = pos5 + calcLength(tokens[CIGAR]);
    }
    ReadContainer* rc = getReadAt(chr, pos5);
    if (rc == nullptr) {
      rc = new ReadContainer(chr, pos5, pos3-pos5);
      registerRead(rc);
    }
    
    std::cout << (int)chr << "@" << pos5 << "-" << pos3 << "; next: " << (int)nextChr << "@" << nextPos << "; prev: " << (int)prevChr << "@" << prevPos << std::endl;
      
    if (hasNext) {
      rc->addDownstreamRead(*this, nextChr, nextPos);
    }
    if (hasPrev) {
      rc->addUpstreamRead(*this, prevChr, -prevPos);
    }
    std::cout << std::endl;
  } catch (const std::exception e) {
    return false;
  }
  return true;
}

void Genome::registerRead (  ReadContainer* rc ) {
  const chr_num_t chr = rc->chromosome;
  readPos.at(chr)->emplace(rc->fivePrimeEnd, rc);   // Link 5' end to identify as successor
  readPos.at(chr)->emplace(-rc->threePrimeEnd, rc); // Link 3' end to identify as predecessor
}


void Genome::printout( ReadContainer* element ) {
  static uint runID = 1;
  assume(readPos.size() > 0, "Nothing to print...", false);
  if (element == nullptr) {
    int i(0);
    for (auto& chr: readPos) {
      std::cout << "Chromosome " << getChrName(i++) << " has " << chr->size() << " elements" << std::endl;
      for (auto& read: *chr) {
	log("  Next read:");
	  printout(read.second);
	  runID++;
      }
    }
  }
  else {
//    std::cout << "    " << element->fivePrimeRead.size() << " elements upstream\n";
    element->runID = runID; // mark as processed

    for (auto& iter : element->fivePrimeRead) {
      if (iter->runID != runID) printout(iter);
    }
    std::stringstream str;
    str << "    Chr: " << (int)element->chromosome << ", 5': " << (int)element->fivePrimeEnd << ", 3': " << (int)element->threePrimeEnd << std::endl;
    std::string line;
    getline(str, line);
    log(line);
//    std::cout << "    " << element->threePrimeRead.size() << " elements downstream\n";
    for (auto& iter : element->threePrimeRead) {
      if (iter->runID != runID) printout(iter);
    }
  }
}


uint Genome::calcLength ( const std::string cigar ) const {
  int L = 0;
  std::string tmp = "";
  for (size_t i(0); i < cigar.length(); i++) {
    if (cigar[i] == 'M' || cigar[i] == 'I') {
      L += atoi(tmp.c_str());
      tmp = "";
    }
    else if (cigar[i] == 'D') {
      L -= atoi(tmp.c_str());
      tmp = "";
    }
    else {
      tmp += cigar[i];
    }
  }
  return (uint)L;
}


