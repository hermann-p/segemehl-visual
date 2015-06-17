#include "genome.h"
#include "readcontainer.h"
#include "utils.h"

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <fstream>
#include <cstdlib>
#include <memory>

#include <sstream>
#include <cstdio>
#include <iostream>

Genome::Genome ( bool init_m, bool init_c ) 
  : readPos(new chr_map_t),
    multistrand((init_m) ? new link_set_t : nullptr),
    circular((init_c) ? new link_set_t : nullptr)
{
}

Genome::~Genome() {
  for (auto& chromosome : *readPos) {
    delete chromosome;
  }
  if (multistrand != nullptr) {
    delete multistrand;
  }
  if (circular != nullptr) {
    delete circular;
  }
}


chr_num_t Genome::getChrNum ( const std::string name ) const {
  auto lookup = chrNums.find(name);
  assume(lookup != chrNums.end(), "Genome::getChrNum: Chromosome " + name + " not found");
  return lookup->second;
}


std::string Genome::getChrName ( const chr_num_t num ) const {
  assume(!chrNames.empty(), "Genome::getChrName: name list empty");
  assume(num < chrNames.size(),
	 "Genome::getChrName Index " + std::to_string(num) + " out of range (" +
	 std::to_string(chrNames.size()) + ")");
  return chrNames[num];
}

short unsigned int Genome::getLength ( const chr_num_t n ) const {
  //  assume(n < chrLens.size(), "Genome::getLength: index out of range " + std::to_string(n));
  return chrLens.at(n);
}


void Genome::createChromosome ( const std::string name, const uint32_t length ) {
  chr_num_t N = chrNames.size();
  chrNames.push_back(name);
  chrLens.push_back(length);
  chrNums.emplace( name, N );
  readPos->push_back(new std::map<chr_pos_t, p_read_t>());
}


p_read_t Genome::createRead( const chr_num_t chr, const chr_pos_t pos, const short len ) {
  p_read_t rc = getReadAt(chr, pos);
  if (rc == nullptr) {                      // no read on chr:pos exists yet
    rc = std::make_shared<ReadContainer>(chr, pos, len);
    readPos->at(chr)->emplace(pos, rc);          // register read with genome navigation map
  }
  return rc;
}


p_read_t Genome::getReadAt ( const chr_num_t chr, const chr_pos_t pos ) const {
  auto chromosome = readPos->at(chr);
  auto lookup = chromosome->find(pos);
  if (lookup == chromosome->end()) {
    return nullptr;
  }
  return (*lookup).second;
}


p_read_t Genome::getReadAt ( const std::string chr, const chr_pos_t pos ) const {
  return getReadAt( getChrNum(chr), pos );
}


void Genome::read ( std::ifstream& input ) {
  log("Genome: reading file...");
  unsigned int L(1);
  for (std::string line; getline(input, line);) {
    //    std::cout << "Line #" << L << std::endl;
    if (line[0] == '@') {
      assume(parseHeaderLine(line), "Could not parse header in line " + std::to_string(L), false);
    }
    else {
      assume(parseDataLine(line), "Could not parse data in line " + std::to_string(L), false);
    }
    L++;
  }
  if (multistrand != nullptr) {
    log("Found " + std::to_string(multistrand->size()) + " strand switching events");
  }
  if (circular != nullptr) {
    log("Found " + std::to_string(circular->size()) + " circular transcripts");
  }

}


void Genome::read ( std::string& fileName ) {
  std::ifstream input(fileName);
  assume(input.good(), "Could not open file: " + fileName);
  read(input);
}


std::ifstream& operator >> ( std::ifstream& input, Genome& target ) {
  target.read(input);
  return input;
}


bool Genome::parseHeaderLine ( const std::string& line ) {
  try {
    switch (line[1]) {
      case 'H': // Token @HD: General file info
	log(line + "\n" );
	break;
      case 'S': { // Token @SQ: Reference dictionary entry
	  auto tokens = strsplit(line, "\t", false);
	  std::string name = strsplit(tokens[1], ":")[1];
	  unsigned short len = atoi( strsplit(tokens[2], ":")[1].c_str() );
	  createChromosome(name, len);
	}
	break;
      case 'P': // Token @PG: program used
	log(line + "\n");
	assume(strsplit(line, "\t")[1] == "ID:segemehl", "This program is designed for segemehl output only!", false);
      default:  // Ignore other tokens
	break;
    }
  } catch (const std::exception e) { // keep going on errors, but warn
    std::cerr << e.what();
    return false;
  }
  return true;
}


bool Genome::parseDataLine ( const std::string& line ) {
    auto tokens = strsplit(line, "\t");
    unsigned char flags  = 0;
    chr_num_t chr = getChrNum(tokens[RNAME]);
    chr_pos_t pos5 = atoi( tokens[POS].c_str() );
    chr_pos_t pos3;
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
	if ((*iter)[2] != ':' || (*iter)[4] != ':') { // doesn't suffice Xx:y:z, e. g. PHRED-String starting with X
	  continue;
	}
	switch((*iter)[1]) {
	case 'X': // start of current split
	  start = atoi( (*iter).substr(5).c_str() );
	  break;
	case 'Y': // enf of current split
	  end = atoi( (*iter).substr(5).c_str() );
	  break;
	case 'Q': // number of current split
	  flags |= ReadContainer::SPLIT;
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
    
    pos3 = pos5 + calcLength(tokens[CIGAR]);
    p_read_t rc = getReadAt(chr, pos5);
    if (rc == nullptr) {
      rc = std::make_shared<ReadContainer>(chr, pos5, pos3-pos5);
      registerRead(rc);
      rc->flags |= flags;
    }
    
    if (hasNext) {
      rc->addDownstreamRead(*this, nextChr, nextPos);
    }
    if (hasPrev) {
      rc->addUpstreamRead(*this, prevChr, -prevPos);
    }
  return true;
}

void Genome::registerRead (  p_read_t rc ) {
  const chr_num_t chr = rc->chromosome;
  readPos->at(chr)->emplace(rc->fivePrimeEnd, rc);   // Link 5' end to identify as successor
  readPos->at(chr)->emplace(-rc->threePrimeEnd, rc); // Link 3' end to identify as predecessor
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


