/* 
 * Class that contains reassembled reads and their relation to the chromosomal
 * positions of all individual fragments, as found by the segemehl matching tool.
 * Assumption: all elements sharing a given 5'/3' position also share the
 * 3'/5' position. This may 33be terribly wrong.
 */


#include "genome.h"
#include "readcontainer.h"
#include "utils.h"
#include "lockfreequeue.h"

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <fstream>
#include <cstdlib>
#include <memory>

#include <thread>

#include <sstream>
#include <cstdio>
#include <iostream>

Genome::Genome ( bool init_m, bool init_c ) 
  : multistrand((init_m) ? new link_set_t : nullptr),
    circular((init_c) ? new link_set_t : nullptr),
    readPos(new chr_map_t)
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


// returns the internal id of a chromosome with given name string
chr_num_t Genome::getChrNum ( const std::string name ) const {
  auto lookup = chrNums.find(name);
  assume(lookup != chrNums.end(), "Genome::getChrNum: Chromosome " + name + " not found");
  return lookup->second;
}


// returns the string name for a chromosome's given internal id
std::string Genome::getChrName ( const chr_num_t num ) const {
  assume(!chrNames.empty(), "Genome::getChrName: name list empty");
  assume(num < chrNames.size(),
	 "Genome::getChrName Index " + std::to_string(num) + " out of range (" +
	 std::to_string(chrNames.size()) + ")");
  return chrNames[num];
}


// returns the length of a specific chromosome
short unsigned int Genome::getLength ( const chr_num_t n ) const {
  //  assume(n < chrLens.size(), "Genome::getLength: index out of range " + std::to_string(n));
  return chrLens.at(n);
}


// create a new chromosome and add it to appropriate lists
void Genome::createChromosome ( const std::string name, const uint32_t length ) {
  chr_num_t N = chrNames.size();
  chrNames.push_back(name);
  chrLens.push_back(length);
  chrNums.emplace( name, N );
  readPos->push_back(new std::map<chr_pos_t, p_read_t>());
}


// creates a new read and registers it with the genome map
/*p_read_t Genome::createRead( const chr_num_t chr, const chr_pos_t pos, const short len ) {
  p_read_t rc = getReadAt(chr, pos);
  if (rc == nullptr) {                      // no read on chr:pos exists yet
    rc = std::make_shared<ReadContainer>(chr, pos, len);
    readPos->at(chr)->emplace(pos, rc);          // register read with genome navigation map
  }
  return rc;
}*/


// searches for an existing read by internal chromosome id and position
p_read_t Genome::getReadAt ( const chr_num_t chr, const chr_pos_t pos ) const {
  auto chromosome = readPos->at(chr);
  auto lookup = chromosome->find(pos);
  if (lookup == chromosome->end()) {
    return nullptr;
  }
  return lookup->second;
}


// allows to search for an existing read by chromosome name and position
p_read_t Genome::getReadAt ( const std::string chr, const chr_pos_t pos ) const {
  return getReadAt( getChrNum(chr), pos );
}


// parses an input stream linewise for data
// uses separate threads for reading&tokenizing and reassembly
// backend called by all other implementations of read()
void Genome::read ( std::ifstream& input ) {
  log("Genome: reading file...");
  unsigned int L(1);                         // line counter
  LockFreeQueue<std::vector<std::string>> q;
  
  auto readerTask = [&q,&input]() {
    debug("Starting reader thread...");
      for (std::string line; getline(input, line);) {
	q.push(strsplit(line, "\t"));
      }
      q.done();
      debug("All data read and tokenized, closing reader rhread");
  };
  
  auto parserTask = [&q,&L,this]() {
    debug("Starting parser thread...");
    std::vector<std::string> splits;
    while(q.pop(splits)) {
      assume(parseDataLine(splits), "Could not parse data in line " + std::to_string(L), false);
      L++;
    }
    debug("All tokens processed, closing parser thread");
  };
  
  // read and parse header lines sequentially
  for (std::string line; getline(input, line);) {
    //    std::cout << "Line #" << L << std::endl;
    if (line[0] == '@') {
      assume(parseHeaderLine(line), "Could not parse header in line " + std::to_string(L), false);
    }
    else {
      // process current data line, else it would get lost
      auto tokens = strsplit(line, "\t");
      parseDataLine(tokens);
      break;  // header lines parsed, break and begin threaded processing
    }
    L++;
  }
  
  // init and start threads
  std::thread readerThread(readerTask);
  std::thread parserThread(parserTask);
  
  // wait for threads to finish
  readerThread.join();
  parserThread.join();
  
  if (multistrand != nullptr) {
    log("Found " + std::to_string(multistrand->size()) + " strand switching events");
  }
  if (circular != nullptr) {
    log("Found " + std::to_string(circular->size()) + " circular transcripts");
  }

}


// opens a stream to input file of given name and tries to read it
void Genome::read ( std::string& fileName ) {
  std::ifstream input(fileName);
  assume(input.good(), "Could not open file: " + fileName);
  read(input);
}


// reads a genome from an existing input stream
std::ifstream& operator >> ( std::ifstream& input, Genome& target ) {
  target.read(input);
  return input;
}


// extracts information from a header string
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


// extracs read data from a tokenized line to reassemble RNA
bool Genome::parseDataLine ( std::vector<std::string>& tokens ) {
    unsigned char flags  = 0;
    chr_num_t chr = getChrNum(tokens[RNAME]);
    chr_pos_t pos5 = atoi( tokens[POS].c_str() );
//    chr_pos_t pos3;
    bool hasNext = false;
    chr_num_t nextChr = 0;
    chr_pos_t nextPos = 0;
    bool hasPrev = false;
    chr_num_t prevChr = 0;
    chr_pos_t prevPos = 0;
    auto iter = tokens.begin();
    for (int i(0); i < QUAL; i++) iter++;             // skip default SAM tags
    for (;iter != tokens.end(); iter++) {             // iterate through custom tags
      if ((*iter)[0] == 'X') {                        // segemehl's tags start with X
	if ((*iter)[2] != ':' || (*iter)[4] != ':') { // doesn't suffice Xx:y:z, e. g. PHRED-String starting with X
	  continue;
	}  // guardian to filter non-segemehl tags starting with X
	switch((*iter)[1]) {
/*	You can't trust start and end information to be set correctly, or, at all
 * 	case 'X': // start of current split
	  start = atoi( (*iter).substr(5).c_str() );
	  break;
	case 'Y': // end of current split
	  end = atoi( (*iter).substr(5).c_str() );
	  break; */
	case 'Q': // number of current split
	  flags |= ReadContainer::SPLIT;              // read consists of multiple splits
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

//    pos3 = pos5 + calcLength(tokens[CIGAR]); // recalculate original position on chromosomes by undoing assumed indels 
    p_read_t rc = getReadAt(chr, pos5);
    if (rc == nullptr) {                     // the read os not yet contained in the genome
      rc = std::make_shared<ReadContainer>(chr, pos5, calcLength(tokens[CIGAR]));
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


// store created read in chromosome map
void Genome::registerRead (  p_read_t rc ) {
  const chr_num_t chr = rc->chromosome;
  readPos->at(chr)->emplace(rc->fivePrimeEnd, rc);   // Link 5' end to identify as successor
  readPos->at(chr)->emplace(-rc->threePrimeEnd, rc); // Link 3' end to identify as predecessor
}


// create original size of the read, before matching algorithm added indels
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
