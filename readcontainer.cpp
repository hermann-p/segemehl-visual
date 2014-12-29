#include "readcontainer.h"
#include "genome.h"
#include "utils.h"

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <memory>

ReadContainer::~ReadContainer() {
  if (threePrimeRead != nullptr) delete threePrimeRead;
  if (threePrimeRefs != nullptr) delete threePrimeRefs;
  if (fivePrimeRead != nullptr) delete fivePrimeRead;
}


bool operator == ( const ReadContainer& a, const ReadContainer& b ) {
  return (a.chromosome == b.chromosome &&
          a.fivePrimeEnd == b.fivePrimeEnd &&
          a.threePrimeEnd == b.threePrimeEnd);
}


short int ReadContainer::length() const {
  return threePrimeEnd - fivePrimeEnd;
}

void ReadContainer::addDownstreamRead ( const Genome& genome, const chr_num_t chr, const chr_pos_t pos ) {
  std::shared_ptr<ReadContainer> dnstr(genome.getReadAt(chr, pos));
  if (dnstr != nullptr) {
    if (fivePrimeRead == nullptr) {
      fivePrimeRead = new link_list_t;
    }
    int link = dnstr->findLink(this);
    if (link >= 0) { // this read supports an existing splice site
      dnstr->threePrimeRefs->at(link)++;
    }
    else { // this splicing was not yet known
      dnstr->threePrimeRead->push_back(shared_from_this());
      dnstr->threePrimeRefs->push_back(1);
    }
  }
}


void ReadContainer::addUpstreamRead ( const Genome& genome, const chr_num_t chr, const chr_pos_t pos ) {
  assume(pos <= 0, "3' ends need to be given as negative numbers!");
  std::shared_ptr<ReadContainer> upstr(genome.getReadAt(chr, pos)); // 3' ends are linked at negative indices
  if (upstr != nullptr) {
    if (threePrimeRefs == nullptr) {
      threePrimeRead = new link_list_t();
      threePrimeRefs = new std::vector<unsigned short>;
    }
    int link = upstr->findLink(this, false);
    if (link >= 0) {
      threePrimeRefs->at(link)++;
    }
    else {
      upstr->fivePrimeRead->push_back(shared_from_this());
      threePrimeRead->push_back(upstr);
      threePrimeRefs->push_back(1);
    }
  }
}


int ReadContainer::findLink ( ReadContainer* partner, const bool fwd ) {
  if (fwd && threePrimeRead == nullptr) {
    threePrimeRead = new link_list_t;
    threePrimeRefs = new std::vector<unsigned short>;
    return -1;
  }
  else if (!fwd && fivePrimeRead == nullptr) {
    fivePrimeRead = new link_list_t;
    return -1;
  }
  auto candidates = (fwd) ? threePrimeRead : fivePrimeRead;
  for (size_t i(0); i < candidates->size(); i++) {
    if (*candidates->at(i) == *partner) {
      return i;
    }
  }
  return -1;
}


ReadContainer::ReadContainer ( ) 
: runID(0),
threePrimeRead(nullptr),
fivePrimeRead(nullptr),
threePrimeRefs(nullptr)
{
}



ReadContainer::ReadContainer ( const chr_num_t chromosome, const chr_pos_t p5, const short int len ) :
  chromosome(chromosome),
  fivePrimeEnd(p5),
  threePrimeEnd(p5 + len),
  runID(0),
  threePrimeRead(nullptr),
  fivePrimeRead(nullptr),
  threePrimeRefs(nullptr)
{
  if (len < 0) {
    flags &= REVERSE;
  }
}