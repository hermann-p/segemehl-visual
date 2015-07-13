#include "readcontainer.h"
#include "genome.h"
#include "utils.h"

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <memory>

#include <iostream>

ReadContainer::~ReadContainer() {
  if (threePrimeRead != nullptr) delete threePrimeRead;
  if (threePrimeRefs != nullptr) delete threePrimeRefs;
  if (fivePrimeRead != nullptr) delete fivePrimeRead;
  if (moreData != nullptr) delete moreData;
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
  auto dnstr = genome.getReadAt(chr, pos);
  if (dnstr != nullptr) {
    // check if this is a multistrand splice junction
    if (dnstr->chromosome != chromosome && !(dnstr->flags & MULTISTRAND) ) {
      flags |= MULTISTRAND;
      genome.multistrand->insert(shared_from_this());    // adding one multistrand seed suffices
    }

    int link = findLink(dnstr); // search link this => dnstr
    int backLink = dnstr->findLink(shared_from_this(), false); // search link this <= dnstr

    if (link == -1) { // this => dnstr unknown yet
      threePrimeRead->push_back(dnstr);
      threePrimeRefs->push_back(1);
    }
    else { // this => dnstr known
      threePrimeRefs->at(link)++;
    }
    if (backLink == -1) { // this <= dnstr unknown yet
      dnstr->fivePrimeRead->push_back(shared_from_this());
    }
  }
}


void ReadContainer::addUpstreamRead ( const Genome& genome, const chr_num_t chr, const chr_pos_t pos ) {
  assume(pos <= 0, "3' ends need to be given as negative numbers!");
  int ipos = pos - 1; // fix difference in one-based offsets
  auto upstr = genome.getReadAt(chr, ipos); // 3' ends are linked at negative indices
  if (upstr != nullptr) {

    // check if this is a multistrand splice junction
    if (genome.multistrand != nullptr &&
	(upstr->chromosome != chromosome || (upstr->flags & MULTISTRAND) == MULTISTRAND)) {
      upstr->flags |= MULTISTRAND;
      flags |= MULTISTRAND;
      genome.multistrand->insert(shared_from_this());
      genome.multistrand->insert(upstr);
    }
    
    int link = findLink(upstr, false); // search link  upstr <= this
    int backLink = upstr->findLink(shared_from_this()); // search link upstr => this

    if (link == -1) {    // upstr <= this unknown yet
      fivePrimeRead->push_back(upstr);
    }

    if (backLink >= 0) { // upstr => this known
      upstr->threePrimeRefs->at(backLink)++;
    }
    else {               // upstr => this new
      upstr->threePrimeRead->push_back(shared_from_this());
      upstr->threePrimeRefs->push_back(1);
    }
  }
}


int ReadContainer::findLink ( std::shared_ptr<ReadContainer> partner, const bool downstream ) {
  if (downstream && threePrimeRead == nullptr) {
    threePrimeRead = new link_list_t;
    threePrimeRefs = new std::vector<unsigned short>;
    return -1;
  }
  else if (!downstream && fivePrimeRead == nullptr) {
    fivePrimeRead = new link_list_t;
    return -1;
  }
  auto candidates = (downstream) ? threePrimeRead : fivePrimeRead;
  for (size_t i(0); i < candidates->size(); i++) {
    if (*candidates->at(i) == *partner) {
      return i;
    }
  }
  return -1;
}


ReadContainer::ReadContainer ( ) 
: threePrimeRead(nullptr),
  fivePrimeRead(nullptr),
  threePrimeRefs(nullptr),
  flags(0),
  moreData(nullptr)
{
}


ReadContainer::ReadContainer ( const chr_num_t chromosome, const chr_pos_t p5, const short int len ) :
  threePrimeRead(nullptr),
  fivePrimeRead(nullptr),
  threePrimeRefs(nullptr),
  chromosome(chromosome),
  fivePrimeEnd(p5),
  threePrimeEnd(p5 + len),
  flags(0),
  moreData(nullptr)
{
  if (len < 0) {
    flags |= REVERSE;
  }
}


std::ostream& operator << ( std::ostream& output, ReadContainer& rc ) {
  if (rc.flags & ReadContainer::DUMMY) {
    output << "dummy read";
  }
  else {
    output << std::to_string(rc.chromosome) << ":" << std::to_string(rc.fivePrimeEnd) << "-" <<
      std::to_string(rc.threePrimeEnd);
  } 
  return output;
}
