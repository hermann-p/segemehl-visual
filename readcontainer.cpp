#include "readcontainer.h"
#include "genome.h"
#include "utils.h"

#include <vector>
#include <string>
#include <map>
#include <unordered_map>

bool operator == ( const ReadContainer& a, const ReadContainer& b ) {
  return (a.chromosome == b.chromosome &&
          a.fivePrimeEnd == b.fivePrimeEnd &&
          a.threePrimeEnd == b.threePrimeEnd);
}


short int ReadContainer::length() const {
  return threePrimeEnd - fivePrimeEnd;
}


void ReadContainer::addDownstreamRead ( const Genome& genome, const chr_num_t chr, const chr_pos_t pos ) {
  ReadContainer* dnstr = genome.getReadAt(chr, pos);
  if (dnstr != nullptr) {
    int link = dnstr->findLink(this);
    if (link >= 0) { // this read supports an existing splice site
      dnstr->threePrimeRefs[link]++;
    }
    else { // this splicing was not yet known
      dnstr->threePrimeRead.push_back(this);
      dnstr->threePrimeRefs.push_back(1);
    }
  }
}


void ReadContainer::addUpstreamRead ( const Genome& genome, const chr_num_t chr, const chr_pos_t pos ) {
  assume(pos <= 0, "3' ends need to be given as negative numbers!");
  ReadContainer* upstr = genome.getReadAt(chr, pos); // 3' ends are linked at negative indices
  if (upstr != nullptr) {
    int link = upstr->findLink(this, false);
    if (link >= 0) {
      threePrimeRefs[link]++;
    }
    else {
      upstr->fivePrimeRead.push_back(this);
      threePrimeRead.push_back(upstr);
      threePrimeRefs.push_back(1);
    }
  }
}


int ReadContainer::findLink ( const ReadContainer* partner, const bool fwd ) const {
  auto candidates = (fwd) ? threePrimeRead : fivePrimeRead;
  for (size_t i(0); i < candidates.size(); i++) {
    if (*candidates[i] == *partner) {
      return i;
    }
  }
  return -1;
}


ReadContainer::ReadContainer ( ) :
  runID(0)
{
}



ReadContainer::ReadContainer ( const chr_num_t chromosome, const chr_pos_t p5, const short int len ) :
  chromosome(chromosome),
  fivePrimeEnd(p5),
  threePrimeEnd(p5 + len),
  runID(0)
{
  if (len < 0) {
    flags &= REVERSE;
  }
}