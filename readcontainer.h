#ifndef READCONTAINER_H
#define READCONTAINER_H

#include "genome.h"

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <memory>

class ReadContainer : public std::enable_shared_from_this<ReadContainer> {
public:
  ReadContainer ();
  ReadContainer( const chr_num_t chromosome, const chr_pos_t p5, const short len );
  std::vector<std::shared_ptr<ReadContainer>> fivePrimeRead;
  std::vector<std::shared_ptr<ReadContainer>> threePrimeRead;
  std::vector<unsigned short> threePrimeRefs;
  chr_num_t chromosome;
  chr_pos_t fivePrimeEnd;
  chr_pos_t threePrimeEnd;
  unsigned char flags;
  short length() const;
  uint runID;
  
  void addUpstreamRead( const Genome& genome, const chr_num_t chr, const chr_pos_t pos );
  void addDownstreamRead( const Genome& genome, const chr_num_t chr, const chr_pos_t pos );

private:  
  int findLink( const ReadContainer* partner, const bool fwd = true ) const;
  enum FLAGS {
    REVERSE = 1,
    CIRCULAR = 2
  };
};

bool operator == ( const ReadContainer& a, const ReadContainer& b );

#endif // READCONTAINER_H