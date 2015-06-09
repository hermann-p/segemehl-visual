#ifndef READCONTAINER_H
#define READCONTAINER_H

#include "genome.h"

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <memory>


// additional data required to construct visible plot elements and
// arrange them in the output
struct PlotInfo {
  uint id;
  int layer; // only for linear plots
  // coordinates      linear | circular
  float c1;        // left   | start angle  
  int c2;        // y      | distance
  float c3;        // right  | end angle
};


typedef std::shared_ptr<ReadContainer> p_read_t;


class ReadContainer : public std::enable_shared_from_this<ReadContainer> {
public:
  // create empty container
  ReadContainer ();
  ReadContainer( const chr_num_t chromosome, const chr_pos_t p5, const short len );
  ~ReadContainer();
  
  void addUpstreamRead( const Genome& genome, const chr_num_t chr, const chr_pos_t pos );
  void addDownstreamRead( const Genome& genome, const chr_num_t chr, const chr_pos_t pos );
  friend std::ostream& operator<<( std::ostream& output, ReadContainer& rc );

  link_list_t* fivePrimeRead;
  link_list_t* threePrimeRead;
  std::vector<unsigned short>* threePrimeRefs;
  chr_num_t chromosome;
  chr_pos_t fivePrimeEnd;
  chr_pos_t threePrimeEnd;
  unsigned char flags;
  short length() const;
//  uint runID;

  static enum {
    REVERSE = 1,
    CIRCULAR = 2,
    SPLIT = 4,
    MULTISTRAND = 8,
    DUMMY = 64,
    PROCESSED = 128
  } FLAGS;

  PlotInfo* moreData;
  int findLink( std::shared_ptr<ReadContainer> partner, const bool downstream = true );
};

bool operator == ( const ReadContainer& a, const ReadContainer& b );

#endif // READCONTAINER_H
