#ifndef PLOTCHROMOSOME_H
#define PLOTCHROMOSOME_H

#include "genome.h"
#include "readcontainer.h"

#include <string>
#include <map>
#include <memory>
#include <iostream>

struct Rect {
  int x, y, w, h;
};

class PlotChromosome {
 public:
  PlotChromosome ( const size_t, const chr_pos_t, const std::string );
  void addExon ( const p_read_t, int layer = 0 );

  int nExons() const;
  void printout();
  std::shared_ptr<Rect> boundingRect();
  void writeEps ( std::ostream& out, const Rect dim, const int dx, const int dy, const float scale, const float color[3] );
  
  chr_pos_t len;
  
 private:
  void assignIds ();
  
  std::map<chr_pos_t, p_read_t> exons;
  std::string name;
};		  

#endif // PLOTCHROMOSOME_H
