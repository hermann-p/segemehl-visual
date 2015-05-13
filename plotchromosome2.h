#ifndef PLOTCHROMOSOME_H
#define PLOTCHROMOSOME_H

#include "genome.h"

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
  void addExon ( chr_pos_t p1, chr_pos_t p2 );

  void printout();
  std::shared_ptr<Rect> boundingRect();
  void writeEps ( std::ostream& out, const Rect dim, const int dx, const int dy, const float scale, const float color[3] ) const;
  
  chr_pos_t len;
  
 private:
  std::map<chr_pos_t, chr_pos_t> exons;
  std::string name;
};		  

#endif // PLOTCHROMOSOME_H
