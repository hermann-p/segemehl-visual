#ifndef LINPLOT_H
#define LINPLOT_H

#include "vplot2.h"
#include "utils.h"

#include <memory>
#include <vector>

struct Exon {
  int lt;
  int rt;
  int y;
  char id;
  chr_num_t chr;
};

struct LineEnds {
  std::shared_ptr<Exon> a;
  std::shared_ptr<Exon> b;
  int num;
};
    

class LinearPlot : public vPlot {
 public:
  LinearPlot(int dx=40, int dy=40);

  void fromRead ( std::shared_ptr<ReadContainer> seed, Genome* genome );
  void writeEps ( const std::string& fileName ) const;

 private:
  void fromRead ( std::shared_ptr<ReadContainer>, Genome*, const int x, const int y, std::shared_ptr<Exon> pred, const int nLinks = -1 );
  std::vector<LineEnds> connections;
  std::vector< std::shared_ptr<Exon> > positions; // [lt, rt, y]
  std::shared_ptr<Rect> boundingRect() const;
  
  int dx, dy;
  int minx, maxx, miny, maxy;
  uint nFilter;

  char nextID;
};



#endif // LINPLOT_H
