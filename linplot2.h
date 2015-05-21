#ifndef LINPLOT_H
#define LINPLOT_H

#include "vplot2.h"
#include "utils.h"

#include <memory>
#include <vector>
#include <map>


class LinearPlot : public vPlot {
 public:
  LinearPlot(int dx=40, int dy=40);

  void fromRead ( std::shared_ptr<ReadContainer> seed, Genome* genome );
  void writeEps ( const std::string& fileName );
  void createPlotCoords ();

 private:
  void insertDummies ();
  std::vector<std::vector<bool>> transitiveReduction ();
  void barycenterCoords ();
  
  Rect boundingBox;
  std::map<chr_pos_t, std::vector<std::shared_ptr<ReadContainer>>> layeredGraph;
  std::vector<std::shared_ptr<ReadContainer>> flatGraph;
  std::shared_ptr<Rect> boundingRect() const;

  int dx, dy;
  int minx, maxx, miny, maxy;
  uint nFilter;

  char nextID;
};



#endif // LINPLOT_H
