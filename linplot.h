#ifndef LINPLOT_H
#define LINPLOT_H

#include "vplot.h"
#include <QGraphicsRectItem>
#include <memory>

class LinearPlot : public vPlot {
public:
  LinearPlot();
  
  void fromRead ( std::shared_ptr<ReadContainer> seed, const Genome* genome );
  void writeEps ( const std::string& fileName ) const;
  void scale ( const float factor );  
  void fitTo ( const int w, const int h );
  
private:
  
  QGraphicsRectItem* fromRead ( std::shared_ptr<ReadContainer> seed, const Genome* genome, const   int x, const int y, const int dir = BOTH );
  uint runID;
  int x_dist;
  
  static enum {
    BACKWARDS = -1, BOTH = 0, FORWARD = 1
  } EXTENSION_DIRECTIONS;
};

#endif
