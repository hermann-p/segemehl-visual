#ifndef LINPLOT_H
#define LINPLOT_H

#include "vplot.h"
#include <QGraphicsRectItem>
#include <memory>

class LinearPlot : public vPlot {
public:
  LinearPlot();
  
  void fromRead ( std::shared_ptr<ReadContainer> seed, const std::shared_ptr<Genome> genome );
  void writeEps ( const std::string& fileName ) const;
  void scale ( const float factor );  
  
private:
  
//  QGraphicsRectItem* fromRead ( ReadContainer* seed, const Genome& genome, const int x, const int y );
  QGraphicsRectItem* fromRead ( std::shared_ptr<ReadContainer> seed, const std::shared_ptr<Genome> genome, const int x, const int y, const int dir = BOTH );
  uint runID;
  int x_dist;
  
  enum EXTENSION_DIRECTIONS {
    BACKWARDS = -1, BOTH = 0, FORWARD = 1
  };
};

#endif