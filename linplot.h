#ifndef LINPLOT_H
#define LINPLOT_H

#include "vplot.h"
#include <QGraphicsRectItem>

class LinearPlot : public vPlot {
public:
  LinearPlot();
  
  void fromRead( ReadContainer* seed, const Genome& genome );
  void writeEps( const std::string& fileName ) const;
  void scale( const float factor );  
  
private:
  
  QGraphicsRectItem* fromRead( ReadContainer* seed, const Genome& genome, const int x, const int y );
  uint runID;
  int x_dist;
};

#endif