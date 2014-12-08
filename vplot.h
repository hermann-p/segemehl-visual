/**
* Common plotting elements and base class
*/

#ifndef VPLOT_H
#define VPLOT_H

#include "genome.h"
#include <QGraphicsScene>
#include <QGraphicsItemGroup>


/**
 * Representation of a chromosomes
 */
class vPlot;

class PlotChromosome : public QGraphicsItemGroup {
public:
  PlotChromosome ( vPlot* parent, const uint color_hint, const unsigned short len );
  
  QColor color;
  QPen pen;
  QBrush brush;
  
  std::string name;
  std::vector<std::string> exonNames;
  uint32_t length;
  QRect boundingRect;
  
  int fivePrime;
};


/**
 * Abstract base class for plots
 */


class vPlot : public QGraphicsScene {
public:
  vPlot();
  
  virtual PlotChromosome* addChromosome( const Genome& genome, const chr_num_t num );
  virtual void fromRead( ReadContainer* seed, const Genome& genome ) = 0;
  virtual void writeEps( const std::string& fileName ) const = 0;
  virtual void scale( const float factor ) = 0;

  std::map<chr_num_t, PlotChromosome*> chromosomes;
  uint pCol;
  int chrX;
  int chrY;
  int y_dist;
};


#endif // VPLOT_H