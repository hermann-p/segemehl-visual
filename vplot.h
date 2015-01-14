/**
* Common plotting elements and base class
*/

#ifndef VPLOT_H
#define VPLOT_H

#include "genome.h"
#include "plotchromosome.h"

#include <QGraphicsScene>
#include <QGraphicsItemGroup>

#include <memory>




/**
 * Abstract base class for plots
 */


class vPlot : public QGraphicsScene {
public:
  vPlot();
  
  virtual PlotChromosome* addChromosome( const Genome* genome, const chr_num_t num );
  virtual void fromRead( std::shared_ptr<ReadContainer> seed, const Genome* genome ) = 0;
  virtual void writeEps( const std::string& fileName ) const = 0;
  virtual void scale( const float factor ) = 0;
  virtual void fitTo( const int w, const int h ) = 0;

  std::map<chr_num_t, PlotChromosome*> chromosomes;
  uint pCol;
  int chrX;
  int chrY;
  int y_dist;
};


#endif // VPLOT_H
