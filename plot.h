#ifndef PLOT_H
#define PLOT_H

#include "genome.h"
#include "readcontainer.h"

#include <QRect>
#include <QLine>
#include <QColor>
#include <QGraphicsScene>
#include <QGraphicsItemGroup>
#include <QGraphicsRectItem>
#include <QLine>
#include <string>
#include <vector>
#include <fstream>



/**
 * Data content for plots
 **/

class PlotChromosome : public QGraphicsItemGroup {
public:
  QColor color;
  QPen pen;
  QBrush brush;
  
  std::string name;
  std::vector<std::string> exonNames;
  uint32_t length;
  QRect boundingRect;
};


/**
* Abstract interface class for plots
**/
class vPlot : public QGraphicsScene {
public:
  virtual PlotChromosome* addChromosome( const Genome& genome, const chr_num_t num );
  virtual void fromRead( ReadContainer* seed, const Genome& genome ) = 0;
  virtual void writeEps( const std::string& fileName ) const = 0;
  virtual void scale( const float factor ) = 0;

  std::map<chr_num_t, PlotChromosome> chromosomes;
};




/**
* Content helper classes for LinearPlot
**/


class LinearConnector : QLine {
public:
  LinearConnector ( QGraphicsRectItem* const r1, QGraphicsRectItem* const r2 );
  QGraphicsRectItem* a;
  QGraphicsRectItem* b;
  void recalculate();
};


/**
* Class to plot transcript isoforms in a tree-like graph.
**/
class LinearPlot : public vPlot {
public:
  void fromRead ( ReadContainer* seed, const Genome& genome ); // create a plot by extending a single read to both directions
  void writeEps ( const std::string& fileName ) const;
  void scale ( const float f );

private:
  QGraphicsRectItem* fromRead ( ReadContainer* seed, const Genome& genome, const int xpos = 0, const int ypos = 0 );
  const int minDist = 30;
  const int chrHeight = 30;
  std::vector<LinearConnector> connectors;
  uint runID;
  
//  PlotChromosome* addChromosome ( const Genome& genome, const chr_num_t num );
  QRect getBoundingRect () const;
};


#endif // PLOT_H