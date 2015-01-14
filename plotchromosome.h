/**
 * Representation of a chromosomes
 */

#ifndef PLOTCHROMOSOME_H
#define PLOTCHROMOSOME_H

#include <QColor>
#include <QPen>
#include <QRect>
#include <QBrush>
#include <QGraphicsItemGroup>
#include <QGraphicsRectItem>

#include "genome.h"

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
  char nextExon;

  QGraphicsRectItem* addExon ( vPlot* parent, const chr_pos_t p1, const chr_pos_t p2 );
  void fitTo ( const int width, const bool full=true );
  
  int fivePrime;
};


#endif // PLOTCHROMOSOME_H
