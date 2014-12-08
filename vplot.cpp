#include "vplot.h"

#include <iostream>
#include <QGraphicsLineItem>
#include <QGraphicsTextItem>
#include <QString>

vPlot::vPlot()
  : QGraphicsScene(),
  pCol(0),
  y_dist(20)
{
}

// Chromosome constructor, trying to determine a color automatically
PlotChromosome::PlotChromosome ( vPlot* parent, const uint color_hint, const unsigned short len ) 
  : QGraphicsItemGroup(),
  color(color_hint * 128 % 256,
        color_hint * 64 % 256,
	255 - ((color_hint * 64) % 256)
  ),
  pen(),
  brush(color)
{  
  int x = parent->width() / 2;
  fivePrime = x - len/2;
  int y = (1 + color_hint) * 2 * parent->y_dist;
  QGraphicsLineItem* theChr = parent->addLine(fivePrime, y, x + len/2, y, pen);
  addToGroup(theChr);
  parent->addItem(this);
}


// Add a chromosome based on information from genome
PlotChromosome* vPlot::addChromosome ( const Genome& genome, const chr_num_t num ) {
    auto lookup = chromosomes.find(num);
  if (lookup != chromosomes.end()) {
    return (lookup->second);
  }
  PlotChromosome* chr = new PlotChromosome(this, pCol++, genome.getLength(num));
  chromosomes.emplace(num, chr);
  chr->length = genome.getLength(num);
  QGraphicsTextItem* name = addText(genome.getChrName(num).c_str());
  name->setPos(chr->childrenBoundingRect().x() - 50, chr->childrenBoundingRect().center().y());
  name->moveBy(0, -12);
//  chr->addToGroup(name);
//  addItem(name);
  return chr;
}
