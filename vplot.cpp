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

// Add a chromosome based on information from genome
PlotChromosome* vPlot::addChromosome ( const Genome* genome, const chr_num_t num ) {
    auto lookup = chromosomes.find(num);
  if (lookup != chromosomes.end()) {
    return (lookup->second);
  }
  PlotChromosome* chr = new PlotChromosome(this, pCol++, genome->getLength(num));
  chromosomes.emplace(num, chr);
  chr->length = genome->getLength(num);
  QGraphicsTextItem* name = addText(genome->getChrName(num).c_str());
  name->setPos(chr->childrenBoundingRect().x() - 50, chr->childrenBoundingRect().center().y());
  name->moveBy(0, -12);
//  chr->addToGroup(name);
//  addItem(name);
  return chr;
}
