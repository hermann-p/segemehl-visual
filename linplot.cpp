#include "linplot.h"
#include "readcontainer.h"
#include "genome.h"

#include <QGraphicsRectItem>

#include <iostream>

LinearPlot::LinearPlot()
  : vPlot(),
  runID(1),
  x_dist(10)
{
}


void LinearPlot::fromRead ( std::shared_ptr<ReadContainer> seed, const Genome* genome ) {
  fromRead(seed, genome, width()/2, height()/2);
//  int biggest = 0;
  for (auto& chr : chromosomes) {
//    auto w = chr.childrenBoundingRect().width();
//    biggest = (w > biggest) ? w : biggest;
    chr.second->fitTo( 1024 );
  }
  
  runID++;
}


// Improved version that avoids zig-zagging through forward- and backward-links
// Which could lead to splice variants without the original seed. Call with dir = 0 (default)
// to start in both directions.
// Currently direction discrimination has been temporary has been disabled to activate zig-zag.
QGraphicsRectItem* LinearPlot::fromRead ( std::shared_ptr<ReadContainer> seed, const Genome* genome, const int x, const int y, const int dir ) {
  seed->flags |= ReadContainer::PROCESSED;
  
  // Create graphics item for this exon
  PlotChromosome* chr = addChromosome(genome, seed->chromosome);
  char exonId = chr->nextExon;
  QGraphicsRectItem* r1 = addRect(x - seed->length()/2, y - y_dist / 4,
	  seed->length(), y_dist/2,
	  chr->pen, chr->brush
 	);
  QGraphicsRectItem* r2 = chr->addExon(this, seed->fivePrimeEnd, seed->length());
  
  // Traverse through predecessors
  if (dir <= 0 && seed->fivePrimeRead != nullptr) {
    // determine vertical position
    int y0 = y - (seed->fivePrimeRead->size() / 2) * y_dist;
    if (seed->fivePrimeRead->size() % 2 == 0) {
      y0 += y_dist / 2;
    }
    for (auto& el : *(seed->fivePrimeRead)) {
      if (!(el->flags & ReadContainer::PROCESSED)) {
	int xpos = x + x_dist + seed->length()/2 + el->length()/2;
	fromRead(el, genome, xpos, y0, BOTH);
	y0 += y_dist;
      }
    }
  }
    
  // Traverse through successors
  if (dir >= 0 && seed->threePrimeRead != nullptr) {
    // determine vertical position
    int y0 = y - (seed->threePrimeRead->size() / 2) * y_dist;
    if (seed->threePrimeRead->size() % 2 == 0) {
      y0 += y_dist / 2;
    }
    for (auto& el : *(seed->threePrimeRead)) {
      if (!(el->flags & ReadContainer::PROCESSED)) {
	int xpos = x - x_dist - seed->length()/2 - el->length()/2;	
	fromRead(el, genome, xpos, y0, BOTH);
	y0 += y_dist;
      }
    }
  }
 
  return r1;
}

void LinearPlot::fitTo ( const int w, const int h ) {
  
}

void LinearPlot::writeEps ( const std::string& fileName ) const {
}


void LinearPlot::scale ( const float factor ) {
}
