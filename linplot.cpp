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


void LinearPlot::fromRead ( ReadContainer* seed, const Genome& genome ) {
  fromRead(seed, genome, width()/2, height()/2);
  runID++;
}

QGraphicsRectItem* LinearPlot::fromRead ( ReadContainer* seed, const Genome& genome, const int x, const int y ) {
  seed->runID = runID;
  
  // Create graphics item for this exon
  PlotChromosome* chr = addChromosome(genome, seed->chromosome);
  QGraphicsRectItem* r1 = addRect(x - seed->length()/2, y - y_dist / 4,
	  seed->length(), y_dist/2,
	  chr->pen, chr->brush
 	);
  QGraphicsRectItem* r2 = addRect(chr->fivePrime + seed->fivePrimeEnd, chr->childrenBoundingRect().center().y() - y_dist / 4,
				  seed->length(), y_dist / 2,
				  chr->pen, chr->brush
			  );
  chr->addToGroup(r2);
  
  // Traverse through predecessors
  if (!seed->fivePrimeRead.empty()) {
    // determine vertical position
    int y0 = y - (seed->fivePrimeRead.size() / 2) * y_dist;
    if (seed->fivePrimeRead.size() % 2 == 0) {
      y0 += y_dist / 2;
    }
    for (auto& el : seed->fivePrimeRead) {
      if (el->runID != runID) {
	int xpos = x + x_dist + seed->length()/2 + el->length()/2;
	fromRead(el, genome, xpos, y0);
	y0 += y_dist;
      }
    }
  }
    
  // Traverse through successors
  if (!seed->threePrimeRead.empty()) {
    // determine vertical position
    int y0 = y - (seed->threePrimeRead.size() / 2) * y_dist;
    if (seed->threePrimeRead.size() % 2 == 0) {
      y0 += y_dist / 2;
    }
    for (auto& el : seed->threePrimeRead) {
      if (el->runID != runID) {
	int xpos = x - x_dist - seed->length()/2 - el->length()/2;	
	fromRead(el, genome, xpos, y0);
	y0 += y_dist;
      }
    }
  }
 
  return r1;
}



void LinearPlot::writeEps( const std::string& fileName ) const {
}


void LinearPlot::scale( const float factor ) {
}