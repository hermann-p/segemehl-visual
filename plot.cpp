#include "plot.h"
#include "utils.h"

#include <QRect>
#include <vector>

// Construct a line and connect it with 2 exons
LinearConnector::LinearConnector ( QGraphicsRectItem* const r1, QGraphicsRectItem* const r2 ) :
  a(r1),
  b(r2)
{
  recalculate();
}


// Move line endpoints to centers of connected exons
void LinearConnector::recalculate() {
  setLine( a->boundingRect().center().x(),
	   a->boundingRect().center().y(),
	   b->boundingRect().center().x(),
	   b->boundingRect().center().y()
  );
}



QRect LinearPlot::getBoundingRect() const {
  assume(!chromosomes.empty(), "Can not determine bounding rect, plot data empty");
  QRect bounds = chromosomes.at(0).boundingRect;
//  for (auto& chr : chromosomes) {
//    bounds |= chr.boundingRect;
//  }
  return bounds;
}


void LinearPlot::fromRead ( ReadContainer* seed, const Genome& genome ) {
  if (seed == nullptr) {
    log("No seed given, can't plot");
    return;
  }
  fromRead(seed, genome, 0, 0);
  runID++;
}


void LinearPlot::scale ( const float f ) {
}


QGraphicsRectItem* LinearPlot::fromRead ( ReadContainer* seed, const Genome& genome, const int xpos, const int ypos ) {
  
  // Create a chromosome (aka exon container)
  PlotChromosome* chr = addChromosome(genome, seed->chromosome);
  
  // Add this exon to the chromosome
  QGraphicsRectItem* chrExon = addRect(chr->x() + seed->fivePrimeEnd,
				       chr->y(),
				       seed->length,
				       chrHeight,
				       chr->pen, chr->brush
			       );
  chr->addToGroup(chrExon);
  
  // Add this exon
  QGraphicsRectItem* exon = addRect(xpos - seed->length / 2.0, // left
				    chr->y(),                  // bottom
				    seed->length, 
				    chrHeight, 
				    chr->pen, chr->brush
			    );
  seed->runID = runID;
  
  if (!seed->fivePrimeRead.empty()) {
    float dy = seed->fivePrimeRead.size() / 2;
    if (seed->fivePrimeRead.size() % 2 == 0) {
      dy -= 0.5f;
    }
    for (auto &prev : seed->fivePrimeRead) {
      if (prev->runID != runID) {
	LinearConnector lc( exon,
			    fromRead(prev, genome, 
				     xpos - minDist - (seed->length + prev->length) / 2.0,
				     ypos + dy * 3 * chrHeight
 				    )
			  );
	dy += 1.0f;
	connectors.push_back(lc);
      }
    }
  }
  
  if (!seed->threePrimeRead.empty()) {
    for (auto &next : seed->threePrimeRead) {
      float dy = seed->threePrimeRead.size() / 2;
      if (seed->threePrimeRead.size() % 2 == 0) {
	dy -= 0.5f;
      }
      if (next->runID != runID) {
	LinearConnector lc( exon,
			    fromRead(next, genome, 
				     xpos + minDist + (seed->length + next->length) / 2.0,
				     ypos + dy * 3 * chrHeight
				    )
			  );
	dy += 1.0f;
      }
    }
  }
  
  return exon;
}


PlotChromosome* vPlot::addChromosome ( const Genome& genome, const chr_num_t num ) {
  auto lookup = chromosomes.find(num);
  if (lookup != chromosomes.end()) {
    return &(lookup->second);
  }
  PlotChromosome* chr = new PlotChromosome();
  chr->length = genome.getLength(num);
  
  return chr;
}


void LinearPlot::writeEps ( const std::string& fileName ) const {
}
