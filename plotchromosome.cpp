#include "plotchromosome.h"
#include "genome.h"
#include "vplot.h"
#include "utils.h"

#include <QTransform>

// Chromosome constructor, trying to determine a color automatically
PlotChromosome::PlotChromosome ( vPlot* parent, const uint color_hint, const unsigned short len ) 
  : QGraphicsItemGroup(),
    color(color_hint * 128 % 256,
	  color_hint * 64 % 256,
	  255 - ((color_hint * 64) % 256)
	  ),
    pen(),
    brush(color),
    nextExon('a')
{  
  int x = parent->width() / 2;
  fivePrime = x - len/2;
  int y = (1 + color_hint) * 2 * parent->y_dist;
  QGraphicsLineItem* theChr = parent->addLine(fivePrime, y, x + len/2, y, pen);
  addToGroup(theChr);
  parent->addItem(this);
}


// add an graphical representation for an exon
QGraphicsRectItem* PlotChromosome::addExon ( vPlot* parent, const int p1, const chr_pos_t L ) {
  nextExon = (++nextExon <= 'z') ? nextExon : 'A';
  auto rect1 = parent->addRect( fivePrime + p1,
				childrenBoundingRect().center().y() - parent->y_dist / 4,
				L, parent->y_dist / 2,
				pen, brush
				);
  addToGroup(rect1);
  return rect1;
}


// rescale chromosome representation to fit on screen
void PlotChromosome::fitTo ( const int width, const bool full ) {
  assume(full, "shortening not supported yet", false);
  float scale = 1.0f * width / childrenBoundingRect().width();
  float dx = (childrenBoundingRect().left() + 50);
  QTransform matrix(scale, 0.0f, 0.0f,
		    0.0f, 1.0f, 0.0f,
		    dx*scale, 0.0f, 1.0f); // combine translation and x-scaling
  //  matrix.reset();
  setTransform( matrix );
}
