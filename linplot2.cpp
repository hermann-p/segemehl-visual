#include "linplot2.h"
#include "readcontainer.h"

#include <fstream>
#include <math.h>
#include <queue>

LinearPlot::LinearPlot ( int dx, int dy ) 
  : vPlot(),
    dx(dx), dy(dy),
    minx(0), maxx(0), miny(0), maxy(0),
    nextID('a'),
    nFilter(0)
{
}


void LinearPlot::fromRead ( std::shared_ptr<ReadContainer> seed, Genome* genome ) {
  class ClosureCnt { // Container to keep a read and its on-screen positions together
  public:
    int xpos;
    float ypos;
    std::shared_ptr<ReadContainer> data;
    std::shared_ptr<Exon> exon;
    ClosureCnt(LinearPlot* parent, std::shared_ptr<ReadContainer> d, int x, int y) :
      data(d),
      xpos(x),
      ypos(y),
      exon(new Exon({xpos - data->length()/2, xpos + data->length()/2, y, 0, data->chromosome}))
    {
      parent->positions.push_back(exon);
      log("-- creating at " + std::to_string(x) + "," + std::to_string(y));
    }
  };

  std::queue<ClosureCnt*> q;
  q.push(new ClosureCnt(this, seed, 0, 0.0f));
  seed->flags |= ReadContainer::PROCESSED;
  
  while (!q.empty()) {
    auto nodeCnt = q.front();
    auto node = nodeCnt->data;
    q.pop();

    auto chr = addChromosome(genome, node->chromosome);
    chr->addExon(node->threePrimeEnd, node->fivePrimeEnd);

    int xThis = nodeCnt->xpos;
    int yThis = nodeCnt->ypos;
    int rt = xThis + node->length() / 2;
    int lt = xThis - node->length() / 2;
    
    if (node->fivePrimeRead != nullptr) {
      int nPreds = node->fivePrimeRead->size();
      float y0 = -nPreds / 2;
      y0 = (nPreds % 2 == 0) ? y0 + 0.5 : y0;
      for (auto pred : *(node->fivePrimeRead)) {
	if (!(pred->flags & ReadContainer::PROCESSED)) {
	  int nThis = pred->findLink(node);
	  log("--> linknum: " + std::to_string(nThis) +
	      " of " + std::to_string(nPreds) );
	  q.push(new ClosureCnt(this, pred, lt - dx - pred->length()/2, yThis + (y0 + nThis) * dy));
	  pred->flags |= ReadContainer::PROCESSED;
	}
      }
    }

    if (node->threePrimeRead != nullptr) {
      int nSuccs = node->threePrimeRead->size();
      float  y = -nSuccs / 2;
      y = (nSuccs % 2 == 0) ? y + 0.5 : y;
      y *= dy;
      for (int i(0); i < nSuccs; ++i) {
	log("--> working on succ #" + std::to_string(i) +
	    " of " + std::to_string(nSuccs) );
	auto succ = node->threePrimeRead->at(i);
	int succX = rt + dx + succ->length() / 2;
	if (!(succ->flags & ReadContainer::PROCESSED)) {
	  ClosureCnt* succCnt = new ClosureCnt(this, succ, succX, yThis + y);
	  q.push(succCnt);
	  succ->flags |= ReadContainer::PROCESSED;
	}
	std::shared_ptr<Exon> lineEnd( new Exon({succX - succ->length()/2, succX + succ->length()/2, (int)(yThis + y), 0, succ->chromosome}) );
	LineEnds connectionLine = {nodeCnt->exon, lineEnd, node->threePrimeRefs->at(i)};
	connections.push_back(connectionLine);
	y += dy;
      }
    }
  }

  for (auto& chr : chromosomes) {
    chr.second->printout();
  }
}


std::shared_ptr<Rect> LinearPlot::boundingRect() const {
  std::shared_ptr<Rect> rect(new Rect);
  
  int x0(0), x1(0), y0(0), y1(0), x2, y2;
  for (auto& ex : positions) {            // calculate space for main plot
    if (x0 == 0 || ex->lt < x0) x0 = ex->lt;
    if (x1 == 0 || ex->rt > x1) x1 = ex->rt;
    if (y0 == 0 || ex->y < y0) y0 = ex->y;
    if (y1 == 0 || ex->y > y1) y1 = ex->y;
  }

  rect->x = x0 - 2*dx;
  rect->y = y0 - 2*dy;
  rect->w = x1 - x0 + 4*dx;
  rect->h = y1 - y0 + 4*dy;
  rect->h += 3 * chromosomes.size() * dy; // reserve space for chromosomes
  
  return rect;
}


void LinearPlot::writeEps ( const std::string& fileName ) const {
  std::ofstream out(fileName);
  assume(out.good(), "Error writing to: " + fileName, false);
  if (!out.good()) return;

  // file header + draw used chromosomes
  auto rect = writeEpsHeader(out, dx, dy);
  auto innerRect = boundingRect();

  log("-- obounds: " + std::to_string(rect->x) + "/" + std::to_string(rect->y) + ", " +
      std::to_string(rect->w) + "x" + std::to_string(rect->h) );

  log("-- ibounds: " + std::to_string(innerRect->x) + "/" + std::to_string(innerRect->y) + ", " +
      std::to_string(innerRect->w) + "x" + std::to_string(innerRect->h) );

  int xOffs(-innerRect->x);
  int yOffs((chromosomes.size() + 1) * dy - innerRect->y);

  log("-- offsets: " + std::to_string(xOffs) + "/" + std::to_string(yOffs) );

  // draw connection lines
  for (auto& conn : connections) {
    int x0(xOffs);
    int x1(x0);
    int y0(yOffs + 0.25*dy);
    int y1(y0);
    if (conn.a->lt < conn.b->lt) { // a left of b
      x0 += conn.a->rt;
      x1 += conn.b->lt;
      y0 += conn.a->y;
      y1 += conn.b->y;
    }
    else {   // a right of b, should not happen anymore
      x0 += conn.b->rt;
      x1 += conn.a->lt;
      y0 += conn.b->y;
      y1 += conn.a->y;
    }
    float dx = x1 - x0;
    float dy = y1 - y0;
    float dydx = dy/dx;
    float L = sqrt(dx*dx + dy*dy);
    float cy = y0 + 0.5 * L * dydx - 0.15 * dy;
    out << "(" << conn.num << ") " << x0 + dx/2 << " " << cy << " "; // information for label
    out << x1 << " " << y1 << " " << x0 << " " << y0 << " conn\n";   // information for line
  }

  // draw treeplot
  for (auto& exon : positions) {
    auto color = PALETTE[exon->chr];
    int width = exon->rt - exon->lt;
    int xPos = exon->lt + xOffs;
    int yPos = exon->y + yOffs;
    out << color[0] << " " << color[1] << " " << color[2] << " "; // define current color
    out << width << " " << xPos << " " << yPos << " exon\n";
  }

  out << "%%EOF";
  out.flush();
  out.close();
}
