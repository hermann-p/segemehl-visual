#include "linplot2.h"
#include "readcontainer.h"

#include <fstream>

LinearPlot::LinearPlot(int dx, int dy) 
  : vPlot(),
    dx(dx), dy(dy),
    minx(0), maxx(0), miny(0), maxy(0),
    nextID('a')
{
}

void LinearPlot::fromRead ( std::shared_ptr<ReadContainer> seed, Genome* genome ) {
  fromRead(seed, genome,  0, 300, nullptr);
  for (auto& chr : chromosomes) {
    chr.second->printout();
  }
}

void LinearPlot::fromRead ( std::shared_ptr<ReadContainer> seed, Genome* g, const int x, const int y, std::shared_ptr<Exon> pred ) {
  seed->flags |= ReadContainer::PROCESSED;
  auto chr = addChromosome(g, seed->chromosome); // pick or create
  chr->addExon(seed->fivePrimeEnd, seed->threePrimeEnd);
  int len = seed->length();
  int x1 = x - len / 2;
  int x2 = x + len / 2;
  std::shared_ptr<Exon> exon(new Exon({x1, x2, y, nextID, seed->chromosome}));
  if (pred != nullptr) {  // not the initial seed
    LineEnds le = {exon, pred};                // parent line positions with exons
    connections.push_back(le);
  }
  positions.push_back(exon);
  nextID = (++nextID <= 'z') ? nextID : 'A'; // count up and check if in range

  if (seed->fivePrimeRead != nullptr) {      // process predecessors
    int y0 = y - (seed->fivePrimeRead->size() / 2) * dy;
    if (seed->fivePrimeRead->size() % 2 == 0) {
      y0 += dy/2;
    }
    for (auto& ex : *(seed->fivePrimeRead)) {
      if (!(ex->flags & ReadContainer::PROCESSED)) {
	int x1 =   x - dx - (len + ex->length()) / 2;
	fromRead(ex, g, x1, y0, exon);
      }
      y0 += dy;
    }
  }

  if (seed->threePrimeRead != nullptr) {     // process successors
    int y0 = y - (seed->threePrimeRead->size() / 2) * dy;
    if (seed->threePrimeRead->size() % 2 == 0) {
      y0 += dy/2;
    }
    for (auto& ex : *(seed->threePrimeRead)) {
      if (!(ex->flags & ReadContainer::PROCESSED)) {
	int x1 =   x + dx + (len + ex->length()) / 2;
	fromRead(ex, g, x1, y0, exon);
      }
      y0 += dy;
    }
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

  int xOffs(-rect->x);
  int yOffs(-(chromosomes.size() + 1) * dy);
  
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
    else {
      x0 += conn.b->rt;
      x1 += conn.a->lt;
      y0 += conn.b->y;
      y1 += conn.a->y;
    }
    out << x0 << " " << y0 << " " << x1 << " " << y1 << " conn\n";
  }

  // draw treeplot
  for (auto& exon : positions) {
    auto color = PALETTE[exon->chr];
    out << color[0] << " " << color[1] << " " << color[2] << " "; // define current color
    out << exon->rt - exon->lt << " " << exon->lt + xOffs << " " << exon->y + yOffs << " exon\n";
  }

  // don't forget text where appropriate

  out << "%%EOF";
  out.flush();
  out.close();
}
