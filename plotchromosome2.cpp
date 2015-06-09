#include "plotchromosome2.h"
#include "readcontainer.h"
#include "vplot2.h"
#include "utils.h"

#include <iostream>

PlotChromosome::PlotChromosome ( const size_t N, const chr_pos_t L, const std::string name )
  : len(L),
    name(name)
{
  //color = COLORS[L];
}


int PlotChromosome::nExons () const {
  return exons.size();
}


void PlotChromosome::addExon ( const std::shared_ptr<ReadContainer> exon, int layer ) {

  auto p1 = exon->fivePrimeEnd;
  if (exons.find(p1) == exons.end()) { // already existss
    exons.emplace(p1, exon);
  }
  if (!exon->moreData) {
    exon->moreData = new PlotInfo;
    ((PlotInfo*)exon->moreData)->layer = layer;
  }
  else {
  }
}


// number all displayed fragments in order of position
void PlotChromosome::assignIds() {
  uint id(0);
  for (auto &exon_it : exons) {
    PlotInfo* pi = (PlotInfo*)(exon_it.second->moreData);
    pi->id = id++;
  }
}


void PlotChromosome::printout() {
  std::cout << "chromosome " << name << std::endl;
  for (auto& ex : exons) {
    std::cout << ex.first << "-" << ex.second << std::endl;
  }
}


std::shared_ptr<Rect> PlotChromosome::boundingRect() {
  std::shared_ptr<Rect> rect(new Rect);
  *rect = {0,0,0,0}; // not used yet
  return rect;
}


void PlotChromosome::writeEps ( std::ostream& out, const Rect dim, const int dx, const int dy, const float scale, const float col[3] ) {
  assignIds();
  int cy = dim.h + .5 * dy;

  // Chromosome name
  out << dx << " " << cy - 0.25 * dy << " moveto\n";
  out << "(" << name << ") show\n";

  // Background line
  out << len*scale << " " // w
      << 2*dx << " " // x
      << cy << " cL\n"; // y

  // Rectangle for exons
  for (auto& exon : exons) {
    out << col[0] << " " << col[1] << " " << col[2] << " "  // r g b
	 << (exon.second->threePrimeEnd-exon.first) * scale << " " // w
	 << 2*dx + scale * exon.first << " "   // x
	 << cy - 0.25 * dy << " exon\n";
  }
  out << std::endl;
}
