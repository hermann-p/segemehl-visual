#include "plotchromosome2.h"
#include "vplot2.h"
#include "utils.h"

#include <iostream>

PlotChromosome::PlotChromosome ( const size_t N, const chr_pos_t L, const std::string name )
  : len(L),
    name(name)
{
  //color = COLORS[L];
}


void PlotChromosome::addExon ( chr_pos_t p1, chr_pos_t p2 ) {
  if (p1 > p2) { // sort for display
    auto tmp = p1;
    p1 = p2;
    p2 = tmp;
  }

  auto lookup = exons.find(p1);
  if (lookup != exons.end()) { // already exists
    return;
  }
  exons.emplace(p1, p2);
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


void PlotChromosome::writeEps ( std::ostream& out, const Rect dim, const int dx, const int dy, const float scale, const float col[3] ) const {
  int cy = dim.y + 1.5 * dy;

  // Background line
  out << "\n" << len*scale << " " // w
      << 2*dx << " " // x
      << cy << " cL\n"; // y

  // Rectangle for exons
  for (auto& exon : exons) {
    out << col[0] << " " << col[1] << " " << col[2] << " "  // r g b
	 << (exon.second-exon.first) * scale << " " // w
	 << 2*dx + scale * exon.first << " "   // x
	 << cy - 0.25 * dy << " exon\n";
  }
}
