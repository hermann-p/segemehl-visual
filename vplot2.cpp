#include "vplot2.h"
#include "utils.h"

vPlot::vPlot() {
}

PlotChromosome* vPlot::addChromosome ( const Genome* genome, const chr_num_t id ) {
  auto lookup = chromosomes.find(id);
  if (lookup != chromosomes.end()) { // exists already
    return lookup->second;
  }
  PlotChromosome* chr = new PlotChromosome(chromosomes.size(),
					   genome->getLength(id),
					   genome->getChrName(id));
  chromosomes.emplace(id, chr);
  return chr;
}


std::shared_ptr<Rect> vPlot::writeEpsHeader ( std::ostream& out, const int dx, const int dy ) const {
  // determine longest chromosome to scale output chromosomes
  int longest = 0;
  for (auto& chr : chromosomes) {
    if (chr.second->len > longest) {
      longest = chr.second->len;
    }
  }

  out << "%!PS-Adobe-3.0 EPSF-3.0\n";

  auto br = boundingRect();
  out << "%%BoundingBox: " << 0 << " " << 0 << " " << br->w << " " << br->h << "\n";

  br->x += dx;
  br->y += dy;
  br->w -= 2*dx;
  br->h -= 2*dy;

  float chr_x_scale = (br->w - 2*dx) * 1.0 / longest; // scale chromosomes to fit to plot

  // ps-functions to save file size

  // line from x y straight right by w
  out << "/cL { %w x y\n newpath\n 0 setgray\n moveto\n 0 rlineto\n 1 setlinewidth\n stroke\n} def\n\n";

  // line from x0 y0 to x1 y1
  out << "/conn {\n newpath\n moveto\n lineto\n 0 setgray\n stroke\n} def\n\n";

  // draw exon of color r g b with width w on position x y
  out << "/dY " << dy*0.4 << " def\n\n";
  out << "/exon{ % r, g, b, w, x, y\n";
  out << " newpath\n";
  out << " moveto % x, y\n";
  out << " 0 dY rlineto\n";
  out << " 0 rlineto % w\n";
  out << " 0 dY neg rlineto\n closepath\n gsave\n";
  out << " setrgbcolor fill\n % r g b\n";
  out << " grestore \n 0 setgray\n stroke \n} def\n\n";
  // write chromosomes 
  for (auto& chr : chromosomes) { // write to file and adjust free space
    chr.second->writeEps(out, *br, dx, dy, chr_x_scale, PALETTE[chr.first]);
    //    auto h = chr.second->boundingRect()->h + dy;
    int h = dy;
    br->y += h;
    br->h -= h;
  }
  br->h += dy;

  return br;
}
