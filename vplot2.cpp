#include "vplot2.h"
#include "utils.h"
#include "readcontainer.h"

vPlot::vPlot()
{
}

int vPlot::exonCount () {
  int n = 0;
  for (auto& chr : chromosomes) {
    n += chr.second->nExons();
  }
  return n;
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


std::shared_ptr<Rect> vPlot::writeEpsHeader ( std::ostream& out, const int dx, const int dy, Rect& contentBounds ) {
  // determine longest chromosome to scale output chromosomes
  int longest = 0;
  for (auto& chr : chromosomes) {
    if (chr.second->len > longest) {
      longest = chr.second->len;
    }
  }

  out << "%!PS-Adobe-3.0 EPSF-3.0\n";

  auto br = boundingRect();
  float scale = 1.0 * br->w / contentBounds.w;
  int h = br->h + scale * contentBounds.h;
  out << "%%BoundingBox: " << 0 << " " << 0 << " " << br->w << " " << h << "\n";

  br->x += dx;
  br->y += dy;
  br->w -= 2*dx;
  br->h -= 2*dy;

  // setup font
  out << "/Helvetica findfont\n" << dy * 0.5 << " scalefont\n0 setgray\nsetfont\n\n";

  float chr_x_scale = (br->w - 2*dx) * 1.0 / longest; // scale chromosomes to fit to plot

  // ps-functions to save file size

  // line from x y straight right by w
  out << "/cL { %w x y\n";
  out << " newpath\n";
  out << " 0 setgray\n";
  out << " moveto\n";
  out << " 0 rlineto\n";
  out << " 1 setlinewidth\n";
  out << " stroke\n";
  out << "} def\n\n";

  // line from x0 y0 to x1 y1
  out << "/conn { %link_num, ln_x, ln_y, x0, y0, x1, y1\n";
  out << " newpath\n";
  out << " moveto\n";
  out << " lineto\n";
  out << " 0 setgray\n";
  out << " stroke\n";
  out << " /Helvetica findfont\n ";
  out << dy * 0.3 << " scalefont\n";
  out << " 1 0 0 setrgbcolor setfont\n";
  out << " moveto\n";
  out << " show\n";
  out << "}def\n\n";

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
  int chrColId(0);                // pick colors from palette
  for (auto& chr : chromosomes) { // write to file and adjust free space
    chr.second->writeEps(out, *br, dx, dy, chr_x_scale, PALETTE[chrColId++]);
    //    auto h = chr.second->boundingRect()->h + dy;
    int h = dy;
    //    br->y += h;
    br->h += h;
    debug("chromosome painted, bbox now: "
	  + std::to_string(br->x) + "/" + std::to_string(br->y)
	  + ": " + std::to_string(br->w) + "x" + std::to_string(br->h));
  }
  br->h += dy;

  return br;
}


void vPlot::connectExons ( p_read_t lt, p_read_t rt ) {
  assume(lt != nullptr, "Left exon doesn't exist");
  assume(rt != nullptr, "Right exon doesn't exist");
  std::cout << "---- connecting " << *lt << " to " << *rt << std::endl;
  if (!lt->threePrimeRead) {
    lt->threePrimeRead = new link_list_t;
  }
  if (!rt->fivePrimeRead) {
    rt->fivePrimeRead = new link_list_t;
  }

  lt->threePrimeRead->push_back(rt);
  rt->fivePrimeRead->push_back(lt);
}
