#ifndef VPLOT_H
#define VPLOT_H

#include "genome.h"
#include "plotchromosome2.h"

#include <memory>
#include <fstream>

class vPlot {
 public:
  vPlot();
  
  static const int WIDTH = 1024;

  virtual PlotChromosome* addChromosome ( Genome* genome, const chr_num_t id );
  virtual void fromRead ( p_read_t seed, Genome* genom ) = 0;
  virtual void writeEps ( const std::string& fileName ) = 0;
  virtual std::shared_ptr<Rect> boundingRect () = 0;
  virtual void createPlotCoords () = 0;
  virtual void addToSummary ( std::ostream& out, std::string title ) = 0;
  
  std::shared_ptr<Rect> writeEpsHeader ( std::ostream& out, const int dx, const int dy, Rect& contentBounds );

  std::map<chr_num_t, PlotChromosome*> chromosomes;

  const float PALETTE[20][3] = { {0.886,0.898,0.043},
				 {0.243,0.243,0},
				 {0.604,0.612,0.039},
				 {0.988,0.996,0.192},
				 {0.992,1,0.322},
				 {0.561,0.835,0.039},
				 {0.149,0.227,0},
				 {0.384,0.569,0.035},
				 {0.682,0.945,0.18},
				 {0.745,0.973,0.314},
				 {0.725,0.035,0.412},
				 {0.2,0,0.11},
				 {0.494,0.031,0.286},
				 {0.847,0.161,0.537},
				 {0.925,0.298,0.643},
				 {0.416,0.059,0.6},
				 {0.11,0.008,0.165},
				 {0.286,0.043,0.408},
				 {0.541,0.169,0.733},
				 {0.682,0.31,0.875}};
				 
  uint minLinks;
  uint maxLinks;
  
 protected:
  int __N;

  int exonCount ();
  void connectExons ( p_read_t lt, p_read_t rt );
  
  // TODO: remove this hack!
  Genome* genome;
};

#endif // VPLOT_H
