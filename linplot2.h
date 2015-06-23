#ifndef LINPLOT_H
#define LINPLOT_H

#include "vplot2.h"
#include "utils.h"
#include "readcontainer.h"

#include <memory>
#include <vector>
#include <map>
#include <unordered_set>


class LinearPlot : public vPlot {
 public:
  LinearPlot(int dx=40, int dy=40);    // constructor with plot parameters

  void fromRead ( p_read_t seed, Genome* genome );  // construct tree from seed element
  void writeEps ( const std::string& fileName );
  void addToSummary ( std::ostream & out, std::string title );

 private:
  void assignLayers ();                // apply longest path algorithm
  //  void correctDepths ();               // correct layer assignments from multisplits
  void insertDummies ();               // neccessary for layer-spanning nodes (see Sugiyama)
  std::vector<std::vector<bool>> transitiveReduction (); // to avoid visual skew by duplicate links
  void barycenterCoords ();            // heuristic to place nodes nicely
  void createPlotCoords ();            // transform logical layer/barycenter to plottable x- and y- cordinates
  
  std::map<chr_pos_t, std::vector<p_read_t>> layeredGraph;
  std::vector<p_read_t> flatGraph;
  std::shared_ptr<Rect> boundingRect();

  Rect br;
  int dx, dy;
  //  int minx, maxx, miny, maxy;
  uint nFilter;

  char nextID;
};



#endif // LINPLOT_H
