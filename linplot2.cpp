#include "linplot2.h"
#include "readcontainer.h"

#include <fstream>
#include <math.h>
#include <queue>
#include <algorithm>

LinearPlot::LinearPlot ( int dx, int dy ) 
  : vPlot(),
    dx(dx), dy(dy),
    minx(0), maxx(0), miny(0), maxy(0),
    nextID('a'),
    nFilter(0)
{
}


void LinearPlot::fromRead ( std::shared_ptr<ReadContainer> seed, Genome* genome ) {
  if (seed->flags & ReadContainer::PROCESSED) { // Is already node in a graph
    return;
  }
  std::queue<std::shared_ptr<ReadContainer>> q;
  q.push(seed);
  int layer = 0;

  // Lambda function to avoid copy-paste for bi-directional traversal
  auto traverse = [&q,layer](std::vector<std::shared_ptr<ReadContainer>> list, int diff) {
    for (auto el : list) {
      if (el->flags & ReadContainer::PROCESSED) {
	continue;
      }
      if (!el->moreData) {
	el->moreData = new PlotInfo;
	el->moreData->layer = layer + diff;
      }
      else {
	int Lnew = layer + diff;
	int Lold = el->moreData->layer;
	bool update = (diff > 1) ? Lnew > Lold : Lnew < Lold; // successors (diff > 1) need to be on higher layers, predecessors on lower ones
	if (update) el->moreData->layer = Lnew;
      }
      q.push(el);
    }
  };
  
  while (!q.empty()) {
    auto node = q.front();
    q.pop();

    auto chromosome = addChromosome(genome, node->chromosome);
    chromosome->addExon(node, layer);
    flatGraph.push_back(node);

    layer = node->moreData->layer;
    node->flags |= ReadContainer::PROCESSED;
    node->moreData->id = flatGraph.size();
    layeredGraph[layer].push_back(node);

    if (node->fivePrimeRead) traverse(*node->fivePrimeRead, -1);
    if (node->threePrimeRead) traverse(*node->threePrimeRead, +1);
  }
}


void LinearPlot::insertDummies () {
  auto& lg = layeredGraph;
  for (auto& LL : lg) { // all layers in graph
    for (auto& node : LL.second) { // all nodes in layer
      if (!node->threePrimeRead || (node->flags & ReadContainer::DUMMY)) { // skip node if no successors
	continue;
      }
      int layer1 = node->moreData->layer;
      for (auto succ : *(node->threePrimeRead)) { // all successors of node
	if (succ->flags & ReadContainer::DUMMY) { // skip dummy nodes
	  continue;
	}
	int layer2 = succ->moreData->layer;
	int diff = layer2 - layer1;
	if (diff > 1) { // no direct neighbours, need dummy padding
	  auto lt = node; // pointer to current dummy's predecessor
	  for (int i(layer1 + 1); i <= layer2 - 1; ++i) { // all layers between nodes
	    std::shared_ptr<ReadContainer> dummy(new ReadContainer());
	    dummy->flags |= ReadContainer::DUMMY;
	    lg.at(i).push_back(dummy);
	    dummy->moreData = new PlotInfo;
	    dummy->moreData->id = flatGraph.size();
	    flatGraph.push_back(dummy);
	    connectExons(lt, dummy);
	    lt = dummy;
	  }
	  connectExons(lt, succ);
	} // padding detected
      } // successors of node
    } // nodes in layer
  } // layers in graph
}


std::vector<std::vector<bool>> LinearPlot::transitiveReduction () {
  int N = flatGraph.size();
  auto& v = flatGraph;
  std::vector<std::vector<bool>> d(N, std::vector<bool>(N)); // connectivity array

  // step 1: mark all existing connections
  for (auto& node : v) {
    int i(node->moreData->id);
    if (node->threePrimeRead) for (auto rt : *(node->threePrimeRead) ) {
      int j(rt->moreData->id);
      d[j][i] = d[i][j] = true;
    }
  }

  // step 2: transitive reduction
  for (auto& x : v) {
    auto i = x->moreData->id;
    for (auto& y : v) {
      auto j = y->moreData->id;
      for (auto& z : v) {
	auto k = z->moreData->id;
	if (d[i][j] && d[j][k]) { // transitive link i->j->k exists
	  d[k][i] = d[i][k] = false; // remove direct i->k one
	}
      }
    }
  }

  return d;
}


void LinearPlot::barycenterCoords () {
  auto& lg = layeredGraph;
  auto d = transitiveReduction(); // calculate reduced link matrix

  for (auto line : d) {
    for (auto row : line) {
      std::cout << row << "  ";
    }
    std::cout << std::endl;
  }

  // function to order elements in a layer according to their y-coordinate (PlotInfo->c2)
  auto yComparator = [](std::shared_ptr<ReadContainer> a, std::shared_ptr<ReadContainer> b) {
    return a->moreData->c2 < b->moreData->c2; // weak ordering of y-coordinate c2
  };
  
  int size0 = lg.begin()->second.size();  // number of elements on leftmost position
  for (auto& LL : lg) { // for each layer
    for (auto el : LL.second) { // each node in layer
      // calculate barycenter of predecessors
      el->moreData->c2 = 0.0f;
      if (el->fivePrimeRead) {
	for (auto pre : *(el->fivePrimeRead)) { // all predecessors
	  el->moreData->c2 += pre->moreData->c2;
	}
	el->moreData->c2 /= el->fivePrimeRead->size();
      }
    }

    // sort elements
    sort(LL.second.begin(), LL.second.end(), yComparator);

    // get offset to center current layer elements around first layer
    int offset = (LL.second.size() > size0) ? size0 - LL.second.size() : 0;

    // apply offset to all elements in layer
    for (int i(1); i < LL.second.size(); ++i) {
      if (LL.second.at(i-1)->moreData->c2 == LL.second.at(i)->moreData->c2) {
	offset += 2; // move elements on same y position apart
      }
      LL.second.at(i)->moreData->c2 += offset;
    }
  }
}


void LinearPlot::createPlotCoords () {
}


std::shared_ptr<Rect> LinearPlot::boundingRect() const {
  return std::make_shared<Rect>(boundingBox);
}


void LinearPlot::writeEps ( const std::string& fileName ) {
  insertDummies();     // neccessary for correct placement
  barycenterCoords();  

/*
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
  out.close(); */
}
