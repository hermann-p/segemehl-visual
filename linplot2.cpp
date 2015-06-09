#include "linplot2.h"
#include "readcontainer.h"

#include <fstream>
#include <math.h>
#include <queue>
#include <algorithm>

LinearPlot::LinearPlot ( int dx, int dy ) 
  : vPlot(),
    dx(dx), dy(dy),
    //    minx(0), maxx(0), miny(0), maxy(0),
    nextID('a'),
    nFilter(0),
    br( {0,0,0,0} )
{
}


void LinearPlot::fromRead ( p_read_t seed, Genome* genome ) {
  if (seed->flags & ReadContainer::PROCESSED) {
    return;
  } // guardian

  std::queue<p_read_t> q;
  q.push(seed);

  while (!q.empty()) {
    auto& node = q.front();
    q.pop();
    
    if (node->flags & ReadContainer::PROCESSED) {
      continue;
    } // guardian

    auto chromosome = addChromosome(genome, node->chromosome);
    chromosome->addExon(node);
    node->flags |= ReadContainer::PROCESSED;
    node->moreData = new PlotInfo;
    node->moreData->id = flatGraph.size();
    node->moreData->c2 = 0;
    flatGraph.push_back(node);
    if (node->threePrimeRead) for (auto& tpr : *node->threePrimeRead) {
      if (!(tpr->flags & ReadContainer::PROCESSED)) {
	q.push(tpr);
      }
    } // add all three prime linked nodes
    if (node->fivePrimeRead) for (auto& fpr : *node->fivePrimeRead) {
      if (!(fpr->flags & ReadContainer::PROCESSED)) {
	q.push(fpr);
      }
    } // add all five prime linked nodes
  }
}


void LinearPlot::assignLayers () {
  std::unordered_set<uint> U;  // ids of assigned nodes
  std::unordered_set<uint> Z;  // ids of nodes assigned below current layer
  int current(1);

  // lambda to select next unassigned node with all known predecessors below current
  auto select_next_node = [&]() {
    for (auto& candidate : flatGraph) {
      if (U.find(candidate->moreData->id) == U.end()) {
	if (!candidate->fivePrimeRead) { // no predecessors => all predecessors below current
	  return candidate;
	}
	bool isgood(true);
	for (auto& pred : *candidate->fivePrimeRead) {
	  if (Z.find(pred->moreData->id) == Z.end()) { // unknown predecessors > 0
	    isgood = false;
	  }
	}
	if (isgood) {
	  return candidate;
	}
      }
    }
    return (p_read_t)nullptr; // no node found
  };

  while (U.size() < flatGraph.size()) {
    p_read_t v = select_next_node();
    if (v) { // element for this layer found
      layeredGraph[current].push_back(v);
      U.insert(v->moreData->id);
    }
    else { // no more elements for this layer
      ++current;
      Z.insert(U.begin(), U.end());
    }
  }
}

/*
void fromRead ( p_read_t seed, Genome* genome ) {
  if (seed->flags & ReadContainer::PROCESSED) { // Is already node in a graph
    return;
  }
  std::queue<p_read_t> q;
  q.push(seed);
  int layer = 0;

  // Lambda function to avoid copy-paste for bi-directional traversal
  auto traverse = [&q,layer](std::vector<p_read_t> list, int diff) {
    for (auto el : list) {
      if (el->flags & ReadContainer::PROCESSED) {
	continue;
      }
      if (!el->moreData) {
	el->moreData = new PlotInfo;
	el->moreData->layer = layer + diff;
	el->moreData->c2 = 0;
      }
      else {
	int Lnew = layer + diff;
	int Lold = el->moreData->layer;
	bool update = (diff > 1) ? (Lnew > Lold) : (Lnew < Lold); // successors (diff > 1) need to be on higher layers, predecessors on lower ones
	if (update) el->moreData->layer = Lnew;
      }
      q.push(el);
    }
  };
  
  while (!q.empty()) {
    auto node = q.front();
    q.pop();

    std::cout << "Working on: " << *node << std::endl;

    if (!(node->flags & ReadContainer::DUMMY)) {
      auto chromosome = addChromosome(genome, node->chromosome);
      chromosome->addExon(node, layer);
    }
    flatGraph.push_back(node);

    layer = node->moreData->layer;
    node->flags |= ReadContainer::PROCESSED;
    node->moreData->id = flatGraph.size();
    layeredGraph[layer].push_back(node);

    if (node->fivePrimeRead) traverse(*node->fivePrimeRead, -1);
    if (node->threePrimeRead) traverse(*node->threePrimeRead, +1);
  }

  correctDepths();
} */


/*void LinearPlot::correctDepths () {
  debug("#### CORRECTING LAYER ASSIGNMENTS ####");
  bool correct(false);
  do {
    correct = true;
    for (auto& layer : layeredGraph) {
      auto L = layer.first;

      for (int i(0); i < layer.second.size(); ++i) {
	auto& node = layer.second.at(i);
	if (node->moreData->layer != L) {
	  debug( "Found level " + std::to_string(node->moreData->layer) + " node in level " +
		 std::to_string(L) + ", -> moving");
	  correct = false;
	  layeredGraph[node->moreData->layer].push_back(node);
	  layer.second.erase( layer.second.begin() + i );
	}
      }
    }
  } while (!correct);
  }*/


void LinearPlot::insertDummies () {
  auto& lg = layeredGraph;
  for (auto& LL : lg) { // all layers in graph
    for (auto& node : LL.second) { // all nodes in layer
      if (!node->threePrimeRead || (node->flags & ReadContainer::DUMMY)) { // skip node if no successors
	continue;
      }
      int layer1 = node->moreData->layer;
      if (node->threePrimeRead) for (auto succ : *(node->threePrimeRead)) { // all successors of node
	if (!succ || (succ->flags & ReadContainer::DUMMY)) { // skip dummy nodes
	  continue;
	}
	int layer2 = succ->moreData->layer;
	int diff = layer2 - layer1;
	if (diff > 1) { // no direct neighbours, need dummy padding
	  auto lt = node; // pointer to current dummy's predecessor
	  for (int i(layer1 + 1); i <= layer2 - 1; ++i) { // all layers between nodes
	    p_read_t dummy(new ReadContainer());
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
/*
    for (auto& node : v) {
    int i(node->moreData->id);
    if (node->threePrimeRead) for (auto rt : *(node->threePrimeRead) ) {
      int j(rt->moreData->id);
      d[j][i] = d[i][j] = true;
    }
    if (node->fivePrimeRead) for (auto rt : *(node->fivePrimeRead) ) {
      int j(rt->moreData->id);
      d[j][i] = d[i][j] = true;
    }
  }
*/    
  
  std::queue<p_read_t> q;
  q.push(v[0]);

  while (!q.empty()) {
    auto node = q.front();
    q.pop();
    auto i = node->moreData->id;
    
    node->flags ^= ReadContainer::PROCESSED;
    if (node->threePrimeRead) for (auto n : *node->threePrimeRead) {
	if (n->flags & ReadContainer::PROCESSED) {
	  q.push(n);
	}
        auto j = n->moreData->id;
	d[i][j] = d[j][i] = true;
    }
    if (node->fivePrimeRead) for (auto n : *node->fivePrimeRead) {
	if (n->flags & ReadContainer::PROCESSED) {
	  q.push(n);
	}
        auto j = n->moreData->id;
	d[i][j] = d[j][i] = true;
    }
  }

  
  #ifdef DEBUG
  std::cout << "before reduction:\n";
  for (auto line : d) {
    for (auto row : line) {
      std::cout << row << "  ";
    }
    std::cout << std::endl;
  }
  #endif

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

  
  #ifdef DEBUG
  std::cout << "after reduction:\n";
  for (auto line : d) {
    for (auto row : line) {
      std::cout << row << "  ";
    }
    std::cout << std::endl;
  }
  #endif

  return d;
}


void LinearPlot::barycenterCoords () {
  auto& lg = layeredGraph;
  auto d = transitiveReduction(); // calculate reduced link matrix

  // function to order elements in a layer according to their y-coordinate (PlotInfo->c2)
  auto yComparator = [](p_read_t a, p_read_t b) {
    return a->moreData->c2 < b->moreData->c2; // weak ordering of y-coordinate c2
  };
  
  int size0 = lg.begin()->second.size();  // number of elements on leftmost position
  for (auto& LL : lg) { // for each layer
    for (auto el : LL.second) { // each node in layer
      // calculate barycenter of predecessors
      el->moreData->c2 = 0;
      if (el->fivePrimeRead) {
	uint count(1);
	for (auto pre : *(el->fivePrimeRead)) { // all predecessors
	  if (d[el->moreData->id][pre->moreData->id]) {
	    el->moreData->c2 += pre->moreData->c2;
	    ++count;
	  }
	}
	std::cout << "node " << *el << " has " << count << " predecessors.\n";
	if (count > 0) {
	  el->moreData->c2 = el->moreData->c2 / count;
	}
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
  debug("createPlotCoords() -- graph has " + std::to_string(layeredGraph.size()) + " layers");
  if (layeredGraph.empty()) {   // hierarchy was not constructed yet
    insertDummies();
    barycenterCoords();
  }

  // find size required per layer; calculate bounding box during progress
  std::map<chr_pos_t, int> llargest; // largest read per layer
  for (auto ii(layeredGraph.begin()); ii != layeredGraph.end(); ++ii) {
    auto i(ii->first);
    llargest[i] = 0;
    for (auto& node : layeredGraph.at(i)) {
      if (!(node->flags & ReadContainer::DUMMY)) {
	int size = node->length();
	llargest[i] = (size > llargest[i]) ? size : llargest[i];
      }
    }
    br.w += llargest[i];
    br.h = (layeredGraph.at(i).size() > br.h) ? layeredGraph.at(i).size() : br.h;
  }
  br.w += (layeredGraph.size() - 1) * dx;
  br.h = (br.h + (br.h - 1) / 2) * dy;

  debug("-- createPlotCoords: creating coords");
  debug("dx = " + std::to_string(dx) + ", dy = " + std::to_string(dy));
  int x0(0);
  uint count(0);
  for (auto ii(layeredGraph.begin()); ii != layeredGraph.end(); ++ii) {
    auto i(ii->first);
    if (i > layeredGraph.begin()->first) { // not the first element
      std::cout << "-- comparing " << i << " and " << (i-1) << std::endl;
      x0 += llargest.at(i-1) + dx;
    }
    for (auto& node : layeredGraph.at(i)) {
      node->moreData->c1 = x0;
      node->moreData->c3 = x0 + node->length();

      node->moreData->c2 *= (1.5 * dy);
      if (node->moreData->c2 < br.y) {
	br.y = node->moreData->c2;
      }

      std::cout << "layer " << i << ", node " << *node << " -- y: " << std::to_string(node->moreData->c2) << ", x: " << node->moreData->c1 << "-" << node->moreData->c3 << std::endl;
    }
    ++count;
  }
}


std::shared_ptr<Rect> LinearPlot::boundingRect() {
  if (br.h == 0 && br.w == 0) { // bounding box not calculated yet
    createPlotCoords();
  }
  return std::make_shared<Rect>(br);
}



void LinearPlot::writeEps ( const std::string& fileName ) {
  assignLayers();
  insertDummies();     // neccessary for correct placement
  barycenterCoords();
  createPlotCoords();

  std::ofstream out(fileName);
  assume(out.good(), "Error writing to: " + fileName, false);
  if (!out.good()) return;

  // file header + draw used chromosomes
  auto innerRect = boundingRect();
  auto chrRect  = writeEpsHeader(out, dx, dy, *innerRect);

  float s = 1.0 * (chrRect->w - 2 * dx) / innerRect->w;
  debug("Scaling factor " + std::to_string(s));

  int xOffs(dx);
  int yOffs(chrRect->y + dy - innerRect->y);

  debug("Offsets: x=" + std::to_string(xOffs) + ", y=" + std::to_string(yOffs) );

  for (auto& layer : layeredGraph) {
    for (auto& node : layer.second) {
      if (!(node->flags & ReadContainer::DUMMY)) { // draw non-dummy nodes
	auto color = PALETTE[node->chromosome];
	int width = node->moreData->c3 - node->moreData->c1;
	int xpos = node->moreData->c1;
	int ypos = node->moreData->c2;
	out << color[0] << " " << color[1] << " " << color[2] << " ";
	out << s*width << " " << xOffs + s*xpos << " " << yOffs + s*ypos << " exon\n";
	node->flags |= ReadContainer::PROCESSED;

	if (node->threePrimeRead) for (int i(0); i < node->threePrimeRead->size(); ++i) { // draw connection lines
          auto& succ = node->threePrimeRead->at(i);
	  if (succ && !(succ->flags & ReadContainer::DUMMY)) {
	    int y2 = succ->moreData->c2;
	    int x2 = succ->moreData->c1;
	    int cx = (xpos + width + x2) / 2;
	    int cy = (ypos + y2) / 2;
	    out << "(" << node->threePrimeRefs->at(i) << ") " << xOffs + s*cx << " " << yOffs + s*cy + 0.15 * dy << " ";
	    out << xOffs + s*(xpos+width) << " " << yOffs+s*ypos << " " << xOffs + s*x2 << " " << yOffs + s*y2 << " conn\n";
	  }
	}
      }
    }
  }
  
  out << "%%EOF";
  out.flush();
  out.close();
}
