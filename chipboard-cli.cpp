// My headers
#include "genome.h"
#include "readcontainer.h"
#include "utils.h"
#include "linplot.h"

// System headers
#include <iostream>
#include <vector>
#include <unistd.h>
#include <string>
#include <memory>

using namespace std;


// Store command line options
struct Options {
  bool valid = false;
  bool circular = false;
  bool multistrand = false;
  vector< pair<chr_num_t, chr_pos_t> > reads;
  string fileName = "";
  int xres = 800;
  int yres = 600;
  string outFileName = "cb_out.pdf";
};


// Print usage information
void printUsage (char** args) {
  cout << "\n\nUsage: \n    " << args[0] << " -c|m|p chr:pos -f <path-to-file> [-r WxH] [-o <filename>]" << endl;
  cout << "\n\
Arguments:\n\
    -h            print this help and exit\n\
    -r WxH        set output width and height in pixels (won't affect vector graphics)\n\
    -o filename   set output file name\n\
At least one required:\n\
    -c            detect all circular transcripts (not implemented yet)\n\
    -m            detect all multistrand spliced transcripts\n\
    -p chr:pos    find all transcripts on chromosome chr, position pos\n\
                  (multiple -p chr:pos allowed)\n\
Required:\n\
    -f filename   path to segemehl-generated .sam file\n\n" << endl;
}


// parse command line options
void parseOptions ( int argc, char** argv, Options* options ) {
  int o;
  // parse all arguments and mark the configuration as valid if both a filename and
  // at least one optional-mandatory selection parameter were given
  while ((o = getopt(argc, argv, "hcmo:r:p:f:")) != -1) {
    switch (o) {
    case 'h': // invalid options => show help and exit
      options->valid = false;
      return;
    case 'c': // circular not implemented yet
      options->circular = true;
      break;
    case 'm': // select strand-switching events
      options->multistrand = true;
      break;
    case 'f': // set filename
      options->fileName = optarg;
      break;
    case 'o':
      options->outFileName = optarg;
      break;
    case 'r': {
      auto dims = strsplit(optarg, "x", false);
      if (dims.size() == 2) {
	unsigned int w = atoi(dims[0].c_str());
	unsigned int h = atoi(dims[0].c_str());
	if (w + h == 0) { // erroneous arguments
	  w = 800;
	  h = 600;
	}
	else {
	  w = (w != 0) ? w : h * 4 / 3;
	  h = (h != 0) ? h : w * 3 / 4;
	}
      }
      break;
    }
    case 'p':
      auto values = strsplit(optarg, ":", false);
      if (values.size() != 2) {
	options->valid = false;
	return;
      }
      auto chr = (chr_num_t)atoi(values[0].c_str());
      auto pos = (chr_pos_t)atoi(values[1].c_str());
      pair<chr_num_t, chr_pos_t> readPos(chr, pos);
      options->reads.push_back(readPos);
      if (options->fileName != "") {
	options->valid = true;
      }
      break;
    }
  }
  options->valid = (options->fileName != "" &&
		    (options->circular || options->multistrand || !options->reads.empty())
		    );
}


int main ( int argc, char** argv ) {
  if (argc < 2) {
    printUsage(argv);
    return 1;
  }
  Options options;
  parseOptions(argc, argv, &options);
  if (!options.valid) {
    printUsage(argv);
    return 1;
  }
  cout << "Options successfully set" << endl;
  for (auto& chrpos : options.reads) {
    cout <<  "    " << to_string(chrpos.first) << ":" << chrpos.second << endl;
  }

  Genome* g = new Genome(options.multistrand, options.circular);
  g->read(options.fileName);
  auto* plots(new vector< shared_ptr<vPlot> >);
  
  if (options.circular) {
    cout << "Skipping circular detection: not implemented yet" << endl;
    return 1;
  }

  if (options.multistrand) {
    for (auto& seed : *(g->multistrand)) {
      if (seed->flags & ReadContainer::MULTISTRAND) {
	shared_ptr<vPlot> plot(new LinearPlot());
	plot->fromRead(seed, g);
	plots->push_back(plot);
      }
    }
  }

  if (!options.reads.empty()) {
    for (auto& tpl : options.reads) {
      try {
	shared_ptr<vPlot> plot(new LinearPlot());
	plot->fromRead(g->getReadAt(tpl.first, tpl.second), g);
	plots->push_back(plot);
      } catch (std::exception e) {
      }
    }
  }

  delete g;
  return 0;
}
