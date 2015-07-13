// My headers
#include "genome.h"
#include "readcontainer.h"
#include "utils.h"
#include "linplot2.h"

// System headers
#include <iostream>
#include <vector>
#include <unistd.h>
#include <string>
#include <memory>

using namespace std;

#define WIDTH 1024
#define HEIGHT 768

// Store command line options
struct Options {
  bool valid = false;
  bool circular = false;
  bool multistrand = false;
  vector< pair<string, chr_pos_t> > reads;
  string fileName = "";
  int xres = WIDTH;
  int yres = HEIGHT;
  string outFileName = "cb_out";
  string reportFileName = "";
};


// Print usage information
void printUsage (char** args) {
  cout << 
"\n"
"\n"
"Usage: \n    " << args[0] << " -c|m|p <chr:pos> -f <path-to-file> [-r WxH] [-o <[path/]filename>] [-s <[path/]filename>]\n" 
"\n"
"Optional:\n"
"    -h            print this help and exit\n"
// "    -r WxH        set output width and height in pixels (won't affect vector graphics)\n"
"    -o filename   set output path and filename\n"
"    -s filename   print summary to file\n"
"At least one required:\n"
"    -c            detect all circular transcripts (not implemented yet)\n"
"    -m            detect all multistrand spliced transcripts\n"
"    -p chr:pos    find all transcripts on chromosome chr, position pos\n"
"                  (multiple -p chr:pos allowed)\n"
"Required:\n"
"    -f filename   path to segemehl-generated .sam input file\n\n";
}


// parse command line options
void parseOptions ( int argc, char** argv, Options* options ) {
  int o;
  // parse all arguments and mark the configuration as valid if both a filename and
  // at least one optional-mandatory selection parameter were given
  while ((o = getopt(argc, argv, "hcmo:r:p:f:s:")) != -1) {
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
    case 's': // set filename for summary
      options->reportFileName = optarg;
    case 'r': {
      auto dims = strsplit(optarg, "x", false);
      if (dims.size() == 2) {
	unsigned int w = atoi(dims[0].c_str());
	unsigned int h = atoi(dims[0].c_str());
	if (w + h == 0) { // erroneous arguments
	  w = WIDTH;
	  h = HEIGHT;
	}
	else { // one argument given: set aspect ratio
	  w = (w != 0) ? w : h * 1.33;
	  h = (h != 0) ? h : w * 0.75;
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
      auto pos = (chr_pos_t)atoi(values[1].c_str());
      pair<string, chr_pos_t> readPos(values[0], pos);
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
    cout <<  "    " << chrpos.first << ":" << chrpos.second << endl;
  }

  Genome* g = new Genome(options.multistrand, options.circular);
  g->read(options.fileName);

  ostream* report = nullptr;
  if (options.reportFileName != "") {
    report = new ofstream(options.reportFileName);
    if (assume(report->good(), "Couldn't create file for report, reports will be skipped, warnings will be shown.")) {
      *report << "filename_event\t" << "chromosome_positions\t" << "min_read_depth\t" << "max_read_depth\n";
    }
  }
  
  if (options.circular) {
    cout << "Skipping circular detection: not implemented yet" << endl;
    uint N(0);
    for (auto& seed : *(g->circular)) {
      cout << "Circular seed: " << *seed << " flags: " << to_string(seed->flags) << endl;
      if (!(seed->flags & ReadContainer::PROCESSED)) {
	LinearPlot plot;
 	plot.fromRead(seed, g);
	plot.writeEps(options.outFileName + "_circ_" + to_string(++N) + ".eps");
	if (report) {
	  plot.addToSummary(*report, options.fileName + "_circular_" + to_string(N) );
	}
      }
    }
  }

  if (options.multistrand) {
    uint N(0);
    for (auto& seed : *(g->multistrand)) {
      if (!(seed->flags & ReadContainer::PROCESSED)) {
	cout << "MS seed: " << *seed << " flags: " << to_string(seed->flags) << endl;
	LinearPlot plot;
 	plot.fromRead(seed, g);
	plot.writeEps(options.outFileName + "_multi_" + to_string(++N) + ".eps");
	if (report) {
	  plot.addToSummary(*report, options.fileName + "_multi_" + to_string(N) );
	}
      }
    }
  }
  
  if (!options.reads.empty()) {
    uint N(0);
    for (auto& tpl : options.reads) {
      shared_ptr<vPlot> plot(new LinearPlot(40, 20));
      auto seed = g->getReadAt(tpl.first, tpl.second);
      if (seed != nullptr && !(seed->flags & ReadContainer::PROCESSED)) { // only existing, new plots
	LinearPlot plot;
 	plot.fromRead(seed, g);
	plot.writeEps(options.outFileName + "_pos_" + to_string(++N) + ".eps");
	if (report) {
	  plot.addToSummary(*report, options.fileName + "_manual_" + to_string(N) );
	}
      }
    }
  }

  delete g;
  return 0;
}
