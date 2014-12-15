#include "utils.h"
#include "genome.h"
#include "readcontainer.h"

#include <memory>

using namespace std;

int main(int argc, char** argv) {
  shared_ptr<Genome> g(new Genome());
  string fileName;
  if (argc > 1) {
    fileName = argv[1];
  }
  else {
    fileName = "out.sam";
  }
  g->read(fileName);
  return 0;
}
