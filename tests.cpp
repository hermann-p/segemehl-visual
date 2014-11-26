#include "utils.h"
#include "genome.h"
#include "readcontainer.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <map>

using namespace std;

int main() {

  Genome g;
  ifstream file("out.sam");
  file >> g;
  g.printout();
  return 0;
}