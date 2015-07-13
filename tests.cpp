#include "utils.h"
#include "genome.h"
#include "readcontainer.h"
#include "lockfreequeue.h"

#include <memory>
#include <iostream>
#include <vector>
#include <fstream>
#include <thread>

using namespace std;

uint errors(0);
uint lineNum(0);
string testmode;
vector<string> errorList;


void set_mode ( const string prg, const string module, const string options = "" ) {
  cout << "Testing " << prg << " from " << module << ((options == "") ? "" : " "+options) << endl;
  testmode = module + "  " + prg + "  " + options;
}


void test_assert( bool success, string strTest, string strError ) {
  cout << to_string(++lineNum) << "  Testing: " << strTest << "... \t";
  if (success) {
    cout << "ok\n";
  }
  else {
    errorList.push_back(to_string(lineNum) + ": " + testmode);
    cout << "ERROR: " << strError << endl;
  }
}

void queueTest();
void threadingTest();

int main( int argc, char** argv ) {
  /* ######################################################################
   * utils.h   strsplit
   * ###################################################################### */
  set_mode("strsplit", "utils");
  {
    string testIn = "This is a test string";
    uint nExpected = 5;
    auto result = strsplit(testIn, " ");
    string strExpected[] = {"This", "is", "a", "test", "string"};
    test_assert(result.size() == nExpected, "number of splits", to_string(result.size()) + " tokens instead of " + to_string(nExpected));
    for (uint i(0); i < result.size(); ++i) {
      test_assert(strExpected[i] == result.at(i), "value #" + to_string(i), "strsplit: " + result.at(i) + " should be " + strExpected[i]);
    }
  }
  set_mode("strsplit", "utils", "without empty tokens");
  {
    string testIn = "This is  a test      string";
    uint nExpected = 5;
    auto result = strsplit(testIn, " ", false);
    string strExpected[] = {"This", "is", "a", "test", "string"};
    test_assert(result.size() == nExpected, "number of splits", to_string(result.size()) + " tokens instead of expected " + to_string(nExpected));
    for (uint i(0); i < result.size(); ++i) {
      test_assert(strExpected[i] == result.at(i), "value #" + to_string(i), "strsplit: " + result.at(i) + " should be " + strExpected[i]);
    }
  }
  set_mode("strsplit", "utils", "with empty tokens");
  {
    string testIn = "This is   a test  string";
    uint nExpected = 8;
    auto result = strsplit(testIn, " ", true);
    string strExpected[] = {"This", "is", "", "", "a", "test", "", "string"};
    test_assert(result.size() == nExpected, "number of splits", to_string(result.size()) + " tokens instead of " + to_string(nExpected));
    for (uint i(0); i < result.size(); ++i) {
      test_assert(strExpected[i] == result.at(i), "value #" + to_string(i), "strsplit: " + result.at(i) + " should be " + strExpected[i]);
    }
  }
  
  test_assert(ifstream("testfile.txt").good(), "Looking for testfile", "need any file called 'testfile.txt' in cwd, queue tests now run with 0 elements, so ignore");
  queueTest();
  threadingTest();

  /* ######################################################################
   * Done, print results
   * ###################################################################### */
  
  cout << "Unit testing done, " << errors << " error(s)\n";
  for (auto& err_str : errorList) {
    cout << err_str << endl;
  }
  return 0;
}


void queueTest() {
  set_mode("Queue", "lockfreequeue", "(sequential access)");
  ifstream infile("testfile.txt");
  uint conflicts(0), counts(0);
  LockFreeQueue<string> q;
  
  for (string line; getline(infile, line); counts++) {
    q.push(line);
  }
  q.done();
  
  infile.close();
  infile.open("testfile.txt");
  
  string testLine, originalLine;
  while (q.pop(testLine) && getline(infile, originalLine)) {
    counts--;
    if (testLine != originalLine) {
      conflicts++;
    }
  }
  
  test_assert(conflicts == 0 && counts == 0, "Basic push & pop", "Mismatch in " + to_string(conflicts) + " lines, " + to_string(counts) + " leftovers");
}

void threadingTest() {
  string testFileName("testfile.txt");
  set_mode("Queue", "lockfreequeue", "(multithread access)");
  uint conflicts(0), counts(0);
  LockFreeQueue<string> q;
  
  ifstream infile(testFileName);
  vector<string> fileCnt;
  for (string line; getline(infile, line);) {
    fileCnt.push_back(line);
  }
  infile.close();
  
  auto creatorTask = [&q,&counts,&testFileName]() {
    ifstream infile(testFileName);
    for (string line; getline(infile, line);) {
      q.push(line);
      counts++;
    }
    q.done();
    infile.close();
  };
  
  auto accessorTask = [&q,&counts,&conflicts,&fileCnt]() {
    uint n(0);
    string testLine;
    while(q.pop(testLine)) {
      if (fileCnt.at(n++) != testLine) {
	conflicts++;
      }
      counts--;
    }
  };
  
  thread creator(creatorTask);
  thread accessor(accessorTask);
  
  creator.join();
  accessor.join();
  
  test_assert(conflicts == 0 && counts == 0, "Multithreaded push & pop", "Mismatch in " + to_string(conflicts) + " lines, " + to_string(counts) + " leftovers");
}

