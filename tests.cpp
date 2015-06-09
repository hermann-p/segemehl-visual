#include "utils.h"
#include "genome.h"
#include "readcontainer.h"

#include <memory>
#include <iostream>
#include <vector>

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
  cout << to_string(++lineNum) << "  Testing: " << strTest << "... ";
  if (success) {
    cout << "ok\n";
  }
  else {
    errorList.push_back(to_string(lineNum) + ": " + testmode);
    cout << "ERROR: " << strError << endl;
  }
}

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

  /* ######################################################################
   * Done, print results
   * ###################################################################### */
  
  cout << "Unit testing done, " << errors << " error(s)\n";
  for (auto& err_str : errorList) {
    cout << err_str << endl;
  }
  return 0;
}
