#include "utils.h"
#include "genome.h"

#include <iostream>
#include <vector>
#include <string>


void log ( const std::__1::string msg ) {
  std::cout << msg << std::endl;
}


void assume ( const bool isGood, const std::__1::string msg, bool fatal ) {
  if (!isGood) {
    if (fatal) {
      std::cerr << "Error: " << msg << ", halting.";
      std::exit(-1);
    }
    else {
      std::cerr << "Warning: " << msg << std::endl;
    }
  }
}


std::vector< std::string > strsplit ( const std::__1::string input, const std::__1::string& delim, bool keepEmpty ) {
  std::string token, theStr = input;
  int L = delim.length();
  std::vector< std::string > result;
  
  while (token != theStr) {
    auto end = theStr.find_first_of(delim);
    token = theStr.substr(0, end);
    theStr = theStr.substr(end + L);
    if (keepEmpty || token.length() > 0) {
      result.push_back(token);
    }
  }
  return result;
}

/*
 * For numerical values. May be inexact for floating point numbers.
 * If one of {b,c} < a and the other one > a, then
 *       b-a <= 0 ^ c-a >= 0
 *   or  b-a >= 0 ^ c-a <= 0
 *   =>  (b-a) * (c-a) <= 0
 */
template <typename T>
bool isBetweenBC ( T a, T b, T c ) {
  return (b-a) * (c-a) <= 0;
}