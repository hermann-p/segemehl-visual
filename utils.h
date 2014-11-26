#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>

// Logging
void log ( const std::string msg );
void assume ( const bool isGood, const std::string msg, const bool fatal = true );

// Split a string at given delimiters
std::vector< std::string > strsplit ( const std::string input, const std::string& delim, bool keepEmpty = true );

// Mathhelpers
template <typename T>
bool isBetweenBC ( T a, T b, T c );

#endif // UTILS_H