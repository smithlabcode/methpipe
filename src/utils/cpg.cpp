
#include <cmath>
#include <string>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <iostream>

#include "cpg.hpp"

using std::string; using std::floor;
using std::istringstream;

Cpg::Cpg(string encoding) {
  char sign;
  string name;
  double level;

  istringstream encoding_stream(encoding);	

	if (!(encoding_stream >> chrom_ >> locus_ >> sign >> name >> level >> total_))
		throw (std::logic_error("Couldn't parse a line \"" + encoding + "\".\n"));

	if ( !(0 <= level && level <= 1) )
		throw (std::logic_error("Methylation level must be between 0 and 1: "
														+ encoding + "\n"));
  
  meth_ = floor(level*total_ + 0.5);
}
