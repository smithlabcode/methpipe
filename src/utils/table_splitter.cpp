
#include <iostream>
#include <vector>
#include <string>
//#include <sstream>

#include "table_splitter.hpp"

using std::istream;
using std::ostream;
using std::vector;
using std::string;
//using std::istringstream;

/*
static vector<string> 
split(string input) {
  istringstream iss(input);
  string token;
  vector<string> tokens;
  
  while (iss >> token)
    tokens.push_back(token);
  
  return tokens;
}*/

TableSplitter::TableSplitter(istream &full_table, size_t rows_per_table)
	: full_table_(&full_table), rows_per_table_(rows_per_table) {
	// Read the first line of the count table which must contain names of the
	// samples.
	//string sample_names_encoding;
	getline(*full_table_, sample_names_);

}

bool
TableSplitter::get_table(std::ostream &chunk_table) {

	chunk_table << sample_names_;

	string row;
	for (size_t row_ind = 0; row_ind < rows_per_table_; ++row_ind) {
		if( !getline(*full_table_, row) )
			return false;

		chunk_table << std::endl << row;
	}

	if( full_table_->rdbuf()->in_avail() )
		return true;

	return false;
}