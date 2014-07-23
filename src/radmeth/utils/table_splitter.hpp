
#ifndef TABLE_SPLITTER_HPP_
#define TABLE_SPLITTER_HPP_

#include <iostream>
#include <string>

class TableSplitter {
public:
	TableSplitter(std::istream &full_table, size_t rows_per_table);
	bool get_table(std::ostream &chunk_table);
private:
	std::string sample_names_;
	std::istream *full_table_;
	size_t rows_per_table_ ; 
};

#endif