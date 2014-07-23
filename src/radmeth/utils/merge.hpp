
#ifndef MERGE_HPP_
#define MERGE_HPP_

#include <vector>
#include <sstream>

void merge_methylomes(std::vector<std::string> names, 
           						std::vector<std::istream*> methylomes, 
           						std::ostream &count_table);

#endif //MERGE_HPP_
