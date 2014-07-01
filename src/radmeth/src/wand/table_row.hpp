/*    Copyright (C) 2013 University of Southern California and
 *                       Egor Dolzhenko
 *                       Andrew D Smith
 *
 *    Authors: Andrew D. Smith and Egor Dolzhenko
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 */

// TODO: Merge all auxiliary functions into the TableRow class.

// This file contains the definition of the TableRow structure designed to store
// rows of the proportion table. There are also the following auxiliary 
// functions:
// read_row -- parses an encoding of the proportion table row into an instance 
//             of the TableRow structure.
// assert_compatibility -- throws an exception if row is not compatible with the
//                         design.
// has_low_coverage -- returns true if coverages are zero either in all samples  
//                     associated with the test factor or in all samples that 
//                     do not.

#ifndef COUNT_TABLE_HPP_
#define COUNT_TABLE_HPP_

#include <string>
#include <vector>
#include <sstream>

#include "design.hpp"

struct TableRow {
  std::string chrom;
  size_t begin;
  size_t end;
  std::vector<size_t> meth_counts;
  std::vector<size_t> total_counts;
  
  friend std::ostream& operator<<(std::ostream& os, const TableRow& row) {
    os << row.chrom << ":" << row.begin << ":" << row.end;
    
    for (size_t sample = 0; sample < row.total_counts.size(); ++sample)
      os << " " << row.total_counts[sample]
         << " " << row.meth_counts[sample];
    
    return os;
  }
};

void read_row(std::string row_encoding, TableRow &table_row);

void assert_compatibility(Design design, TableRow &row);

bool has_low_coverage(const Design &design, size_t test_factor, 
                      const TableRow &row);

bool has_extreme_counts(const Design &design, const TableRow &row);

#endif //COUNT_TABLE_HPP_
