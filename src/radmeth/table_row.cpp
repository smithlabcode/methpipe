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



#include <algorithm>

#include "smithlab_utils.hpp"
#include "table_row.hpp"

using std::string; using std::istringstream; 
using std::istream; 

static size_t 
parse_natural_number(string encoding) {
  istringstream iss(encoding);
  size_t number;
  iss >> number;
  if (!iss)
    throw SMITHLABException(encoding + " does not encode a natural number");
  return number;
}

static void
parse_row_name(string row_name_encoding, TableRow &row) {
  const size_t num_colon = 
            std::count(row_name_encoding.begin(), row_name_encoding.end(), ':');
            
  if (num_colon != 2)
    throw SMITHLABException("Each row in the count table must start with "
                            "a line chromosome:start:end. Got \"" + 
                            row_name_encoding + "\" instead" );
  
  istringstream iss(row_name_encoding);
  
  getline(iss, row.chrom, ':');
  
  if (row.chrom.empty())
    throw SMITHLABException("Error parsing " + row_name_encoding + 
                            ": chromosome name is missing.");
  
  string token;
  
  getline(iss, token, ':');
  row.begin = parse_natural_number(token);
  
  iss >> token;
  row.end = parse_natural_number(token);
}

static void
parse_row_read_counts(istringstream &row_name_encoding, TableRow &row) {
  size_t total_count, meth_count;
  
  while (row_name_encoding >> total_count >> meth_count) {
    row.total_counts.push_back(total_count);
    row.meth_counts.push_back(meth_count);
  }
  
  if (!row_name_encoding.eof())
    throw SMITHLABException("some row entries are not natural numbers: " +  
                            row_name_encoding.str());
}

void
read_row(std::string row_encoding, TableRow &row) {
  istringstream iss(row_encoding);
  
  string row_name_encoding;
  
  iss >> row_name_encoding;
  
  parse_row_name(row_name_encoding, row);
  parse_row_read_counts(iss, row);
}

void
assert_compatibility(Design design, TableRow &row) {
  const size_t num_samples = design.num_samples();
  if (num_samples != row.total_counts.size() || 
      num_samples != row.meth_counts.size())
    throw SMITHLABException("There is a row with incorrect number of samples.");
}

bool 
has_low_coverage(const Design &design, size_t test_factor, 
              const TableRow &row) {
  
  bool is_covered_in_test_factor_samples = false;
  bool is_covered_in_other_samples = false;
    
  for (size_t sample = 0; sample < design.num_samples(); ++sample) {
    if (design(sample, test_factor) == 1) {
      if (row.total_counts[sample] != 0)
        is_covered_in_test_factor_samples = true;
    } else {
      if (row.total_counts[sample] != 0)
        is_covered_in_other_samples = true;
    }
  }

  return !is_covered_in_test_factor_samples || !is_covered_in_other_samples;
}

bool 
has_extreme_counts(const Design &design, const TableRow &row) {
  
  bool is_maximally_methylated = true;
  bool is_unmethylated = true;
  
  for (size_t sample = 0; sample < design.num_samples(); ++sample) {
    if (row.total_counts[sample] != row.meth_counts[sample])
      is_maximally_methylated = false;
    
    if (row.meth_counts[sample] != 0)
      is_unmethylated = false;
  }
  
  return is_maximally_methylated || is_unmethylated;
}
