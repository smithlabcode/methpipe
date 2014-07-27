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
 
#include "design.hpp"

#include "smithlab_utils.hpp" //for SMITHLAB_exception

using std::istream; using std::string;
using std::vector; using std::istringstream;

Design::Design(istream &is) {
  string header_encoding;
  getline(is, header_encoding);

  istringstream header_is(header_encoding);
  string header_name;
  while (header_is >> header_name)
    factor_names_.push_back(header_name);
  
  string row;
  while (getline(is, row)) {

    if (row.empty())
      continue;

    istringstream row_is(row);
    string token;
    row_is >> token;
    sample_names_.push_back(token);
    
    vector<double> matrix_row;
    while (row_is >> token) {
      if (token.length() == 1 && (token == "0" || token == "1"))
        matrix_row.push_back(token == "1");
      else
        throw SMITHLABException("only binary factor levels are allowed:\n" 
                                + row);
    }

    if (matrix_row.size() != num_factors())
      throw SMITHLABException("each row must have as many columns as "
			      "factors:\n" + row);

    matrix_.push_back(vector<double>());
    swap(matrix_.back(), matrix_row);
  }
}

double
Design::operator() (size_t sample, size_t factor) const {
  return matrix_[sample][factor];
}

std::ostream&
operator<<(std::ostream& os, const Design &design) {

  for(size_t factor = 0; factor < design.num_factors(); ++factor) {
    os << design.factor_names_[factor];
    if (factor + 1 != design.num_factors())
      os << "\t";
  }
  os << std::endl;

  for(size_t sample = 0; sample < design.num_samples(); ++sample) {
    os << design.sample_names_[sample] << "\t";
    for(size_t factor = 0; factor < design.num_factors(); ++factor) {
      os << design.matrix_[sample][factor];
      if (factor + 1 != design.num_factors())
        os << "\t";
    }
      os << "\n";
  }
  return os;
}

void 
Design::remove_factor(size_t factor) {
  factor_names_.erase(factor_names_.begin() + factor);
  
  for (size_t sample = 0; sample < num_samples(); ++sample)
    matrix_[sample].erase(matrix_[sample].begin() + factor);
}

void 
Design::remove_factor_name(std::string name) {
  vector<string>::const_iterator name_it = 
                std::find(factor_names_.begin(), factor_names_.end(), name);
                
  if (name_it == factor_names_.end())
    throw SMITHLABException(name + " is not one of the factor names");

  size_t factor =  name_it - factor_names_.begin();

  remove_factor(factor);
}

bool 
Design::operator== (const Design &other_design) const {
  return (factor_names_ == other_design.factor_names_) &&
          (sample_names_ == other_design.sample_names_) &&
          (matrix_ == other_design.matrix_);
}

bool 
Design::operator!= (const Design &other_design) const {
  return !(*this == other_design);
}
