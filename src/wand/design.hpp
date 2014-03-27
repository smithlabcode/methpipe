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

// Contains declaration of the Design class whose objects store design 
// matrices.

#ifndef DESIGN_HPP_
#define DESIGN_HPP_

#include <sstream>
#include <vector>

// Objects of this class represent design matrices.
//
// A Design object can be created like this:
//
//  string encoding = "f1\tf2\ns1\t1\t1\ns2\t1\t0";
//  istringstream iss(encoding);
//  Design design(iss);
//
// To drop the second factor from the desin:
//
//  design.remove_factor(1);

class Design {
public:
  Design() {}
  Design(std::istream &is);
  size_t num_factors() const { return factor_names_.size(); }
  size_t num_samples() const { return sample_names_.size(); }
  std::vector<std::string> factor_names() const { return factor_names_; }
  std::vector<std::string> sample_names() const { return sample_names_; }
  std::vector<std::vector<double> > matrix() const { return matrix_; }
  double operator() (size_t sample, size_t factor) const;
  void remove_factor_name(std::string name);
  void remove_factor(size_t factor);
  bool operator== (const Design &other_design) const;
  bool operator!= (const Design &other_design) const;
  friend std::ostream& operator<<(std::ostream& os, const Design &design);
private:
  std::vector<std::string> factor_names_;
  std::vector<std::string> sample_names_;
  std::vector<std::vector<double> > matrix_;
};

#endif // DESIGN_HPP_
