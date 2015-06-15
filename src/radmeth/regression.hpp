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

// Declares the Regression class which implements beta-binomial regression.

#ifndef REGRESSION_HPP_
#define REGRESSION_HPP_

// Std headers.
#include <vector>
#include <string>

struct Design {
  std::vector<std::string> factor_names;
  std::vector<std::string> sample_names;
  std::vector<std::vector<double> > matrix;
};

std::istream& operator>> (std::istream &is, Design &design);
std::ostream& operator<< (std::ostream &os, const Design &design);
void remove_factor(Design &design, size_t factor);

struct SiteProportions {
  std::string chrom;
  size_t position;
  std::string strand;
  std::string context;
  std::vector<size_t> total;
  std::vector<size_t> meth;
};

std::istream&
operator>>(std::istream &table_encoding, SiteProportions &props);

struct Regression {
  Design design;
  SiteProportions props;
  double max_loglik;
};

bool fit(Regression &r,
          std::vector<double> initial_parameters = std::vector<double>());

#endif // REGRESSION_HPP_
