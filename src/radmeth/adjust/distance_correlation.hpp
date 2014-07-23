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
 
#ifndef DISTANCE_CORRELATION_HPP_
#define DISTANCE_CORRELATION_HPP_

#include <vector>
#include <sstream>
#include <string>
#include <algorithm>

#include "pvallocus.hpp"

#include "bin_for_distance.hpp"

class DistanceCorrelation {
public:
  DistanceCorrelation(BinForDistance bin_for_dist) 
    : bin_for_dist_(bin_for_dist) {};
  std::vector<double> correlation_table(const std::vector<PvalLocus> &loci);

private:
  double correlation(const std::vector<double> &x, 
                      const std::vector<double> &y);
  void bin(const std::vector<PvalLocus> &loci);
  std::vector< std::vector<double> > x_pvals_for_bin_;
  std::vector< std::vector<double> > y_pvals_for_bin_;
  const BinForDistance bin_for_dist_;
};

#endif //DISTANCE_CORRELATION_HPP_
