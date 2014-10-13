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

#ifndef COMBINE_PVALS_HPP_
#define COMBINE_PVALS_HPP_

#include <vector>

#include "pvallocus.hpp"
#include "bin_for_distance.hpp"

class ProximalLoci {
public:
  ProximalLoci(std::vector<PvalLocus> &loci, size_t max_distance)
    : loci_(loci), max_distance_(max_distance), next_pos_(loci.begin()) {};
  bool get(std::vector<PvalLocus> &neighbors);
  PvalLocus cur_region() {return *(next_pos_ - 1);}

private:
  const std::vector<PvalLocus> &loci_;
  size_t max_distance_;
  std::vector<PvalLocus>::const_iterator next_pos_;
};

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

void
combine_pvals(std::vector<PvalLocus> &loci,
                   const BinForDistance &bin_for_distance);

#endif //COMBINE_PVALS_HPP_
