#ifndef DISTANCE_CORRELATION_HPP_
#define DISTANCE_CORRELATION_HPP_

#include <vector>
#include <sstream>
#include <string>
#include <algorithm>

#include "locus.hpp"

#include "bin_for_distance.hpp"

class DistanceCorrelation {
public:
  DistanceCorrelation(BinForDistance bin_for_dist) 
    : bin_for_dist_(bin_for_dist) {};
  std::vector<double> correlation_table(const std::vector<LocusIterator> &loci);

private:
  double correlation(const std::vector<double> &x, 
                      const std::vector<double> &y);
  void bin(const std::vector<LocusIterator> &loci);
  std::vector< std::vector<double> > x_pvals_for_bin_;
  std::vector< std::vector<double> > y_pvals_for_bin_;
  const BinForDistance bin_for_dist_;
};

#endif //DISTANCE_CORRELATION_HPP_
