
#include <vector>

#include "gmock/gmock.h"

#include "locus.hpp"
#include "bin_for_distance.hpp"
#include "distance_correlation.hpp"

using std::vector;

using ::testing::Eq; using ::testing::DoubleNear;

const double max_abs_error = 0.001;

TEST(DistancewiseCorrelation, ComputeCorrelationPerBinNoEmptyBins) {
  vector<Locus> input_loci; 
  input_loci.push_back(Locus("chr1 1 2 c 0.01"));
  input_loci.push_back(Locus("chr1 5 6 c 0.03"));
  input_loci.push_back(Locus("chr1 10 11 c 0.04"));
  input_loci.push_back(Locus("chr1 14 15 c 0.05"));
  input_loci.push_back(Locus("chr1 70 71 c 0.15"));
  input_loci.push_back(Locus("chr1 80 81 c 0.55"));
  
  BinForDistance bin_for_dist("1:15:5");
  DistanceCorrelation corr(bin_for_dist);
  
  vector<LocusIterator> good_loci_iterators;
  get_iterators_to_good_loci(input_loci, good_loci_iterators);
  
  
  
  vector<double> cor_for_bin = corr.correlation_table(good_loci_iterators);
  
  ASSERT_THAT(cor_for_bin[0], DoubleNear(0.982, max_abs_error));
  ASSERT_THAT(cor_for_bin[1], DoubleNear(0.993, max_abs_error));
}

TEST(DistancewiseCorrelation, ComputeCorrelationPerBinSomeEmptyBins) { 
  vector<Locus> loci; 
  loci.push_back(Locus("chr1 1 2 c 0.01"));
  loci.push_back(Locus("chr1 5 6 c 0.03"));
  loci.push_back(Locus("chr1 10 11 c 0.04"));
  loci.push_back(Locus("chr1 14 15 c 0.05"));
  loci.push_back(Locus("chr1 70 71 c 0.15"));
  loci.push_back(Locus("chr1 80 81 c 0.55"));

  BinForDistance bin_for_dist("1:20:5");
  DistanceCorrelation corr(bin_for_dist);
  
  vector<LocusIterator> good_loci_iterators;
  get_iterators_to_good_loci(loci, good_loci_iterators);
  
  vector<double> cor_for_bin = corr.correlation_table(good_loci_iterators);
  ASSERT_THAT(cor_for_bin.size(), Eq(3));
  ASSERT_THAT(cor_for_bin[0], DoubleNear(0.982, max_abs_error));
  ASSERT_THAT(cor_for_bin[1], DoubleNear(0.993, max_abs_error));
  ASSERT_THAT(cor_for_bin[2], DoubleNear(0, max_abs_error));
}


