
#include <vector>
#include <iostream>

#include "gmock/gmock.h"

#include "pvallocus.hpp"
#include "bin_for_distance.hpp"
#include "distance_correlation.hpp"

using std::vector;
using std::istringstream;

using ::testing::Eq; using ::testing::DoubleNear;

const double max_abs_error = 0.001;

TEST(DistancewiseCorrelation, ComputeCorrelationPerBinNoEmptyBins) {
  istringstream input ( "chr1 1 2 c 0.01\n"
                        "chr1 5 6 c 0.03\n"
                        "chr1 10 11 c 0.04\n"
                        "chr1 14 15 c 0.05\n"
                        "chr1 70 71 c 0.15\n"
                        "chr1 80 81 c 0.55"  );
  
  BinForDistance bin_for_dist("1:15:5");
  
  vector<PvalLocus> pval_loci;
  initialize_pval_loci(input, pval_loci);

  DistanceCorrelation corr(bin_for_dist);  
  
  vector<double> cor_for_bin = corr.correlation_table(pval_loci);
  
  ASSERT_THAT(cor_for_bin[0], DoubleNear(0.969632, max_abs_error));
  ASSERT_THAT(cor_for_bin[1], DoubleNear(0.956296, max_abs_error));
}


TEST(DistancewiseCorrelation, ComputeCorrelationPerBinSomeEmptyBins) {
  istringstream input ( "chr1 1   2 c 0.01\n"
                        "chr1 5   6 c 0.03\n"
                        "chr1 10 11 c 0.04\n"
                        "chr1 14 15 c 0.05\n"
                        "chr1 70 71 c 0.15\n"
                        "chr1 80 81 c 0.55\n" );

  BinForDistance bin_for_dist("1:20:5");
  DistanceCorrelation corr(bin_for_dist);
 
  vector<PvalLocus> pval_loci;
  initialize_pval_loci(input, pval_loci);
  
  vector<double> cor_for_bin = corr.correlation_table(pval_loci);
  ASSERT_THAT(cor_for_bin.size(), Eq(3));
  ASSERT_THAT(cor_for_bin[0], DoubleNear(0.969632, max_abs_error));
  ASSERT_THAT(cor_for_bin[1], DoubleNear(0.956296, max_abs_error));
  ASSERT_THAT(cor_for_bin[2], DoubleNear(0, max_abs_error));
}