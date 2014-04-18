
#include <vector>
#include <sstream>
#include <string>
#include <iostream>

#include "gmock/gmock.h"

#include "pvallocus.hpp"
#include "bin_for_distance.hpp"
#include "distance_correlation.hpp"

#include "combine_pvals.hpp"

using ::testing::Eq;
using std::vector;
using std::istringstream;
using std::string;
using std::istringstream;
using std::ostringstream;

TEST(CombinePvalues, FromSameChromosome) {

  string input = "chr1 1 2 c 0.01\n"
                 "chr1 5 6 c 0.03\n"
                 "chr1 10 11 c 0.04\n"
                 "chr1 14 15 c 0.05\n"
                 "chr1 70 71 c 0.15\n"
                 "chr1 80 81 c 0.55\n";

  istringstream loci_encoding(input);

  vector<PvalLocus> pval_loci;
  initialize_pval_loci(loci_encoding, pval_loci);

  BinForDistance bin_for_dist("1:15:5");
  DistanceCorrelation corr(bin_for_dist);
  
  combine_pvals(pval_loci, bin_for_dist);

  istringstream second_loci_encoding(input);
  ostringstream output_loci_encoding;

  update_pval_loci(second_loci_encoding, pval_loci, output_loci_encoding);  
  
  string exptected_output = "chr1\t1\t2\tc:0.01:0.0197799\t0\n"
                            "chr1\t5\t6\tc:0.03:0.0197799\t0\n"
                            "chr1\t10\t11\tc:0.04:0.0197799\t0\n"
                            "chr1\t14\t15\tc:0.05:0.0197799\t0\n"
                            "chr1\t70\t71\tc:0.15:0.322599\t0\n"
                            "chr1\t80\t81\tc:0.55:0.322599\t0\n";
    
  ASSERT_THAT(output_loci_encoding.str(), Eq(exptected_output));
}
