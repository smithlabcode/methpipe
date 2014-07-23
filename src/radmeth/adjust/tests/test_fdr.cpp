
#include <vector>
#include <sstream>
#include <string>

#include "gmock/gmock.h"

#include "pvallocus.hpp"
#include "fdr.hpp"

using std::vector;
using std::istringstream;
using std::ostringstream;
using ::testing::Eq;
using std::string;

TEST(Fdr, TestOnScoresThatAreValidPvalues) {
  string input = "chr1  1  2 c 0.03  +\n"
                 "chr1  5  6 c 0.004 +\n"
                 "chr1 10 11 c 0.01  +\n"
                 "chr1 14 15 c 0.03  +\n"
                 "chr1 70 71 c 0.85  +\n"
                 "chr1 80 81 c 0.55  +\n";

  istringstream stream_input (input);
  vector<PvalLocus> pval_loci;
  initialize_pval_loci(stream_input, pval_loci);
  pval_loci[0].combined_pval = 0.03;
  pval_loci[1].combined_pval = 0.004;
  pval_loci[2].combined_pval = 0.01;
  pval_loci[3].combined_pval = 0.03;
  pval_loci[4].combined_pval = 0.85;
  pval_loci[5].combined_pval = 0.55;

  fdr(pval_loci);
 
  string expected_output = "chr1\t1\t2\tc:0.03:0.03\t0.045\n"
                           "chr1\t5\t6\tc:0.004:0.004\t0.024\n"
                           "chr1\t10\t11\tc:0.01:0.01\t0.03\n"
                           "chr1\t14\t15\tc:0.03:0.03\t0.045\n"
                           "chr1\t70\t71\tc:0.85:0.85\t0.85\n"
                           "chr1\t80\t81\tc:0.55:0.55\t0.66\n";
  
  istringstream second_stream_input(input);
  ostringstream output_stream;

  update_pval_loci(second_stream_input, pval_loci, output_stream);

  ASSERT_THAT(output_stream.str(), Eq(expected_output));
}

TEST(Fdr, TestOnScoresThatAreValidAndInvalidPvalues) {
  string input = "chr1 1 2 c 0.03 +\n"
                "chr1 5 6 c 0.004 +\n"
                "chr1 10 11 c 0.01 +\n"
                "chr1 14 15 c 0.03 +\n"
                "chr1 17 18 c -0.1 +\n"
                "chr1 70 71 c 0.85 +\n"
                "chr1 72 73 c 1.05 +\n"
                "chr1 80 81 c 0.55 +\n";
 
  istringstream stream_input (input);
  vector<PvalLocus> pval_loci;
  
  initialize_pval_loci(stream_input, pval_loci);

  pval_loci[0].combined_pval = 0.03; 
  pval_loci[1].combined_pval = 0.004;
  pval_loci[2].combined_pval = 0.01;
  pval_loci[3].combined_pval = 0.03;
  pval_loci[4].combined_pval = 0.85;
  pval_loci[5].combined_pval = 0.55;

  fdr(pval_loci);

  istringstream second_stream_input(input);
  ostringstream output_stream;

  update_pval_loci(second_stream_input, pval_loci, output_stream);

 
  string expected_output = "chr1\t1\t2\tc:0.03:0.03\t0.045\n"
                           "chr1\t5\t6\tc:0.004:0.004\t0.024\n"
                           "chr1\t10\t11\tc:0.01:0.01\t0.03\n"
                           "chr1\t14\t15\tc:0.03:0.03\t0.045\n"
                           "chr1\t17\t18\tc\t-0.1\n"
                           "chr1\t70\t71\tc:0.85:0.85\t0.85\n"
                           "chr1\t72\t73\tc\t1.05\n"
                           "chr1\t80\t81\tc:0.55:0.55\t0.66\n";
  
  ASSERT_THAT(pval_loci.size(), Eq(6));
  ASSERT_THAT(output_stream.str(), Eq(expected_output));
}