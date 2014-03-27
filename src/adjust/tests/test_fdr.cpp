
#include <vector>

#include "gmock/gmock.h"

#include "locus.hpp"

#include "fdr.hpp"

using std::vector;

using ::testing::Eq;

TEST(Fdr, TestOnScoresThatAreValidPvalues) {
  vector<Locus> input_loci = { Locus("chr1 1 2 c 0.03 +"),
                               Locus("chr1 5 6 c 0.004 +"),
                               Locus("chr1 10 11 c 0.01 +"),
                               Locus("chr1 14 15 c 0.03 +"),
                               Locus("chr1 70 71 c 0.85 +"),
                               Locus("chr1 80 81 c 0.55 +") };
 
  vector<Locus> expected_output_loci = { Locus("chr1 1 2 c:0.03 0.045 +"),
                                         Locus("chr1 5 6 c:0.004 0.024 +"),
                                         Locus("chr1 10 11 c:0.01 0.03 +"),
                                         Locus("chr1 14 15 c:0.03 0.045 +"),
                                         Locus("chr1 70 71 c:0.85 0.85 +"),
                                         Locus("chr1 80 81 c:0.55 0.66 +") };
  
  vector<LocusIterator> good_loci_iterators;
  
  get_iterators_to_good_loci(input_loci, good_loci_iterators);
  
  fdr(good_loci_iterators);
  
  ASSERT_THAT(input_loci, Eq(expected_output_loci));
}

TEST(Fdr, TestOnScoresThatAreValidAndInvalidPvalues) {
  vector<Locus> input_loci = { Locus("chr1 1 2 c 0.03 +"),
                               Locus("chr1 5 6 c 0.004 +"),
                               Locus("chr1 10 11 c 0.01 +"),
                               Locus("chr1 14 15 c 0.03 +"),
                               Locus("chr1 17 18 c -0.1 +"),
                               Locus("chr1 70 71 c 0.85 +"),
                               Locus("chr1 72 73 c 1.05 +"),
                               Locus("chr1 80 81 c 0.55 +") };
 
  vector<Locus> expected_output_loci = { Locus("chr1 1 2 c:0.03 0.045 +"),
                                         Locus("chr1 5 6 c:0.004 0.024 +"),
                                         Locus("chr1 10 11 c:0.01 0.03 +"),
                                         Locus("chr1 14 15 c:0.03 0.045 +"),
                                         Locus("chr1 17 18 c -0.1 +"),
                                         Locus("chr1 70 71 c:0.85 0.85 +"),
                                         Locus("chr1 72 73 c 1.05 +"),
                                         Locus("chr1 80 81 c:0.55 0.66 +") };
  
  vector<LocusIterator> good_loci_iterators;
  
  get_iterators_to_good_loci(input_loci, good_loci_iterators);
  
  fdr(good_loci_iterators);
  ASSERT_THAT(input_loci, Eq(expected_output_loci));
}
