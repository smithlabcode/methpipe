
#include <vector>

#include "gmock/gmock.h"

#include "locus.hpp"
#include "bin_for_distance.hpp"
#include "distance_correlation.hpp"

#include "combine_pvals.hpp"

using ::testing::Eq;

using std::vector;

TEST(CombinePvalues, FromSameChromosome) {
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
  
  combine_pvals(good_loci_iterators, bin_for_dist);
  
  vector<Locus> exptected_output_loci;
  exptected_output_loci.push_back(Locus("chr1  1  2 c:0.01 7.19577e-05"));
  exptected_output_loci.push_back(Locus("chr1  5  6 c:0.03 7.19577e-05"));
  exptected_output_loci.push_back(Locus("chr1 10 11 c:0.04 7.19577e-05"));
  exptected_output_loci.push_back(Locus("chr1 14 15 c:0.05 7.19577e-05"));
  exptected_output_loci.push_back(Locus("chr1 70 71 c:0.15 1.00"));
  exptected_output_loci.push_back(Locus("chr1 80 81 c:0.55 1.00"));
    
  ASSERT_THAT(input_loci, Eq(exptected_output_loci));
}
