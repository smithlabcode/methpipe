
#include <vector>
#include <sstream>

#include "gmock/gmock.h"

#include "pvallocus.hpp"

#include "proximal_loci.hpp"

using std::vector;
using ::testing::Eq;
using std::istringstream;

TEST(ProximalCpgs, FindsProximalLociOnSingleChromosome) {
  istringstream input ( "chr1 1 2 c 0.01\n"
                        "chr1 5 6 c 0.03\n"
                        "chr1 10 11 c 0.04\n"
                        "chr1 14 15 c 0.05\n"
                        "chr1 70 71 c 0.15\n"
                        "chr1 80 81 c 0.55" );

  vector<PvalLocus> pval_loci;
  initialize_pval_loci(input, pval_loci);

  ProximalLoci proximal_loci(pval_loci, 10);
  
  vector<PvalLocus> neighbors;
  
  proximal_loci.get(neighbors);
//  vector<LocusIterator> true_neighbors_a(good_loci_iterators.begin(),
//                                         good_loci_iterators.begin() + 3);
  ASSERT_THAT(neighbors.size(), Eq(3));

  proximal_loci.get(neighbors);
//  vector<PvalLocus> true_neighbors_b( good_loci_iterators.begin(), 
//                                      good_loci_iterators.begin() + 4);
  
  ASSERT_THAT(neighbors.size(), Eq(4));
  
  for (size_t iter = 0; iter < 4; ++iter)
    proximal_loci.get(neighbors);

//  vector<LocusIterator> true_neighbors_c(good_loci_iterators.end() - 2, 
//                                         good_loci_iterators.end());
  ASSERT_THAT(neighbors.size(), Eq(2));
}


TEST(ProximalCpgs, FindsProximalLociOnMultipleChromosome) {
  istringstream input ( "chr2 3 4 c 0.03\n"
                        "chr3 5 6 c 0.04\n"
                        "chr3 7 8 c 0.05\n"
                        "chr3 9 10 c 0.15\n"
                        "chr4 11 12 c 0.55" );
  
  vector<PvalLocus> pval_loci;
  initialize_pval_loci(input, pval_loci);

  ProximalLoci proximal_loci(pval_loci, 10);
  vector<PvalLocus> neighbors;
  
  proximal_loci.get(neighbors);
//  vector<LocusIterator> true_neighbors_a(good_loci_iterators.begin(), 
//                                        good_loci_iterators.begin() + 1);
  ASSERT_THAT(neighbors.size(), Eq(1));
    
  proximal_loci.get(neighbors);
//  vector<LocusIterator> true_neighbors_b(good_loci_iterators.begin() + 1, 
//                                         good_loci_iterators.begin() + 4);
  ASSERT_THAT(neighbors.size(), Eq(3));

  proximal_loci.get(neighbors);
  proximal_loci.get(neighbors);
  ASSERT_THAT(neighbors.size(), Eq(3));

  proximal_loci.get(neighbors);
//  vector<LocusIterator> true_neighbors_c(good_loci_iterators.begin() + 4, 
//                                          good_loci_iterators.end());
  ASSERT_THAT(neighbors.size(), Eq(1));
}
