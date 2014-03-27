
#include <vector>

#include "gmock/gmock.h"

#include "locus.hpp"

#include "proximal_loci.hpp"

using std::vector;

using ::testing::Eq;

TEST(ProximalCpgs, FindsProximalLociOnSingleChromosome) {
  vector<Locus> loci;
  loci.push_back(Locus("chr1 1 2 c 0.01"));
  loci.push_back(Locus("chr1 5 6 c 0.03"));
  loci.push_back(Locus("chr1 10 11 c 0.04"));
  loci.push_back(Locus("chr1 14 15 c 0.05"));
  loci.push_back(Locus("chr1 70 71 c 0.15"));
  loci.push_back(Locus("chr1 80 81 c 0.55"));

  vector<LocusIterator> good_loci_iterators;
  get_iterators_to_good_loci(loci, good_loci_iterators);

  ProximalLoci proximal_loci(good_loci_iterators, 10);
  
  vector<LocusIterator> neighbors;
  
  proximal_loci.get(neighbors);
  vector<LocusIterator> true_neighbors_a(good_loci_iterators.begin(),
                                         good_loci_iterators.begin() + 3);
  ASSERT_THAT(neighbors, Eq(true_neighbors_a));
  
  proximal_loci.get(neighbors);
  vector<LocusIterator> true_neighbors_b(good_loci_iterators.begin(), 
                                         good_loci_iterators.begin() + 4);
  
  ASSERT_THAT(neighbors, Eq(true_neighbors_b));
  
  for (size_t iter = 0; iter < 4; ++iter)
    proximal_loci.get(neighbors);

  vector<LocusIterator> true_neighbors_c(good_loci_iterators.end() - 2, 
                                         good_loci_iterators.end());
  ASSERT_THAT(neighbors, Eq(true_neighbors_c));
}

TEST(ProximalCpgs, FindsProximalLociOnMultipleChromosome) {
  vector<Locus> loci;
  loci.push_back(Locus("chr2 3 4 c 0.03"));
  loci.push_back(Locus("chr3 5 6 c 0.04"));
  loci.push_back(Locus("chr3 7 8 c 0.05"));
  loci.push_back(Locus("chr3 9 10 c 0.15"));
  loci.push_back(Locus("chr4 11 12 c 0.55"));
  
  vector<LocusIterator> good_loci_iterators;
  get_iterators_to_good_loci(loci, good_loci_iterators);
  
  ProximalLoci proximal_loci(good_loci_iterators, 10);
  vector<LocusIterator> neighbors;
  
  proximal_loci.get(neighbors);
  vector<LocusIterator> true_neighbors_a(good_loci_iterators.begin(), 
                                        good_loci_iterators.begin() + 1);
  ASSERT_THAT(neighbors, Eq(true_neighbors_a));
    
  proximal_loci.get(neighbors);
  vector<LocusIterator> true_neighbors_b(good_loci_iterators.begin() + 1, 
                                         good_loci_iterators.begin() + 4);
  ASSERT_THAT(neighbors, Eq(true_neighbors_b));

  proximal_loci.get(neighbors);
  proximal_loci.get(neighbors);
  ASSERT_THAT(neighbors, Eq(true_neighbors_b));

  proximal_loci.get(neighbors);
  vector<LocusIterator> true_neighbors_c(good_loci_iterators.begin() + 4, 
                                          good_loci_iterators.end());
  ASSERT_THAT(neighbors, Eq(true_neighbors_c));
}
