
#include <vector>

#include "gmock/gmock.h"

#include "locus.hpp"

using std::string; using std::istringstream;
using std::vector;

using testing::ElementsAre;
using testing::Eq;

TEST(ALocus, InitializesFromString) {
  Locus locus("chr1 1 10  c 0.01");
  ASSERT_THAT(locus.chrom(), Eq("chr1"));
  ASSERT_THAT(locus.begin(), Eq(1));
  ASSERT_THAT(locus.end(), Eq(10));
  ASSERT_THAT(locus.name(), Eq("c"));
  ASSERT_THAT(locus.score(), Eq(0.01));
}

TEST(ALocus, ThrowsIfLessThanFiveColumns) {
  ASSERT_ANY_THROW(Locus("chr1 1 10 c"));
}

TEST(ALocus, InitializesFromStringHavingExtraColumns) {
  ASSERT_NO_THROW(Locus("chr1	3000573	3000574	c	0.688883	+"));
}

TEST(ALocus, ThrowsIfBeginSmallerThanEnd) {
  ASSERT_ANY_THROW(Locus("chr1 11 10 c 0.01"));
}

TEST(ReadLoci, LoadsLociFromInputStream) {
  string loc_enc_a = "chr1 1 10  c 0.01\n";
  string loc_enc_b = "chr1 2 3  c 0.1\n";
  string loc_enc_c = "chr1 4 5  c 0.05";
  
  istringstream loci_enc(loc_enc_a + loc_enc_b + loc_enc_c);
  
  vector<Locus> loci;
  
  read_loci(loci_enc, loci);
  
  loc_enc_c = "chr1 4 5  c 0.05";
  
  vector<Locus> expected_loci = {Locus(loc_enc_a), Locus(loc_enc_b), 
                                  Locus(loc_enc_c)};
  
  ASSERT_THAT(loci, Eq(expected_loci));
}

TEST(GetIteratorsToGoodLoci, GetsIteratorsToLociWhoseScoresArePValues) {
  vector<Locus> input_loci = { Locus("chr1 1   2 c 0.03 +"),
                               Locus("chr1 5   6 c 0.004 +"),
                               Locus("chr1 10 11 c 0.01 +"),
                               Locus("chr1 14 15 c 0.03 +"),
                               Locus("chr1 17 18 c -0.1 +"),
                               Locus("chr1 70 71 c 0.85 +"),
                               Locus("chr1 72 73 c 1.05 +"),
                               Locus("chr1 80 81 c 0.55 +") };
  
  vector<LocusIterator> good_loci_iterators;
  
  get_iterators_to_good_loci(input_loci, good_loci_iterators);
  
  ASSERT_THAT(good_loci_iterators.size(), Eq(6));
}
