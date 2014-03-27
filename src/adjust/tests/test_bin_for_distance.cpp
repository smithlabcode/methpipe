
#include "gmock/gmock.h"

#include "bin_for_distance.hpp"

using ::testing::Eq;

TEST(Binning, CorrectlyCalculatedBin) {
  BinForDistance bin_for_dist("1:16:5");
  ASSERT_THAT(bin_for_dist.which_bin(1), Eq(0));
  ASSERT_THAT(bin_for_dist.which_bin(5), Eq(0));
  ASSERT_THAT(bin_for_dist.which_bin(6), Eq(1));
  ASSERT_THAT(bin_for_dist.which_bin(11), Eq(2));
}

TEST(Binning, DetectsInvalidBin) {
  BinForDistance bin_for_dist("1:16:5");
  ASSERT_THAT(bin_for_dist.which_bin(0), Eq(bin_for_dist.invalid_bin()));
  ASSERT_THAT(bin_for_dist.which_bin(16), bin_for_dist.invalid_bin());
}
