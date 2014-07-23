
#include <sstream>
#include <string>

#include "gmock/gmock.h"

#include "cpg.hpp"

using std::string;
using ::testing::Eq;

TEST(Cpg, InitializesFromInputStream) {
  string cpg_encoding = "chr1	3000900	+	CpG	0.833333	6";
  Cpg cpg(cpg_encoding);
  ASSERT_THAT(cpg.chrom(), Eq("chr1"));
  ASSERT_THAT(cpg.locus(), Eq(3000900));
  ASSERT_THAT(cpg.total(), Eq(6));
  ASSERT_THAT(cpg.meth(),  Eq(5));
}

TEST(Cpg, ThrowsIfPositionIsNotInteger) {
  string cpg_encoding = "chr1 1.01 + CpG 0.833333  6";
  ASSERT_ANY_THROW(Cpg cpg(cpg_encoding));
}

TEST(Cpg, ThrowsIfInvalidMethylationLevel) {
  string cpg_encoding = "chr1 3000900 + CpG 1.833333  6";
  ASSERT_ANY_THROW(Cpg cpg(cpg_encoding));
}