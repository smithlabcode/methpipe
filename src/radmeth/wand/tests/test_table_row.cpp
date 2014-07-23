/*    Copyright (C) 2013 University of Southern California and
 *                       Egor Dolzhenko
 *                       Andrew D Smith
 *
 *    Authors: Andrew D. Smith and Egor Dolzhenko
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 */

#include <sstream>
#include <vector>
#include <string>

#include "gmock/gmock.h"

#include "table_row.hpp"

using std::istringstream; using std::vector;
using std::string; using std::istream;

using ::testing::Eq; using ::testing::Pointwise;
using ::testing::Le; using ::testing::Ge;
using ::testing::ElementsAre; using ::testing::Test;

TEST(ReadRow, InitializesTableRowFromString) {
  TableRow row;
  read_row("chr1:10:11 10 5 10 5 10 7 10 7", row);
  ASSERT_THAT(row.chrom, Eq("chr1"));
  ASSERT_THAT(row.begin, Eq(10));
  ASSERT_THAT(row.end, Eq(11));
  ASSERT_THAT(row.total_counts, ElementsAre(10, 10, 10, 10));
  ASSERT_THAT(row.meth_counts, ElementsAre(5, 5, 7, 7));
}

TEST(ReadRow, ThrowsErrorIfChromosomeIsMissing) {
  TableRow row;
  ASSERT_ANY_THROW(read_row(":10:11 10 5 10 5 10 7 10 7", row));
}

TEST(ReadRow, ThrowsErrorIfBeginIsNotNaturalNumber) {
  TableRow row;
  ASSERT_ANY_THROW(read_row("chr1:a:11 10 5 10 5 10 7 10 7", row));
}

TEST(ReadRow, ThrowsErrorIfEndIsNotNaturalNumber) {
  TableRow row;
  ASSERT_ANY_THROW(read_row("chr1:10:a 10 5 10 5 10 7 10 7", row));
}

TEST(ReadRow, ThrowsErrorIfRowNameDoesntHaveTwoColumns) {
  TableRow row;
  ASSERT_ANY_THROW(read_row("chr1:10 10 5 10 5 10 7 10 7", row));
}

TEST(ReadRow, ThrowsErrorIfThereIsTotalThatIsNotNaturalNumber) {
  TableRow row;
  ASSERT_ANY_THROW(read_row("chr1:10:1 10 5 a 5 10 7 10 7", row));
}

TEST(ReadRow, ThrowsErrorIfThereIsMethylationCountThatIsNotNaturalNumber) {
  TableRow row;
  ASSERT_ANY_THROW(read_row("chr1:10:1 10 5 10 a 10 7 10 7", row));
}

TEST(AssertCompatibility, ThrowsIfNumberSamplesAndNumberCountsIsNotSame) {
  istringstream iss("f\n"
                    "s1 1\n"
                    "s2 1");
  Design design(iss);
  TableRow row;
  read_row("chr1:10:11 10 5 10 5 10 7 10 7", row);
  ASSERT_ANY_THROW(assert_compatibility(design, row));
}

TEST(AssertCompatibility, NoThrowIfNumberSamplesAndNumberCountsIsSame) {
  istringstream iss("f\n"
                    "s1 1\n"
                    "s2 1");
  Design design(iss);
  TableRow row;
  read_row("chr1:10:11 10 5 10 5", row);
  ASSERT_NO_THROW(assert_compatibility(design, row));
}

TEST(FlagForLowCoverage, DontFlagIfCoveredInTestFactorSamplesAndOthers) {
  istringstream iss("f test_factor\n"
                    "s1 1 0\n"
                    "s2 1 0\n"
                    "s3 1 0\n"
                    "s4 1 1\n"
                    "s5 1 1\n"
                    "s6 1 1");
  
  Design design(iss);
  TableRow row;
  read_row("chr1:10:11 5 2 0 0 5 2 4 1 0 0 0 0", row);
  
  ASSERT_FALSE(has_low_coverage(design, 1, row));  
}

TEST(FlagForLowCoverage, FlagIfNotCoveredInTestFactorSamples) {
  istringstream iss("f test_factor\n"
                    "s1 1 0\n"
                    "s2 1 0\n"
                    "s3 1 0\n"
                    "s4 1 1\n"
                    "s5 1 1\n"
                    "s6 1 1");
  
  Design design(iss);
  TableRow row;
  read_row("chr1:10:11 5 2 5 1 5 2 0 0 0 0 0 0", row);
  
  ASSERT_TRUE(has_low_coverage(design, 1, row));  
}

TEST(FlagForLowCoverage, FlagIfNotCoveredInOtherSamples) {
  istringstream iss("f test_factor\n"
                    "s1 1 0\n"
                    "s2 1 0\n"
                    "s3 1 0\n"
                    "s4 1 1\n"
                    "s5 1 1\n"
                    "s6 1 1");
  
  Design design(iss);
  TableRow row;
  read_row("chr1:10:11 0 0 0 0 0 0 5 2 5 1 5 2", row);
  
  ASSERT_TRUE(has_low_coverage(design, 1, row));  
}

TEST(FlagForExtremeCounts, FlagIfMaximumMethylationCountsInAllSamples) {
  istringstream iss("f test_factor\n"
                    "s1 1 0\n"
                    "s2 1 0\n"
                    "s3 1 0\n"
                    "s4 1 1\n"
                    "s5 1 1\n"
                    "s6 1 1");
  
  Design design(iss);
  TableRow row;
  read_row("chr1:10:11 5 5 10 10 20 20 15 15 13 13 15 15", row);
  
  ASSERT_TRUE(has_extreme_counts(design, row));  
}

TEST(FlagForExtremeCounts, FlagIfZeroMethylationCountsInAllSamples) {
  istringstream iss("f test_factor\n"
                    "s1 1 0\n"
                    "s2 1 0\n"
                    "s3 1 0\n"
                    "s4 1 1\n"
                    "s5 1 1\n"
                    "s6 1 1");
  
  Design design(iss);
  TableRow row;
  read_row("chr1:10:11 5 0 10 0 20 0 15 0 13 0 15 0", row);
  
  ASSERT_TRUE(has_extreme_counts(design, row));  
}
