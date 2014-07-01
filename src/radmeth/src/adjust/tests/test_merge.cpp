
#include <vector>
#include <string>
#include <sstream>

#include "gmock/gmock.h"

#include "locus.hpp"
#include "merge.hpp"

using std::vector; using std::string;
using std::istringstream; using std::ostringstream;
using std::stringstream; using std::ostream; 
using std::istream;

using ::testing::Eq;


TEST(Merge, MergeSignificantCpgsOnSameChromosome) {
  string loci_encoding = "chr1  1  2 c:-0.4:-0.1:0.01:0.01  0.01\n"
                         "chr1  3  4 c:-0.2:-0.3:0.01:0.01  0\n"
                         "chr1  5  6 c:10:20:0.05:0.05      0.05\n"
                         "chr1  7  8 c:0.4:0.2:0.01:0.01    0.01\n"
                         "chr1  9 10 c:0.1:0.2:0.01:0.01    0.01\n"
                         "chr1 11 12 c:0.4:0.2:0.01:0.01    0.01";
  istringstream loci_stream(loci_encoding);
  stringstream dmr_stream;
  
  merge(loci_stream, dmr_stream, 0.011);
  
  vector<Locus> dmrs;
  read_loci(dmr_stream, dmrs);

  vector<Locus> expected_dmrs = { Locus("chr1  1   4 dmr:-0.3:-0.2 2"),
                                  Locus("chr1  7  12 dmr:0.3:0.2 3") };
                                  
  ASSERT_THAT(dmrs, Eq(expected_dmrs));
}

TEST(Merge, MergeSignificantCpgsOnMultipleChromosomes) {
  string loci_encoding = "chr1  1  2 c:0.4:0.1:0.01:0.01  0.009\n"
                         "chr3  3  4 c:0.3:0.1:0.01:0.01  0.008\n"
                         "chr3  7  8 c:0.4:0.3:0.01:0.01  0.007\n"
                         "chr3  9 10 c:0.1:0.2:0.01:0.01  0.002\n"
                         "chr3 11 12 c:0.4:0.2:0.01:0.01  0.003";
  istringstream loci_stream(loci_encoding);
  stringstream dmr_stream;
  
  merge(loci_stream, dmr_stream, 0.01);
  
  vector<Locus> dmrs;
  read_loci(dmr_stream, dmrs);

  vector<Locus> expected_dmrs = { Locus("chr1  1   2 dmr:0.4:0.1 1"),
                                  Locus("chr3  3  12 dmr:0.3:0.2 4") };
                                  
  ASSERT_THAT(dmrs, Eq(expected_dmrs));
}

TEST(Merge, MergeSignificantCpgsWithInvalidePvalues) {
  string loci_encoding = "chr1  1  2 c:0.4:0.1:0.01:0.01   0.009\n"
                         "chr1  3  4 c:0.3:0.1:0.01:0.01   0.008\n"
                         "chr1  7  8 c:0.4:0.3:0.01:0.01  -0.007\n"
                         "chr1  9 10 c:0.1:0.2:0.01:0.01   0.002\n"
                         "chr1 11 12 c:0.4:0.2:0.01:0.01   0.003";
  istringstream loci_stream(loci_encoding);
  stringstream dmr_stream;
  
  merge(loci_stream, dmr_stream, 0.01);
  
  vector<Locus> dmrs;
  read_loci(dmr_stream, dmrs);

  vector<Locus> expected_dmrs = { Locus("chr1  1   4 dmr:0.35:0.1 2"),
                                  Locus("chr1  9  12 dmr:0.25:0.2 2") };
                                  
  ASSERT_THAT(dmrs, Eq(expected_dmrs));
}
