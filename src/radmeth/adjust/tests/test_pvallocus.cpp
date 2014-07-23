
#include <vector>

#include "gmock/gmock.h"

#include "pvallocus.hpp"

using std::string; using std::istringstream;
using std::vector; using std::ostringstream;

using testing::ElementsAre;
using testing::Eq;

TEST(PvalLocus, DirectInitialization) {
  PvalLocus plocus(0, 3000573, 0.2273);

  ASSERT_THAT(plocus.chrom_ind, Eq(0));
  ASSERT_THAT(plocus.pos, Eq(3000573));
  ASSERT_THAT(plocus.raw_pval, Eq(0.2273));
  ASSERT_THAT(plocus.combined_pval, Eq(0));
  ASSERT_THAT(plocus.corrected_pval, Eq(0));
}

TEST(PvalLoci, GetInitializedFromString) {
  istringstream loci_encoding (
  		"chr1  3000573 3000574 c:0.582447:0.143575 0.2273\n"
        "chr1  3000725 3000726 c:0.246345:0.06131  0.702953\n"
        "chr2  3000900 3000901 c:0.502302:0.11534  0.422782\n" );
  
  vector<PvalLocus> pval_loci;
  initialize_pval_loci(loci_encoding, pval_loci);

  ASSERT_THAT(pval_loci.size(), Eq(3));

  ASSERT_THAT(pval_loci[0].chrom_ind, Eq(0));
  ASSERT_THAT(pval_loci[0].pos, Eq(3000573));
  ASSERT_THAT(pval_loci[0].raw_pval, Eq(0.2273));

  ASSERT_THAT(pval_loci[1].chrom_ind, Eq(0));
  ASSERT_THAT(pval_loci[1].pos, Eq(3000725));
  ASSERT_THAT(pval_loci[1].raw_pval, Eq(0.702953));

  ASSERT_THAT(pval_loci[2].chrom_ind, Eq(1));
  ASSERT_THAT(pval_loci[2].pos, Eq(3000900));
  ASSERT_THAT(pval_loci[2].raw_pval, Eq(0.422782));
}

TEST(PvalLoci, NonPvalueLociSkippedDuringInitialization) {
  istringstream loci_encoding (
  			"chr1  3000573 3000574 c:0.582447:0.143575 0.2273\n"
        "chr1  3000725 3000726 c:0.246345:0.06131  -1\n"
        "chr2  3000900 3000901 c:0.502302:0.11534  0.422782\n" );
  
  vector<PvalLocus> pval_loci;
  initialize_pval_loci(loci_encoding, pval_loci);

  ASSERT_THAT(pval_loci.size(), Eq(2));

  ASSERT_THAT(pval_loci[0].chrom_ind, Eq(0));
  ASSERT_THAT(pval_loci[0].pos, Eq(3000573));
  ASSERT_THAT(pval_loci[0].raw_pval, Eq(0.2273));

  ASSERT_THAT(pval_loci[1].chrom_ind, Eq(1));
  ASSERT_THAT(pval_loci[1].pos, Eq(3000900));
  ASSERT_THAT(pval_loci[1].raw_pval, Eq(0.422782));
}

TEST(PvalLoci, PvalueLociUpdated) {
  
	string input = "chr1  3000573 3000574 c:0.582447:0.143575 0.2273\n"
        		   	 "chr1  3000725 3000726 c:0.246345:0.06131  -1\n"
        		   	 "chr2  3000900 3000901 c:0.502302:0.11534  0.422782\n";

	istringstream loci_encoding(input);
  
	vector<PvalLocus> pval_loci;
	initialize_pval_loci(loci_encoding, pval_loci);
	
	istringstream second_loci_encoding(input);
	ostringstream output_loci_encoding;

	update_pval_loci(second_loci_encoding, pval_loci, output_loci_encoding);

	string output = "chr1\t3000573\t3000574\tc:0.582447:0.143575:0.2273:0\t0\n"
        		   	  "chr1\t3000725\t3000726\tc:0.246345:0.06131\t-1\n"
        		   	  "chr2\t3000900\t3000901\tc:0.502302:0.11534:0.422782:0\t0\n";

  ASSERT_THAT(output_loci_encoding.str(), Eq(output));
}