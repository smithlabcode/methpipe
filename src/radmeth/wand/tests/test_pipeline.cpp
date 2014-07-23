
#include <sstream>
#include <string>

#include "gmock/gmock.h"

#include "pipeline.hpp"

using std::istringstream; using std::ostringstream;
using std::string;

using ::testing::Eq;

// Note: there is no test for placefolder that the want should put if 
// the p-value is nan (i.e. the regression fails).

TEST(WandPipeline, ChecksThatSamplenamesInProportionTableAndDesignMatch) {
  istringstream table_encoding(                   "s1    s2     s3     s4\n"
                                "chr1:1:2       2   1  10 5  10  7   7  7\n" );
                                
  istringstream design_encoding("f1 f2\n"
                                "s1 1 1\n"
                                "s3 1 1\n"
                                "s2 1 0\n"
                                "s4 1 0\n");
  
  ostringstream out;
  ASSERT_ANY_THROW(wand(design_encoding, table_encoding, "f2", out));
}

TEST(WandPipeline, CalculatesPvalues) {
  istringstream table_encoding(                   "s1    s2     s3     s4\n"
                                "chr1:1:2       2   1  10 5  10  7   7  7\n"
                                "chr1:102:103   7   3   7 4   8  7   9  5\n"
                                "chr1:203:204  10   5   4 4  10  7   8  7\n"
                                "chr1:305:306  10   5  10 5   8  7  10  7\n"
                                "chr2:307:308   7   3   7 4   8  7   9  5\n");
                                
  istringstream design_encoding("f1 f2\n"
                                "s1 1 1\n"
                                "s2 1 1\n"
                                "s3 1 0\n"
                                "s4 1 0\n");
  
  ostringstream out;
  
  wand(design_encoding, table_encoding, "f2", out);
  
  string output = "chr1\t1\t2\tc:-1.54047:0.323533\t0.077052\n"
                  "chr1\t102\t103\tc:-0.87547:0.205883\t0.240856\n"
                  "chr1\t203\t204\tc:-0.664898:0.134904\t0.428651\n"
                  "chr1\t305\t306\tc:-1.25276:0.277778\t0.0726977\n"
                  "chr2\t307\t308\tc:-0.87547:0.205883\t0.240856\n";
  
  ASSERT_THAT(out.str(), Eq(output));
}

TEST(WandPipeline, ThrowsIfTestFactorIsNotPartOfDesign) {
  istringstream table_encoding(                   "s1     s2     s3     s4\n"
                                "chr1:1:2       2   1  10  5  10  7   7  7\n"
                                "chr1:102:103   7   3   7  4   8  7   9  5\n"
                                "chr1:203:204  10   5   4  4  10  7   8  7\n"
                                "chr1:305:306  10   5  10  5   8  7  10  7\n"
                                "chr2:307:308   7   3   7  4   8  7   9  5\n");
                                
  istringstream design_encoding("f1 f2\n"
                                "s1 1 1\n"
                                "s2 1 1\n"
                                "s3 1 0\n"
                                "s4 1 0\n");
  
  ostringstream out;
  
  ASSERT_ANY_THROW(wand(design_encoding, table_encoding, "not a factor", out));
}

TEST(WandPipeline, PutsPlaceholderForLowlyCoveredSamples) {
  istringstream table_encoding(                   "s1     s2     s3     s4\n"
                                "chr1:1:2       2   1  10  5  10  7   7  7\n"
                                "chr1:102:103   0   0   0  0   8  7   9  5\n"
                                "chr1:203:204  10   5   4  4   0  0   0  0\n"
                                "chr1:305:306  10   5  10  5   8  7  10  7\n"
                                "chr2:307:308   0   0   0  0   0  0   0  0\n"
                                "chr2:309:310  10  10   5  5   3  3   1  1\n"
                                "chr2:311:312   5   0  10  0  30  0  10  0\n");
                                
  istringstream design_encoding("f1 f2\n"
                                "s1 1 1\n"
                                "s2 1 1\n"
                                "s3 1 0\n"
                                "s4 1 0\n");
  
  ostringstream out;
  
  wand(design_encoding, table_encoding, "f2", out);
  
  string output = "chr1\t1\t2\tc:-1.54047:0.323533\t0.077052\n"
                  "chr1\t102\t103\tc:0:0\t-1\nchr1\t203\t204\tc:0:0\t-1\n"
                  "chr1\t305\t306\tc:-1.25276:0.277778\t0.0726977\n"
                  "chr2\t307\t308\tc:0:0\t-1\n"
                  "chr2\t309\t310\tc:0:0\t-1\n"
                  "chr2\t311\t312\tc:0:0\t-1\n";
  
  ASSERT_THAT(out.str(), Eq(output));
}
