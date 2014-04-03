
#include <sstream>
#include <vector>
#include <string>

#include "gmock/gmock.h"

#include "cpg.hpp"
#include "merge.hpp"

using std::istringstream; 
using std::vector;
using std::string; 
using std::ostringstream;
using std::istream;
using ::testing::Eq;

TEST(Merge, MergesThreeMethylomeEncodings) {
  istringstream methylome_a("chr1 1 + CpG 0.5 2\n"
                            "chr1 2 + CpG 0.3 4");

  istringstream methylome_b("chr1 1 + CpG 0.5 10\n"
                            "chr1 2 + CpG 0.3 10");

  istringstream methylome_c("chr1 1 + CpG 0.5 8\n"
                            "chr1 2 + CpG 0.1 4");
                            
  vector<istream*> methylomes = { &methylome_a, 
                                  &methylome_b, 
                                  &methylome_c};

  vector<string> names = {"methylome_a", "methylome_b", "methylome_c"};
                                    
  ostringstream count_table;
  
  merge_methylomes(names, methylomes, count_table);
  
  ASSERT_THAT(count_table.str(), Eq("methylome_a\tmethylome_b\tmethylome_c\n"
                                    "chr1:1:2\t2\t1\t10\t5\t8\t4\n"
                                    "chr1:2:3\t4\t1\t10\t3\t4\t0"));
}