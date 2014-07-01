#include <sstream>
#include <vector>
#include <string>

#include "gmock/gmock.h"

#include "table_splitter.hpp"

using std::istringstream; 
using std::vector;
using std::string; 
using std::ostringstream;
using std::istream;
using ::testing::Eq;

TEST(SplitTable, SplitsACountTableIntoThree) {

	istringstream full_table (
		"control_a	control_b	control_c	case_a	case_b	case_c\n"
		"chr1:108:109	 9	6	 10	8	  1	 1	 2  2	  2  1	14 1\n"
		"chr1:114:115	17	7	 10	0	 14	 3	 5  1	  9	 1   7 1\n"
		"chr1:160:161	12	8	 10	5	 17	 4	15 14	 13	 6	 4 4\n"
		"chr1:309:310	 1	1	  1	0	 17	12	12  8	  2	 1	19 8\n"
		"chr1:499:500	 8	4   6	5	 15	 6	14 10	 14	11	15 1\n"
		"chr1:510:511	 0	0	  0	0	 14	 8	 4	0	  5	 3	 5 1\n"
		"chr1:641:642	19	8	 15	6	 13	 5	 4	0	  8	 5	 3 3\n"
		"chr1:646:647	13	8	  1	0	 18	 7	 8	6	 10	 6	 2 1\n"
		"chr1:649:650	14	4	  2	2	 18	 5	 1	1	  0	 0	 5 2\n"
		"chr1:679:689	12	2	  1	1	 17	 4	 1	1	  0	 0	 5 1" );

	TableSplitter splitter(full_table, 3);

	ostringstream table_a, table_b, table_c, table_d;

	string expected_table_a = 
		"control_a	control_b	control_c	case_a	case_b	case_c\n"
		"chr1:108:109	 9	6	 10	8	  1	 1	 2  2	  2  1	14 1\n"
		"chr1:114:115	17	7	 10	0	 14	 3	 5  1	  9	 1   7 1\n"
		"chr1:160:161	12	8	 10	5	 17	 4	15 14	 13	 6	 4 4";

	splitter.get_table(table_a);

	ASSERT_THAT(table_a.str(), Eq(expected_table_a));

	string expected_table_b = 
		"control_a	control_b	control_c	case_a	case_b	case_c\n"
		"chr1:309:310	 1	1	  1	0	 17	12	12  8	  2	 1	19 8\n"
		"chr1:499:500	 8	4   6	5	 15	 6	14 10	 14	11	15 1\n"
		"chr1:510:511	 0	0	  0	0	 14	 8	 4	0	  5	 3	 5 1";

	splitter.get_table(table_b);

	ASSERT_THAT(table_b.str(), Eq(expected_table_b));

	string expected_table_c = 
		"control_a	control_b	control_c	case_a	case_b	case_c\n"
		"chr1:641:642	19	8	 15	6	 13	 5	 4	0	  8	 5	 3 3\n"
		"chr1:646:647	13	8	  1	0	 18	 7	 8	6	 10	 6	 2 1\n"
		"chr1:649:650	14	4	  2	2	 18	 5	 1	1	  0	 0	 5 2";

	splitter.get_table(table_c);

	ASSERT_THAT(table_c.str(), Eq(expected_table_c));

	string expected_table_d = 
		"control_a	control_b	control_c	case_a	case_b	case_c\n"
		"chr1:679:689	12	2	  1	1	 17	 4	 1	1	  0	 0	 5 1";
	splitter.get_table(table_d);

	ASSERT_THAT(table_d.str(), Eq(expected_table_d));
}