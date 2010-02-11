
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;

void
parse_bedline(const string &line,
			  string &chrom,
			  size_t &start,
			  size_t &end)
{
		stringstream ss(line);
		ss >> chrom >> start >> end;
		
		return;
}

void
output_region(const string& chrom,
			  const size_t start,
			  const size_t end,
			  const string &line)
{
		if (start < end)
		{
				cout << chrom << "\t"
					 << start << "\t"
					 << end << "\t"
					 << line  << endl;
		}
		return;
}

int
main(int argn, char ** argv)
{
		ifstream infile_1(argv[1], fstream::in);	
		ifstream infile_2(argv[2], fstream::in);

		string  chrom, chrom1, chrom2,
				line1, line2;
		size_t  start(0), end(0),
				start1(0), end1(0),
				start2(0), end2(0);	

		getline(infile_1, line1);
		parse_bedline(line1, chrom1, start1, end1);
		getline(infile_2, line2);
		parse_bedline(line2, chrom2, start2, end2);

		while (!infile_1.eof() && !infile_2.eof())
		{
				if (chrom1 == chrom2)
				{
						start = max(start1, start2);
						end = min(end1, end2);
						
// 						cout << "intersection " << start << "\t" << end << endl;
// 						cout << "region1 " << start1 << "\t" << end1 << endl;
// 						cout << "region2 " << start2 << "\t" << end2 << endl;
						
						
						if (start < end) // interstion
						{
								// output the first segments
								output_region(chrom1, start1, start, line1);
								
								// update the remaining segments
								start1 = end;
								start2 = end;

// 								if (start == start1)
// 										start1 = end;
// 								else
// 										end1 = start;
// 								start2 = end;
								
//  								cout << "new region 1 " << start1 << "\t" << end1 << endl;
						}
						else if (end1 <= start2) // line1 is in front of line2
						{
								output_region(chrom1, start1, end1, line1);
								getline(infile_1, line1);
								if (!line1.empty())
										parse_bedline(line1, chrom1, start1, end1);
								else
										break;
						}
						else if (end2 <= start1) // line2 is in front of line2
						{
								getline(infile_2, line2);
								if (!line1.empty())
										parse_bedline(line2, chrom2, start2, end2);
								else
										break;
						}
				}
				else if (chrom1 < chrom2)  // line1 is in front of line2
				{
						output_region(chrom1, start1, end1, line1);
						getline(infile_1, line1);
						if (!line1.empty())
								parse_bedline(line1, chrom1, start1, end1);
						else
								break;
				} 
				else if (chrom1 > chrom2)  // line2 is in front of line1
				{
						getline(infile_2, line2);
						if (!line1.empty())
								parse_bedline(line2, chrom2, start2, end2);
						else
								break;
				}
		}

		output_region(chrom1, start1, end1, line1);
		
		// add all remaining lines from file 1
		infile_1.peek();
		while (!infile_1.eof())
		{
				getline(infile_1, line1);
				parse_bedline(line1, chrom1, start1, end1);
				output_region(chrom1, start1, end1, line1);
				infile_1.peek();
		}
		
		return 0;
}
