
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
output_region(const string &chrom,
			  const size_t start,
			  const size_t end,
			  const string &elements1,
			  const string &elements2)
{
		string elements;
		if (!elements1.empty() && elements1[0] == ',')
				elements += elements1.substr(1);
		else
				elements += elements1;
		
		elements += ";";

		if (!elements1.empty() && elements2[0] == ',')
				elements += elements2.substr(1);
		else
				elements += elements2;
		
		cout << chrom << "\t"
			 << start << "\t"
			 << end << "\t"
			 << elements << endl;
		
		return;
}


int
main(int argn, char ** argv)
{
		ifstream infile_1(argv[1], fstream::in);	
		ifstream infile_2(argv[2], fstream::in);

		string  chrom, chrom1, chrom2,
				elements1, elements2,
				line1, line2;
		size_t  start(0), end(0),
				start1(0), end1(0),
				start2(0), end2(0);	

		getline(infile_1, line1);
		parse_bedline(line1, chrom1, start1, end1);
		getline(infile_2, line2);
		parse_bedline(line2, chrom2, start2, end2);
		
		if (chrom1 < chrom2 || (chrom1 == chrom2 && start1 < start2))
		{
				chrom = chrom1;
				start = start1;
				end = end1;
		}
		else
		{
				chrom = chrom2;
				start = start2;
				end = end2;
		}

		while (!infile_1.eof() && !infile_2.eof())
		{
				if (chrom1 < chrom2 // line1 is in front of line2
					|| (chrom1 == chrom2 && start1 < start2))
				{
						if (chrom1 == chrom
							&& start1 < end) // current region overlapps with the
											 // union region
						{
								elements1 += ("," + line1);
								end = max(end, end1);
						}
						else
						{
								output_region(chrom, start, end, elements1, elements2);
								
								// start a new joined region
								chrom = chrom1;
								start = start1;
								end = end1;
								elements1 = line1;
								elements2 = "";
						}
						
						// read in a line 
						getline(infile_1, line1);
						parse_bedline(line1, chrom1, start1, end1);
				}
				else // line2 is in front of line1 
				{
						if (chrom2 == chrom
							&& start2 < end) // current region overlapps with the
											 // union region
						{
								elements2 += ("," + line2);
								end = max(end, end2);
						}
						else
						{
								output_region(chrom, start, end, elements1, elements2);
								
								// start a new joined region
								chrom = chrom2;
								start = start2;
								end = end2;
								elements1 = "";
								elements2 = line2;
						}
						
						// read in a line 
						getline(infile_2, line2);
						parse_bedline(line2, chrom2, start2, end2);
				}
		}

		output_region(chrom, start, end, elements1, elements2);

		// add all remaining lines from file 1
		while (!infile_1.eof())
		{
				getline(infile_1, line1);
				parse_bedline(line1, chrom1, start1, end1);
				output_region(chrom1, start1, end1, line1, "");
				infile_1.peek();
		}
		
		
		// add all remaining lines from file 2
		while (!infile_2.eof())
		{
				getline(infile_2, line2);
				parse_bedline(line2, chrom2, start2, end2);
				output_region(chrom2, start2, end2, "", line2);
				infile_2.peek();
		}
		
		return 0;
}
