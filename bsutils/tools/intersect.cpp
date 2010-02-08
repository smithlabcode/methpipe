
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
			  size_t &end,
			  string &other)
{
		stringstream ss(line);
		ss >> chrom >> start >> end;
		other = ss.str();
		return;
}


int
main(int argn, char ** argv)
{
		ifstream infile_1(argv[1], fstream::in);	
		ifstream infile_2(argv[2], fstream::in);

		string chrom1, chrom2, other1, other2, line;
		size_t  start(0), end(0),
				start1(0), end1(0),
				start2(0), end2(0);	

		getline(infile_1, line);
		parse_bedline(line, chrom1, start1, end1, other1);
		getline(infile_2, line);
		parse_bedline(line, chrom2, start2, end2, other2);

		while (!infile_1.eof() && !infile_2.eof())
		{

				start = max(start1, start2);
				end = min(end1, end2);
				if (chrom1 == chrom2 && start < end)
				{
						cout << chrom1 << "\t"
							 << start << "\t"
							 << end << "\t"
							 << other1 << "\t"
							 << other2 << endl;
				}
				if (( chrom1 < chrom2 ) || (chrom1 == chrom2 && end1 < end2))
				{
						getline(infile_1, line);
						parse_bedline(line, chrom1, start1, end1, other1);
				} 
				else
				{
						getline(infile_2, line);
						parse_bedline(line, chrom2, start2, end2, other2);
				}
		}
}
