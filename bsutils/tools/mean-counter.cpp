
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>
#include <cstdlib>


using namespace std;


void
split(const std::string &s,
	  const std::string &delim,
	  vector<string> &result,
	  bool get_empty_fields = false)
{
		result.clear();
		size_t dlen = delim.length();
		size_t slen = s.length();
		size_t start = 0;
		size_t next = 0;
		
		while(start != string::npos)
		{
				next = s.find(delim, start);
				if(next == std::string::npos) next = slen;
				if(start < next || (start == next && get_empty_fields))
						result.push_back(s.substr(start, next - start));
				start = (next + dlen <= slen) ? (next + dlen) : string::npos;
		}
		return;
}


void
split_whitespace(const std::string &s,
				 vector<string> &result)
{
		string local_s = s;
		for(size_t i = 0; i < local_s.length(); i++)
				if(isspace(local_s[i])) local_s[i] = '\t';
		split(local_s, "\t", result);
		return;
}

void
parse_bedline(const string &line,
			  vector<string> &cols)
{
		vector<string> local_strings;
		split_whitespace(line, local_strings);
		swap(cols, local_strings);

		return;
}



int
main(int argn, char ** argv)
{
		string line;
		size_t n(0);
		double total(0);
		vector<string> pre_cols;

		if ( getline(cin, line) )
		{
				parse_bedline(line, pre_cols);
				total = atof(pre_cols.back().c_str());
				n = 1;
		}
		else
		{
				return 0;
		}


		while ( getline(cin, line) )
		{
				vector<string> cols;
				parse_bedline(line, cols);
				
				if (pre_cols[0] == cols[0]
					&& pre_cols[1] == cols[1]
					&& pre_cols[2] == cols[2])
				{
						total += atof(cols.back().c_str());
						++n;
				}
				else
				{
						copy(pre_cols.begin(), pre_cols.end() - 1,
							 ostream_iterator<string>(cout, "\t"));
						cout << "\t" << n
							 << "\t" << total / n << endl;
						swap(pre_cols, cols);
						total = atof(pre_cols.back().c_str());
						n = 1;
				}
		}
		copy(pre_cols.begin(), pre_cols.end() - 1,
			 ostream_iterator<string>(cout, "\t"));
		cout << "\t" << n
			 << "\t" << total / n << endl;
		return 0;
}
