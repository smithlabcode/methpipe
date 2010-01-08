/* StringTool: implementation
 * Song Qiang <qiang.song@usc.edu> 2009
 *
 */

#include <sstream>

#include <cctype>

#include "StringTool.hpp"

void
StringTool::split_whitespace(const std::string &s,
							 vector<string> &result)
{
		string local_s = s;
		for(size_t i = 0; i < local_s.length(); i++)
				if(isspace(local_s[i])) local_s[i] = '\t';
		StringTool::split(local_s, "\t", result);
		return;
}


void
StringTool::split(const std::string &s,
				  const std::string &delim,
				  vector<string> &result,
				  bool get_empty_fields)
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


// strip the whitespaces from the beginning and the end of the string
void
strip(const std::string& s, string &result) 
{
		const size_t len = s.length();
		size_t i = 0;
		while(i < len && isspace(s[i])) ++i;
		size_t j = len - 1;
		while(j >= i && isspace(s[j])) --j;
		++j;
		if (i == len)
				result =  "";
		else 
				result =  s.substr(i, j - i);
}


// type conversion
template <class T> string		// depreciated. for backward compitable only
StringTool::toa(T t)
{
		std::ostringstream s;
		s << t;
		return s.str();
}

template <class T> string 		// recommened
StringTool::tostring(T t)
{
		std::ostringstream oss;
		oss << t;
		return oss.str();
}

int 
StringTool::atoi(const string &str)
{
		return std::atoi(str.c_str());
}

double 
StringTool::atof(const string &str)
{
		return std::atof(str.c_str());
}
