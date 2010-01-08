/* StringTool: header file
 * Song Qiang <qiang.song@usc.edu> 2009
 *
 */

#ifndef STRING_TOOL_QIANGSONG_HPP
#define STRING_TOOL_QIANGSONG_HPP

#include <cstdlib>

#include <vector>
#include <string>

using std::string;
using std::vector;

namespace StringTool
{
		// split a string into substrings seperated by some delimitor
		void split_whitespace(const string &s, vector<string> &result);
		void split(const string &s,
				   const string &delim,
				   vector<string> &result,
				   bool get_empty_fields = false);

		// strip spaces, newline, return from the begining and end of a string
		void strip(const std::string &s, string &result);


		// type conversion
		template <class T> string toa(T t); //depreciated
		template <class T> string tostring(T t); // recommended
		int atoi(const string &str);
		double atof(const string &str);
}

#endif
