/*    climates: a program that reads 2 files with mates in format BED2 (BED + <read> + <qual>)
 *    (this file is result of combine program, sorted by read_id), 
 *	outputs two files: BED-format with start/end changed according to clipped mates and 
 *	FASTQ corresponding reads that can be used in the next steps of the pipeline
 * 
 *    Copyright (C) 2009 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith, Elena Y. Harris
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include<iomanip>
#include<map>
#include<bitset>
#include<ctime>
#include<list>
#include<cmath>
#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "MappedRead.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::ios;


static void
mapped_stat(string input1, long int &reads_count, long int &bases_count, long int &mism_count)
{
   ifstream in1(input1.c_str());
   if (!in1)
      throw RMAPException("cannot open input file " + input1);

   MappedRead one;
   in1.peek();
   while(!in1.eof() && in1 >> one){
     reads_count++;
     bases_count += one.seq.length();
     mism_count += static_cast<long int>(one.r.get_score());
     in1.peek();
   }

   in1.close(); 
}

int main(int argc, const char** argv){

try {
  /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], "a program to count mapping statistics ",
			   "");

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (argc < 2) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    string results1 = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    long int bases_count = 0;
    long int reads_count = 0;
    long int mism_count = 0;
 
    mapped_stat(results1, reads_count, bases_count, mism_count);

      cout << "READS TOTAL:\t" << reads_count << endl
           << "BASES MAPPED:\t" << bases_count << endl
           << "MISMATCHES TOTAL:\t" << mism_count << endl;
    double err_rate = static_cast<double>(mism_count) / static_cast<double>(bases_count);
    cout << "ERROR RATE:\t" << err_rate << endl;
}
  catch (const RMAPException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;

}
