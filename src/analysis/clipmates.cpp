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

static int 
same_name(const MappedRead &a, const MappedRead &b) {
  const string a_name(a.r.get_name());
  const string b_name(b.r.get_name());
  const size_t len = a_name.length();
  return a_name.compare(0, len-1, b_name, 0, len-1);
}

static bool
mapped_correctly(const MappedRead &a, const MappedRead &b)
{
 return ( a.r.same_chrom(b.r) && ((b.r.get_start() >= a.r.get_start() && b.r.get_start() < a.r.get_end()) || (a.r.get_start() >= b.r.get_start() && a.r.get_start() < b.r.get_end())));
}

static bool
is_N(char c)
{
   return c == 'N' || c == 'n';
}

static char
revcompl(char c){
   if(c == 'A' || c == 'a')
     return 'T';
   else if(c == 'T' || c == 't')
     return 'A';
   else if(c == 'C' || c == 'c')
     return 'G';
   else if(c == 'G' || c == 'g')
     return 'C';
   else 
     return 'N';
}

static void
mask_N(MappedRead &one, MappedRead &two, const size_t overlap_start_one, const size_t overlap_start_two,
       const size_t overlap_size, bool one_plus)
{
   //one is mapped to positive strand, and two to negative
   //one_plus is true if first mate is mapped to positive

   size_t j = two.seq.length() - overlap_start_two - 1 ;
   for(size_t i = overlap_start_one; i < overlap_start_one + overlap_size; i++){
     if(one_plus){
       if(is_N(one.seq[i]))
         one.seq[i] = revcompl(two.seq[j]);
     }//if
     else{
       if(is_N(two.seq[j]))
         two.seq[j] = revcompl(one.seq[i]);
     }//else
     j--;
   }//for i
   if(one_plus)
     fill_n(two.seq.begin() + (two.seq.length() - overlap_start_two - overlap_size), overlap_size, 'N');
   else
     fill_n(one.seq.begin() + overlap_start_one, overlap_size, 'N'); 
}

static void
clip(MappedRead &one, MappedRead &two, size_t &clipped_count){
   if(mapped_correctly(one, two)){
     clipped_count++;
     size_t overlap_start_one = 0, overlap_start_two = 0, overlap_size = 0;
     if(one.r.pos_strand()){
       overlap_start_one = (two.r.get_start() > one.r.get_start()) ? two.r.get_start() - one.r.get_start() : 0;
       overlap_start_two = (two.r.get_start() > one.r.get_start()) ?
                           0 : one.r.get_start() - two.r.get_start();
       overlap_size = one.seq.length() - overlap_start_one;
       if(one.r.get_end() > two.r.get_end()) 
         overlap_size -= one.r.get_end() - two.r.get_end();
       mask_N(one, two, overlap_start_one, overlap_start_two, overlap_size, true);
     }
     else{
       overlap_start_two = (one.r.get_start() > two.r.get_start()) ? one.r.get_start() - two.r.get_start() : 0;
       overlap_start_one = (one.r.get_start() > two.r.get_start()) ? 0 : two.r.get_start() - one.r.get_start();
       overlap_size = two.seq.length() - overlap_start_two;
       if(two.r.get_end() > one.r.get_end())
         overlap_size -= two.r.get_end() - one.r.get_end();
       mask_N(two, one, overlap_start_two, overlap_start_one, overlap_size, false);
     }
      
   }//if mapped correctly and overlap
}

static void
mask_overlap(string input1, string input2, string output1, string output2,
          size_t &reads_count, size_t &same_name_count, size_t &clipped_count)
{
   ifstream in1(input1.c_str());
   ifstream in2(input2.c_str());
   if (!in1)
      throw RMAPException("cannot open input file " + input1);
   if(!in2)
      throw RMAPException("cannot open input file " + input2);

   ofstream out1, out2;
   if(output1.empty() || output2.empty()){
     throw RMAPException("Please provide output files");
   }
   out1.open(output1.c_str(), ios::out);
   out2.open(output2.c_str(), ios::out);

   MappedRead one, two;
   in1 >> one;
   in2 >> two;
   reads_count += 2;
   while(true){
     //if same read ids
     int x = same_name(one, two);
     if(x == 0){
       same_name_count++;
       clip(one, two, clipped_count);
       out1 << one << endl;
       out2 << two << endl;
       in1.peek();
       in2.peek();
       if(in1.eof() || in2.eof())
         break;
       in1 >> one;
       in2 >> two;
       reads_count += 2;
      }//one == two
      else if(x < 0){ 
       out1 << one << endl;
       in1.peek();
       if(in1.eof()){
         out2 << two << endl;
         break;
       }  
       in1 >> one;
       reads_count++;
      }//one < two
      else{
       out2 << two << endl;
       in2.peek();
       if(in2.eof()){
         out1 << one << endl;
         break;
       }
       in2 >> two;
       reads_count++;
      }//one > two
   }//while 

   while(!in1.eof() && in1 >> one){
     out1 << one << endl;
     in1.peek();
   }
   while(!in2.eof() && in2 >> two){
     out2 << two << endl;
     in2.peek();
   }

   in1.close(); in2.close();
   out1.close(); out2.close();
}

int main(int argc, const char** argv){

try {
	string results1;//mapping results file with the first mates
	string results2;//mapping results file with the second mates
	string output;//
	string output2;//mapping results with overlapping ends masked
    bool VERBOSE = false;
  /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], "a program to clip off the ends of mates that overlap. "
			   "Input in BED2 format (BED + <read> + <qual>) that is result of combine progr "
				"Output is in two files BED and FASTQ (mixture of mates, same sizes)",
			   "");
    opt_parse.add_opt("output", 'x', "output mates1 (default: stdout)", 
		      false, output);
    opt_parse.add_opt("out_reads", 'b', "output mates2 (default: stdout)", 
		      true , output2);
    opt_parse.add_opt("mates1", 'u', "mapped results with mapped mates1 ", 
		      false , results1);
    opt_parse.add_opt("mates2", 'c', "mapped results with mapped mates2 ", 
		      false , results2);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (argc < 9) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }

    /****************** END COMMAND LINE OPTIONS *****************/

    cout << "Warning: clipmates must be used before revcompl" << endl;
    size_t clipped_count = 0;
    size_t reads_count = 0;
    size_t same_name_count = 0;
 
    mask_overlap(results1, results2, output, output2, reads_count, same_name_count, clipped_count);

    if (VERBOSE)
      cerr << "READS TOTAL (single count):\t" << reads_count << endl
           << "SAME NAME (in pairs):\t" << same_name_count << endl
           << "CLIPPED (in pairs):\t" << clipped_count << endl;

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
