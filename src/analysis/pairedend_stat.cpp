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
  if(len == 0)
    return -1; 
  return a_name.compare(0, len-1, b_name, 0, len-1);
}

static void
mask_overlap(string input1, string input2, 
		size_t &reads_count, size_t &mapped_correctly, size_t &incorrect_chr,
                size_t &incorrect_strand, size_t &incorrect_orient, size_t &incorrect_fragment_size,
                vector<size_t> &fragm_distr, const size_t min_fragm, const size_t max_fragm)
{
   ifstream in1(input1.c_str());
   ifstream in2(input2.c_str());
   if (!in1)
      throw RMAPException("cannot open input file " + input1);
   if(!in2)
      throw RMAPException("cannot open input file " + input2);


   MappedRead one, two;
   in1 >> one;
   in2 >> two;
   reads_count += 2;
   while(true){
     //if same read ids
     int x = same_name(one, two);
     if(x == 0){
       size_t fragm_size = (one.r.get_strand() == '+' ? two.r.get_end() - one.r.get_start() : 
                           one.r.get_end() - two.r.get_start());
       if(one.r.get_chrom() != two.r.get_chrom())
         incorrect_chr++;
       else if(one.r.get_strand() == two.r.get_strand())
         incorrect_strand++;
       else if(fragm_size < 0)
         incorrect_orient++;
       else if(fragm_size < min_fragm || fragm_size > max_fragm)
         incorrect_fragment_size++;
       else{
         fragm_distr.push_back(fragm_size);
         mapped_correctly++; 
       }
       in1.peek();
       in2.peek();
       if(in1.eof() || in2.eof())
         break;
       in1 >> one;
       in2 >> two;
       reads_count += 2;
      }//one == two
      else if(x < 0){ 
       in1.peek();
       if(in1.eof()){
         break;
       }  
       in1 >> one;
       reads_count++;
      }//one < two
      else{
       in2.peek();
       if(in2.eof()){
         break;
       }
       in2 >> two;
       reads_count++;
      }//one > two
   }//while 

   while(!in1.eof() && in1 >> one){
     reads_count++;
     in1.peek();
   }
   while(!in2.eof() && in2 >> two){
     reads_count++;
     in2.peek();
   }

   in1.close(); in2.close();
}

static void
stat_fragm_size(vector<size_t> &fragm_distr, double &ave_fs, double &median,
               double &std_fs, size_t &min_fragm, size_t &max_fragm)
{
   double sum = 0.0;
   size_t i = 0;
   min_fragm = 1000;
   max_fragm = 0;
   size_t asize = fragm_distr.size();
   for(i = 0; i < asize; i++){
     sum += fragm_distr[i];
     if(fragm_distr[i] > max_fragm){
       max_fragm = fragm_distr[i];
     }
     if(fragm_distr[i] < min_fragm)
       min_fragm = fragm_distr[i];
   }//for
   ave_fs = sum / (asize + 0.0);
   sort(fragm_distr.begin(), fragm_distr.end());
   median = fragm_distr[asize / 2];
   if(asize % 2 > 0){
     median += fragm_distr[(asize/2) + 1];
     median /= 2.0;
   }

   sum = 0.0;
   for(i = 0; i < asize; i++){
     double dif = fragm_distr[i] - ave_fs;
     sum += dif * dif;
   }
   
   std_fs = (asize - 1 > 0 ? sum / (asize - 1.0) : 0);

}

static void
print(string outfile, size_t &reads_count, size_t &mapped_correctly, 
     size_t &incorrect_chr, size_t &incorrect_strand, size_t &incorrect_orient, 
     size_t &incorrect_fragment_size, double &ave_fragm_size, double &median,
     double &std_fragm_size, size_t &min_fragm, size_t &max_fragm)
{
   std::ofstream of;
   if (!outfile.empty()) of.open(outfile.c_str());
   std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
   out << "Total mapped reads (single count):\t" << reads_count << endl;
   out << "Mapped correctly pairs (pairs count):\t" << mapped_correctly << endl;
   out << "Mapped incorrectly to different chromosomes (pairs count):\t" << incorrect_chr << endl;
   out << "Mapped incorrectly to the same strand (pairs count):\t" << incorrect_strand << endl;
   out << "Mapped incorrectly by orientation (pairs count):\t" << incorrect_orient << endl;
   out << "Mapped incorrectly due to fragment size (pairs count):\t" << incorrect_fragment_size << endl;
   out << endl;
   out << "Fragment size statistics: " << endl;
   out << "Mean:\t" << ave_fragm_size << endl;
   out << "Median:\t" << median << endl;
   out << "STD:\t" << std_fragm_size << endl;
   out << "Min size:\t" << min_fragm << endl;
   out << "Max size:\t" << max_fragm << endl;
   
}//print

int main(int argc, const char** argv){

try {
	string results1;//mapping results file with the first mates
	string results2;//mapping results file with the second mates
	string output;//
	long int set_min = 0;
        long int set_max = 1000;
    bool VERBOSE = false;
  /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], "a program to calculate paired-end mapped statistics.",
			   "");
    opt_parse.add_opt("output", 'o', "output statistics (default: stdout)", 
		      false, output);
    opt_parse.add_opt("mates1", 'u', "mapped results with mapped mates1 ", 
		      false , results1);
    opt_parse.add_opt("mates2", 'c', "mapped results with mapped mates2 ", 
		      false , results2);
    opt_parse.add_opt("min_fragm_size", 'i', "min fragment size (default: 0)",
                      false , set_min);
    opt_parse.add_opt("max_fragm_size", 'a', "max fragment size (default: 1000)",
                      false , set_max);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (argc < 5) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }

    /****************** END COMMAND LINE OPTIONS *****************/

    cout << "Warning: this program must be used before revcompl" << endl;
    size_t mapped_correctly = 0;//pairs count
    size_t reads_count = 0;//single count of total mapped reads
    size_t incorrect_chr = 0;
    size_t incorrect_strand = 0;
    size_t incorrect_orient = 0;
    size_t incorrect_fragment_size = 0;
    vector<size_t> fragm_distr;
    mask_overlap(results1, results2, reads_count, mapped_correctly, incorrect_chr,
                incorrect_strand, incorrect_orient, incorrect_fragment_size, 
                fragm_distr, set_min, set_max);
    double ave_fragm_size = 0, median = 0, std_fragm_size = 0;
    size_t min_fragm = 1000, max_fragm = 0;
    stat_fragm_size(fragm_distr, ave_fragm_size, median, std_fragm_size, min_fragm, max_fragm);
    print(output, reads_count, mapped_correctly, incorrect_chr,
                incorrect_strand, incorrect_orient, incorrect_fragment_size,
		ave_fragm_size, median, std_fragm_size, min_fragm, max_fragm);


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
