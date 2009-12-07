/*    cpgsmoother: a program for smoothing CpG methylation frequency
 *    over a genome-scale dataset.
 *
 *    Copyright (C) 2009 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
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

#include <string>
#include <vector>
#include <popt.h>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>

#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"

#include "BSUtils.hpp"

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::istream;
using std::max;


istream& 
operator>>(istream &in, GenomicRegion &region){
  static const size_t buffer_size = 10000; // Magic
  char buffer[buffer_size];
  in.getline(buffer, buffer_size);
  region = GenomicRegion(buffer);
  return in;
}


/*  THIS FUNCTION FILLS A BUFFER FOR GenomicRegion OBJECTS
 */
static void 
fill_buffer(std::ifstream &in, const size_t buffer_start, 
	    vector<GenomicRegion> &buffer) {
  GenomicRegion tmp;
  size_t i = buffer_start;
  assert(buffer_start <= buffer.size());
  for (; i != buffer.size() && !in.eof(); ++i) {
    in >> tmp;
    buffer[i].swap(tmp);
    in.peek();
  }
  if (i < buffer.size())
    buffer.erase(buffer.begin() + i, buffer.end());
}


template <class T>
class FileIterator {
public:
  FileIterator(const string f, const size_t bs);
  void increment_first() {
    if (++first == buffer.end()) {
      assert(first <= last);
      refill_buffer();
    }
  }
  void increment_curr() {
    ++curr;
    if (curr > last) {
      cerr << *curr << endl
	   << *last << endl;
      exit(0);
    }
    assert(first <= curr);
    assert(curr <= last);
  }
  void increment_last() {
    if (++last == buffer.end()) {
      assert(first <= last);
      refill_buffer();
    }
  }
  typename vector<T>::const_iterator get_first() const {return first;}
  typename vector<T>::const_iterator get_curr() const {return curr;}
  typename vector<T>::const_iterator get_last() const {return last;}
  bool first_is_good() const {return (!in.eof() || first < buffer.end());}
  bool curr_is_good() const {return (!in.eof() || curr < buffer.end());}
  bool last_is_good() const {return (!in.eof() || last < buffer.end());}
private:
  std::ifstream in;
  vector<T> buffer;
  typename vector<T>::iterator first;
  typename vector<T>::iterator curr;
  typename vector<T>::iterator last;
  void refill_buffer();
};

/* THIS REFILL BUFFER IS USED WHEN INCREMENTS TO EITHER THE FIRST OR
   THE LAST CURRENTLY USED ELEMENTS IN THE BUFFER HIT THE END OF THE
   BUFFER. HOPEFULLY THE FIRST ONE WILL NOT HIT THE END BEFORE THE
   LAST: THE FIRST IS ALWAYS SUPPOSED TO BE LESS THAN OR EQUAL TO THE
   LAST.
 */
template <class T> void
FileIterator<T>::refill_buffer() {
  assert(first <= curr);
  assert(curr <= last);
  const size_t diff = last - first;
  const size_t diff_curr = curr - first;
  copy(first, last, buffer.begin());
  // Not sure if the code below actualy works or is the best way to
  // grow the buffer
  assert(diff < buffer.size());
  if (diff == buffer.size()) {
    vector<T> newbuff(2*buffer.size());
    copy(buffer.begin(), buffer.end(), newbuff.begin());
    buffer.swap(newbuff);
  }
  first = buffer.begin();
  last = first + diff;
  curr = first + diff_curr;
  fill_buffer(in, diff, buffer);
}


template <class T>
FileIterator<T>::FileIterator(const string f, const size_t bs) :
  buffer(vector<T>(bs)) {
  in.open(f.c_str());
  if (!in) throw RMAPException("cannot open input file " + f);
  fill_buffer(in, 0, buffer);
  first = buffer.begin();
  curr = buffer.begin();
  last = buffer.begin();
}


static bool
precedes(const GenomicRegion &r, const size_t offset) {
  return r.get_end() <= offset;
}


static bool
succeeds(const GenomicRegion &r, const size_t offset) {
  return r.get_start() > offset;
}


static inline double
weight(const double dist, double bandwidth) {
  const double u = dist/bandwidth;
  return 0.75*(1.0 - u*u);
}


static void
scan_cpgs(const size_t context_size, FileIterator<GenomicRegion> &regions,
	  std::ostream &out) {
  // iterate
  while (regions.curr_is_good()) {
    // set the current context
    const size_t curr_pos = regions.get_curr()->get_start();
    const size_t context_lhs = (curr_pos >= context_size) ?
      curr_pos - context_size : 0;
    const size_t context_rhs = curr_pos + context_size;
    
    // adjust the context ends
    while (regions.last_is_good() &&
	   regions.get_curr()->same_chrom(*regions.get_last()) &&
	   !succeeds(*regions.get_last(), context_rhs))
      regions.increment_last();

    while (regions.first_is_good() &&
	   !regions.get_curr()->same_chrom(*regions.get_first()))
      regions.increment_first();
    
    while (regions.first_is_good() &&
	   regions.get_curr()->same_chrom(*regions.get_first()) &&
	   precedes(*regions.get_first(), context_lhs))
      regions.increment_first();
    
    // iterate over the context
    double meth = 0.0, unmeth = 0.0;
    
    for (vector<GenomicRegion>::const_iterator i(regions.get_first());
	 i != regions.get_last(); ++i) {
      const double dist = std::abs(static_cast<int>(i->get_start()) - 
				   static_cast<int>(curr_pos));
      const size_t n_reads =
	std::max(1, atoi(rmap::split(i->get_name(), ":").back().c_str()));
      const double meth_freq = i->get_score();
      // add to the context
      const double w = weight(dist, context_size);
      assert(finite(w));
      meth += w*n_reads*meth_freq;
      unmeth += w*n_reads*(1.0 - meth_freq);
    }
    out << regions.get_curr()->get_chrom() << "\t" 
	<< curr_pos << "\t" << curr_pos + 1 << "\t"
	<< "CpG:" << meth + unmeth << "\t" << meth/(meth + unmeth) << "\t+\n";
    
    regions.increment_curr();
  }
}


int 
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    
    string outfile;
    size_t context_size = 50;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("cpgsmoother", "a program for smoothing CpG "
			   "methylation frequency over a genome-scale "
			   "dataset", "<cpgs-bed>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("context", 'c', "size of CpG context to consider", 
		      false, context_size);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string cpgs_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    FileIterator<GenomicRegion> regions(cpgs_file, 100000);
    
    std::ostream *out = (outfile.empty()) ? &cout : new std::ofstream(outfile.c_str());
    scan_cpgs(context_size, regions, *out);
    if (out != &cout) delete out;
    
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
