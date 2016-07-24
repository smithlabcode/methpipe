/*   duplicate-remover:
 *
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith, Ben Decato, Song Qiang
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

#include <sys/types.h>
#include <unistd.h>

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MappedRead.hpp"
#include "bsutils.hpp"

using std::string;
using std::vector;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::unordered_map;


static bool
precedes(const MappedRead &a, const MappedRead &b) {
  return a.r.get_chrom() < b.r.get_chrom() ||
    (a.r.get_chrom() == b.r.get_chrom() && 
     (a.r.get_start() < b.r.get_start() ||
      (a.r.get_start() == b.r.get_start() && 
       (a.r.get_end() < b.r.get_end() ||
	(a.r.get_end() == b.r.get_end() &&
	 (a.r.get_strand() < b.r.get_strand()))))));
}


static bool
equivalent(const MappedRead &a, const MappedRead &b) {
  return a.r.same_chrom(b.r) &&
    a.r.get_start() == b.r.get_start() &&
    a.r.get_end() == b.r.get_end() &&
    a.r.get_strand() == b.r.get_strand();
}


static void
get_cpgs(const vector<MappedRead> &mr, vector<size_t> &cpg_pos) {
  const size_t lim = mr.front().seq.length();
  for (size_t i = 1; i < lim; ++i) {
    size_t j = 0;
    while (j < mr.size() && !is_cpg(mr[j].seq, i - 1)) ++j;
    if (j < mr.size())
      cpg_pos.push_back(i - 1);
  }
}


static void
get_cytosines(const vector<MappedRead> &mr, vector<size_t> &c_pos) {
  const size_t lim = mr.front().seq.length();
  for (size_t i = 0; i < lim; ++i) {
    size_t j = 0;
    while (j < mr.size() && !is_cytosine(mr[j].seq[i])) ++j;
    if (j < mr.size())
      c_pos.push_back(i);
  }
}


static void
get_meth_patterns(const bool ALL_C, vector<MappedRead> &mr) {
  
  vector<size_t> sites;
  if (ALL_C)
    get_cytosines(mr, sites);
  else get_cpgs(mr, sites);
  
  unordered_map<string, vector<size_t> > patterns;
  for (size_t i = 0; i < mr.size(); ++i) {
    string s;
    for (size_t j = 0; j < sites.size(); ++j)
      s += (is_cytosine(mr[i].seq[sites[j]]) ? '1' : '0');
    patterns[s].push_back(i);
  }
  
  std::unordered_set<size_t> keepers;
  for (unordered_map<string, vector<size_t> >::iterator i(patterns.begin());
       i != patterns.end(); ++i)
    keepers.insert(i->second[rand() % i->second.size()]);
  
  size_t j = 0;
  for (size_t i = 0; i < mr.size(); ++i)
    if (keepers.find(i) != keepers.end()) {
      mr[j] = mr[i];
      ++j;
    }
  mr.erase(mr.begin() + j, mr.end());
}


int main(int argc, const char **argv) {

  try {
    bool VERBOSE = false;
    bool USE_SEQUENCE = false;
    bool ALL_C = false;
    bool DISABLE_SORT_TEST = false;
    bool INPUT_FROM_STDIN = false;
    
    string outfile;
    string statfile;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "program to remove "
			   "duplicate reads from sorted mapped reads", 
			   "<mapped-reads>");
    opt_parse.add_opt("output", 'o', "output file for unique reads",
                      false, outfile);
    opt_parse.add_opt("stdin", '\0', "take input from stdin",
                      false, INPUT_FROM_STDIN);
    opt_parse.add_opt("stats", 'S', "statistics output file", false, statfile);
    opt_parse.add_opt("seq", 's', "use sequence info", false, USE_SEQUENCE);
    opt_parse.add_opt("all-cytosines", 'A', "use all cytosines (default: CpG)", 
		      false, ALL_C);
    opt_parse.add_opt("disable", 'D', "disable sort test", 
		      false, DISABLE_SORT_TEST);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
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
    if (!INPUT_FROM_STDIN && leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    string infile;
    if (!leftover_args.empty())
      infile = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    srand(time(0) + getpid());

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    std::ifstream ifs;
    if (!infile.empty()) ifs.open(infile.c_str());
    std::istream in(infile.empty() ? cin.rdbuf() : ifs.rdbuf());
    
    MappedRead mr;
    if (!(in >> mr))
      throw SMITHLABException("error reading file: " + infile);
    
    size_t reads_in = 0;
    size_t reads_out = 0;
    size_t good_bases_in = 0;
    size_t good_bases_out = 0;
    size_t reads_with_duplicates = 0;
    
    vector<MappedRead> buffer(1, mr);
    while (in >> mr) {
      ++reads_in;
      good_bases_in += mr.seq.length();
      if (!DISABLE_SORT_TEST && precedes(mr, buffer.front()))
	throw SMITHLABException("input not properly sorted:\n" + 
				toa(buffer.front()) + "\n" + toa(mr));
      if (!equivalent(buffer.front(), mr)) {
	if (USE_SEQUENCE) {
	  const size_t orig_buffer_size = buffer.size();
	  get_meth_patterns(ALL_C, buffer); // get the CpGs for the buffer
	  copy(buffer.begin(), buffer.end(), 
	       std::ostream_iterator<MappedRead>(out, "\n"));
	  reads_out += buffer.size();
      good_bases_out += buffer.size();
	  reads_with_duplicates += (buffer.size() < orig_buffer_size);
	}
	else {
	  const size_t selected = rand() % buffer.size();
	  out << buffer[selected] << "\n";
	  good_bases_out += buffer[selected].seq.length();
	  ++reads_out;
	  reads_with_duplicates += (buffer.size() > 1);
	}
	buffer.clear();
      }
      buffer.push_back(mr);
    }
    
    if (USE_SEQUENCE) {
      const size_t orig_buffer_size = buffer.size();
      get_meth_patterns(ALL_C, buffer);
      copy(buffer.begin(), buffer.end(), 
	   std::ostream_iterator<MappedRead>(out, "\n"));
      reads_out += buffer.size();
      reads_with_duplicates += (buffer.size() < orig_buffer_size);
    }
    else {
      const size_t selected = rand() % buffer.size();
      out << buffer[selected] << "\n";
      good_bases_out += buffer[selected].seq.length();
      ++reads_out;
      reads_with_duplicates += (buffer.size() > 1);
    }
    
    if (!statfile.empty()) {
      std::ofstream out_stat(statfile.c_str());    
      out_stat << "TOTAL READS IN:\t" << reads_in + 1 << "\n"
	       << "GOOD BASES IN:\t" << good_bases_in << "\n"
	       << "TOTAL READS OUT:\t" << reads_out << "\n"
	       << "GOOD BASES OUT:\t" << good_bases_out << "\n"
	       << "DUPLICATES REMOVED:\t" << reads_in - reads_out << "\n"
	       << "READS WITH DUPLICATES:\t" << reads_with_duplicates << "\n";
    }
  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
