/*
 *    conversion in a bisulfite sequencing experiment
 *
 *    Copyright (C) 2009-2012 University of Southern California and
 *                            Andrew D. Smith
 *
 *    Authors: Song Qiang
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

#include <tr1/cmath>

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MappedRead.hpp"
#include "MethpipeFiles.hpp"


using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::accumulate;
using std::tr1::round;

static bool
read_site(const bool is_new_fmt, std::istream &in, string &chrom,
	  size_t &pos, string &strand, string &seq,
	  double &meth, size_t &coverage) {
  return (is_new_fmt ?
	  methpipe::read_site(in, chrom, pos, strand,
			      seq, meth, coverage) :
	  methpipe::read_site_old(in, chrom, pos, strand,
				  seq, meth, coverage));
}



static void
write_site(const bool is_new_fmt, std::ostream &outf,
	   const string &chrom, const size_t pos, const string &strand,
	   const string &seq, const double meth, const size_t coverage) {
  if (is_new_fmt)
    methpipe::write_site(outf, chrom, pos, strand, seq, meth, coverage);
  else
    methpipe::write_site_old(outf, chrom, pos, strand, seq, meth, coverage);
}



static void
check_consistent_sites(const size_t line_number, const string &file_name,
		       const string &chrom, const size_t pos,
		       const string &strand, const string &context,
		       const string &other_chrom, const size_t other_pos,
		       const string &other_strand, const string &other_context) {
  if (chrom != other_chrom)
    throw SMITHLABException("inconsistent chromosome name "
			    "[line=" + toa(line_number) + ",file="
			    + file_name + "]");
  if (pos != other_pos)
    throw SMITHLABException("inconsistent position "
			    "[line=" + toa(line_number) + ",file="
			    + file_name + "]");

  if (strand != other_strand)
    throw SMITHLABException("inconsistent strand "
			    "[line=" + toa(line_number) + ",file="
			    + file_name + "]");

  if (context != other_context)
    throw SMITHLABException("inconsistent context "
			    "[line=" + toa(line_number) + ",file="
			    + file_name + "]");
}



int
main(int argc, const char **argv) {

  try {

    string outfile;
    string out_stat;
    bool VERBOSE;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "merge multiple methcounts files",
                           "<methcounts-files>");
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("stat", 'S', "file to write statistics",
                      false , out_stat);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
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
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const vector<string> methcounts_files(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    vector<std::ifstream*> infiles(methcounts_files.size());
    for (size_t i = 0; i < methcounts_files.size(); ++i)
      infiles[i] = new std::ifstream(methcounts_files[i].c_str());
    const bool new_methcount_fmt =
      methpipe::is_methpipe_file_single(methcounts_files.front());

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
   
    string chrom, strand, seq;
    size_t pos, coverage;
    double meth;

    size_t line_count = 0;

    while (read_site(new_methcount_fmt, *infiles.front(), chrom, pos,
		     strand, seq, meth, coverage)) {
      ++line_count;

      size_t total_coverage = coverage;
      double total_meth = static_cast<size_t>(round(coverage*meth));

      string other_chrom, other_strand, other_seq;
      size_t other_pos = 0ul;

      for (size_t i = 1; i < infiles.size(); ++i) {
        read_site(new_methcount_fmt, *infiles[i], other_chrom,
		  other_pos, other_strand, other_seq, meth, coverage);

	check_consistent_sites(line_count, methcounts_files[i],
			       chrom, pos, strand, seq,
			       other_chrom, other_pos, other_strand, other_seq);

	total_coverage += coverage;
	total_meth += static_cast<size_t>(round(coverage*meth));
      }

      const double methout = total_meth/std::max(total_coverage, 1ul);
      write_site(new_methcount_fmt, out, chrom, pos, strand,
		 seq, methout, total_coverage);
    }

    for (size_t i = 0; i < infiles.size(); ++i) {
      infiles[i]->close();
      delete infiles[i];
    }
  }
  catch (const SMITHLABException &e)  {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
