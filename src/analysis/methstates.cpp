/*    methstates: a program for converting read sequences in
 *    MappedRead format files into methylation states at CpGs covered
 *    by those reads
 *
 *    Copyright (C) 2011 University of Southern California and
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
#include <iostream>
#include <algorithm>
#include <numeric>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MappedRead.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::unordered_map;


inline static bool
is_cpg(const string &s, const size_t idx) {
  return toupper(s[idx]) == 'C' && toupper(s[idx + 1]) == 'G';
}


static void
collect_cpgs(const string &s, unordered_map<size_t, size_t> &cpgs) {
  const size_t lim = s.length() - 1;
  size_t cpg_count = 0;
  for (size_t i = 0; i < lim; ++i)
    if (is_cpg(s, i)) {
      cpgs[i] = cpg_count;
      ++cpg_count;
    }
}

static bool
convert_meth_states_pos(const string &chrom,
                        const unordered_map<size_t, size_t> &cpgs,
                        MappedRead &mr,
                        size_t &start_pos, string &seq) {
  const size_t width = mr.r.get_width();
  const size_t offset = mr.r.get_start();

  size_t cpg_count = 0;
  string states;
  size_t first_cpg = std::numeric_limits<size_t>::max();
  //size_t last_cpg = first_cpg;
  for (size_t i = 0; i < width; ++i) {
    if (offset + i < chrom.length() && is_cpg(chrom, offset + i)) {
      if (mr.seq[i] == 'C') {
        states += 'C';
        ++cpg_count;
      }
      else if (mr.seq[i] == 'T') {
        states += 'T';
        ++cpg_count;
      }
      else states += 'N';
      if (first_cpg == std::numeric_limits<size_t>::max()) {
        first_cpg = i;
      }
      //last_cpg = i;
    }
  }
  if (first_cpg != std::numeric_limits<size_t>::max()) {
    start_pos = cpgs.find(offset + first_cpg)->second;
    seq = states;
  }
  return cpg_count > 0;
}


static bool
convert_meth_states_neg(const string &chrom,
                        const unordered_map<size_t, size_t> &cpgs,
                        MappedRead &mr,
                        size_t &start_pos, string &seq) {
  const size_t width = mr.r.get_width();
  const size_t offset = mr.r.get_start();

  // NEED TO TAKE REVERSE COMPLEMENT FOR THE NEGATIVE STRAND ONES!
  revcomp_inplace(mr.seq);

  size_t cpg_count = 0;
  string states;
  size_t first_cpg = std::numeric_limits<size_t>::max();
  //size_t last_cpg = first_cpg;
  for (size_t i = 0; i < width; ++i) {
    if (offset + i > 0 && is_cpg(chrom, offset + i - 1)) {
      if (mr.seq[i] == 'G') {
        states += 'C';
        ++cpg_count;
      }
      else if (mr.seq[i] == 'A') {
        states += 'T';
        ++cpg_count;
      }
      else states += 'N';
      if (first_cpg == std::numeric_limits<size_t>::max()) {
        first_cpg = i;
      }
      //last_cpg = i;
    }
  }
  if (first_cpg != std::numeric_limits<size_t>::max()) {
    const unordered_map<size_t, size_t>::const_iterator
      the_cpg(cpgs.find(offset + first_cpg - 1));
    assert(the_cpg != cpgs.end());
    start_pos = the_cpg->second;
    seq = states;
  }
  return cpg_count > 0;
}



int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;

    string chrom_file;
    string outfile;
    string fasta_suffix = "fa";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], "convert read sequences "
                           "in MappedRead format to methylation states "
                           "at CpGs covered by those reads",
                           "<mapped-reads>");
    opt_parse.add_opt("output", 'o', "output file name", false, outfile);
    opt_parse.add_opt("chrom", 'c', "file or dir of chroms (.fa extn)",
                      true , chrom_file);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    unordered_map<string, string> chrom_files;
    identify_chromosomes(chrom_file, fasta_suffix, chrom_files);
    if (VERBOSE) {
      cerr << "CHROMS:\t" << chrom_files.size() << endl;
    }

    std::ifstream in(mapped_reads_file.c_str());
    if (!in)
      throw SMITHLABException("cannot open input file " + mapped_reads_file);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    unordered_map<size_t, size_t> cpgs;
    vector<string> chrom_names, chroms;
    string chrom_name;
    MappedRead mr;
    GenomicRegion chrom_region("chr0", 0, 0);
    while (!in.eof() && in >> mr) {
      // get the correct chrom if it has changed
      if (chroms.empty() || !mr.r.same_chrom(chrom_region)) {
        const unordered_map<string, string>::const_iterator
          fn(chrom_files.find(mr.r.get_chrom()));
        if (fn == chrom_files.end())
          throw SMITHLABException("could not find chrom: " + mr.r.get_chrom());
        chrom_names.clear();
        chroms.clear();
        read_fasta_file(fn->second.c_str(), chrom_names, chroms);
        if (VERBOSE)
          cerr << "PROCESSING: " << chrom_names.front() << endl;
        collect_cpgs(chroms.front(), cpgs);
        chrom_region.set_chrom(chrom_names.front());
        chrom_name = chrom_names.front();
      }
      size_t start_pos = std::numeric_limits<size_t>::max();
      string seq;
      const bool has_cpgs = mr.r.pos_strand() ?
        convert_meth_states_pos(chroms.front(), cpgs, mr, start_pos, seq) :
        convert_meth_states_neg(chroms.front(), cpgs, mr, start_pos, seq);
      if (has_cpgs)
        out << chrom_name << '\t'
            << start_pos << '\t'
            << seq << '\n';
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

