/* methstates: a program for converting read sequences in SAM format
 * files into methylation states at CpGs covered by those reads
 *
 * Copyright (C) 2011-2022 University of Southern California and
 *                         Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <unordered_set>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "htslib_wrapper.hpp"
#include "sam_record.hpp"
#include "cigar_utils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::unordered_map;
using std::unordered_set;
using std::runtime_error;


inline static bool
is_cpg(const string &s, const size_t idx) {
  return toupper(s[idx]) == 'C' && toupper(s[idx + 1]) == 'G';
}


static void
collect_cpgs(const string &s, unordered_map<size_t, size_t> &cpgs) {
  cpgs.clear();
  const size_t lim = s.length() - 1;
  size_t cpg_count = 0;
  for (size_t i = 0; i < lim; ++i)
    if (is_cpg(s, i))
      cpgs[i] = cpg_count++;
}


static bool
convert_meth_states_pos(const string &chrom,
                        const unordered_map<size_t, size_t> &cpgs,
                        const sam_rec &aln,
                        size_t &start_pos, string &states) {
  states.clear();

  const size_t width = cigar_rseq_ops(aln.cigar);
  const size_t offset = aln.pos - 1;

  string the_seq(aln.seq);
  apply_cigar(aln.cigar, the_seq);
  if (the_seq.size() != width)
    throw runtime_error("bad sam record format: " + aln.tostring());

  size_t cpg_count = 0;
  size_t first_cpg = std::numeric_limits<size_t>::max();
  for (size_t i = 0; i < width; ++i) {
    if (offset + i < chrom.length() && is_cpg(chrom, offset + i)) {

      if (the_seq[i] == 'C') {
        states += 'C';
        ++cpg_count;
      }
      else if (the_seq[i] == 'T') {
        states += 'T';
        ++cpg_count;
      }
      else states += 'N';

      if (first_cpg == std::numeric_limits<size_t>::max())
        first_cpg = i;
    }
  }
  if (first_cpg != std::numeric_limits<size_t>::max()) {
    auto the_cpg = cpgs.find(offset + first_cpg);
    if (the_cpg == end(cpgs))
      throw runtime_error("cannot locate site on positive strand: " +
                          aln.tostring());
    start_pos = the_cpg->second;
  }
  return cpg_count > 0;
}


static bool
convert_meth_states_neg(const string &chrom,
                        const unordered_map<size_t, size_t> &cpgs,
                        const sam_rec &aln,
                        size_t &start_pos, string &states) {
  states.clear();

  const size_t width = cigar_rseq_ops(aln.cigar);
  const size_t offset = aln.pos - 1;

  string the_seq(aln.seq);
  revcomp_inplace(the_seq);
  apply_cigar(aln.cigar, the_seq);
  if (the_seq.size() != width)
    throw runtime_error("bad sam record format: " + aln.tostring());

  size_t cpg_count = 0;
  size_t first_cpg = std::numeric_limits<size_t>::max();
  for (size_t i = 0; i < width; ++i) {
    if (offset + i > 0 && is_cpg(chrom, offset + i - 1)) {
      if (the_seq[i] == 'G') {
        states += 'C';
        ++cpg_count;
      }
      else if (the_seq[i] == 'A') {
        states += 'T';
        ++cpg_count;
      }
      else states += 'N';
      if (first_cpg == std::numeric_limits<size_t>::max()) {
        first_cpg = i;
      }
    }
  }
  if (first_cpg != std::numeric_limits<size_t>::max()) {
    auto the_cpg = cpgs.find(offset + first_cpg - 1);
    if (the_cpg == end(cpgs))
      throw runtime_error("cannot locate site on negative strand: " +
                          aln.tostring());
    start_pos = the_cpg->second;
  }
  return cpg_count > 0;
}


static void
get_chrom(const string &chrom_name,
          const vector<string> &all_chroms,
          const unordered_map<string, size_t> &chrom_lookup,
          string &chrom) {

  auto the_chrom = chrom_lookup.find(chrom_name);
  if (the_chrom == end(chrom_lookup))
    throw runtime_error("could not find chrom: " + chrom_name);

  chrom = all_chroms[the_chrom->second];
  if (chrom.empty())
    throw runtime_error("problem with chrom: " + chrom_name);
}


int
main_methstates(int argc, const char **argv) {

  try {

    const string description =
      "Convert mapped reads in SAM format into a format that indicates binary \
      sequences of methylation states in each read, indexed by the identity   \
      of the CpG they cover, along with the chromosome. Only reads that       \
      cover a CpG site are included in the output. All output is relative to  \
      the positive reference strand. This format is used as input to other    \
      tools, and is not intended to be human-interpretable. All chromosome    \
      sequences are loaded at once.";

    bool VERBOSE = false;

    string chrom_file;
    string outfile;
    string fasta_suffix = "fa";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], description, "<sam-file>");
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

    vector<string> chrom_files;
    if (isdir(chrom_file.c_str()))
      read_dir(chrom_file, fasta_suffix, chrom_files);
    else chrom_files.push_back(chrom_file);

    if (VERBOSE)
      cerr << "n_chrom_files:\t" << chrom_files.size() << endl;

    /* first load in all the chromosome sequences and names, and make
       a map from chromosome name to the location of the chromosome
       itself */
    vector<string> all_chroms;
    vector<string> chrom_names;
    unordered_map<string, size_t> chrom_lookup;
    for (auto i(begin(chrom_files)); i != end(chrom_files); ++i) {
      vector<string> tmp_chroms, tmp_names;
      read_fasta_file_short_names(*i, tmp_names, tmp_chroms);
      for (size_t j = 0; j < tmp_chroms.size(); ++j) {
        if (chrom_lookup.find(tmp_names[j]) != end(chrom_lookup))
          throw runtime_error("repeated chromosome or name: " + tmp_names[j]);
        chrom_names.push_back(tmp_names[j]);
        chrom_lookup[chrom_names.back()] = all_chroms.size();
        all_chroms.push_back("");
        all_chroms.back().swap(tmp_chroms[j]);
      }
    }

    if (VERBOSE)
      cerr << "n_chroms: " << all_chroms.size() << endl;

    SAMReader in(mapped_reads_file);
    if (!in)
      throw runtime_error("cannot open input file " + mapped_reads_file);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    // given a chrom, get cpg location from cpg index
    unordered_map<size_t, size_t> cpgs;

    unordered_set<string> chroms_seen;
    string chrom_name;
    string chrom;

    // iterate over records/reads in the SAM file, sequentially
    // processing each before considering the next
    sam_rec aln;
    while (in >> aln) {

      // get the correct chrom if it has changed
      if (aln.rname != chrom_name) {
        // make sure all reads from same chrom are contiguous in the file
        if (chroms_seen.find(aln.rname) != end(chroms_seen))
          throw runtime_error("chroms out of order (check SAM file sorted)");

        chrom_name = aln.rname;
        if (VERBOSE)
          cerr << "processing " << chrom_name << endl;

        get_chrom(chrom_name, all_chroms, chrom_lookup, chrom);
        collect_cpgs(chrom, cpgs);
      }

      size_t start_pos = std::numeric_limits<size_t>::max();
      string seq;

      const bool has_cpgs = check_flag(aln, samflags::read_rc) ?
        convert_meth_states_neg(chrom, cpgs, aln, start_pos, seq) :
        convert_meth_states_pos(chrom, cpgs, aln, start_pos, seq);

      if (has_cpgs)
        out << aln.rname << '\t'
            << start_pos << '\t'
            << seq << '\n';
    }
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
