/* duplicate-remover:
 *
 * Copyright (C) 2013-2018 University of Southern California and
 *                         Andrew D. Smith
 *
 * Authors: Andrew D. Smith, Ben Decato, Song Qiang
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
 */

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "sam_record.hpp"
#include "bsutils.hpp"
#include "zlib_wrapper.hpp"
#include "cigar_utils.hpp"
#include "htslib_wrapper.hpp"

using std::string;
using std::vector;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::unordered_map;
using std::unordered_set;
using std::runtime_error;


inline char
get_strand(const sam_rec &r) {
  return check_flag(r, samflags::read_rc) ? '-' : '+';
}


inline size_t
get_end(const sam_rec &r) {
  return r.pos + cigar_rseq_ops(r.cigar);
}

static bool
precedes_by_start(const sam_rec &a, const sam_rec &b) {
  return (a.rname < b.rname ||
          (a.rname == b.rname && a.pos < b.pos));
}

static bool
precedes_by_end_and_strand(const sam_rec &a, const sam_rec &b) {
  const size_t end_a = get_end(a), end_b = get_end(b);
  return (end_a < end_b ||
          (end_a == end_b && get_strand(a) < get_strand(b)));
}

static bool
equivalent_chrom_and_start(const sam_rec &a, const sam_rec &b) {
  return a.rname == b.rname && a.pos == b.pos;
}

static bool
equivalent_end_and_strand(const sam_rec &a, const sam_rec &b) {
  return get_end(a) == get_end(b) && get_strand(a) == get_strand(b);
}

static void
get_cpgs(const vector<sam_rec> &aln, vector<size_t> &cpg_pos) {
  const size_t lim = aln.front().seq.length();
  for (size_t i = 1; i < lim; ++i) {
    size_t j = 0;
    while (j < aln.size() && !is_cpg(aln[j].seq, i - 1)) ++j;
    if (j < aln.size())
      cpg_pos.push_back(i - 1);
  }
}


static void
get_cytosines(const vector<sam_rec> &aln, vector<size_t> &c_pos) {
  const size_t lim = aln.front().seq.length();
  for (size_t i = 0; i < lim; ++i) {
    size_t j = 0;
    while (j < aln.size() && !is_cytosine(aln[j].seq[i])) ++j;
    if (j < aln.size())
      c_pos.push_back(i);
  }
}


static void
get_meth_patterns(const bool ALL_C,
                  vector<sam_rec> &aln, vector<size_t> &hist) {

  vector<size_t> sites;
  if (ALL_C)
    get_cytosines(aln, sites);
  else get_cpgs(aln, sites);

  unordered_map<string, vector<size_t> > patterns;
  for (size_t i = 0; i < aln.size(); ++i) {
    string s;
    for (size_t j = 0; j < sites.size(); ++j)
      s += (is_cytosine(aln[i].seq[sites[j]]) ? '1' : '0');
    patterns[s].push_back(i);
  }

  unordered_set<size_t> keepers;
  for (auto i(patterns.begin()); i != end(patterns); ++i) {
    const size_t n_dups = i->second.size();
    keepers.insert(i->second[rand() % n_dups]);
    if (hist.size() <= n_dups)
      hist.resize(n_dups + 1);
    hist[n_dups]++;
  }

  size_t j = 0;
  for (size_t i = 0; i < aln.size(); ++i)
    if (keepers.find(i) != end(keepers)) {
      aln[j] = aln[i];
      ++j;
    }
  aln.erase(begin(aln) + j, end(aln));
}


template<class T>
static void
process_inner_buffer(const bool USE_SEQUENCE,
               const bool ALL_C,
               size_t &reads_out,
               size_t &good_bases_out,
               size_t &reads_with_duplicates,
               vector<size_t> &hist,
               vector<sam_rec> &buffer,
               T &out) {
  if (USE_SEQUENCE) {
    const size_t orig_buffer_size = buffer.size();
    get_meth_patterns(ALL_C, buffer, hist);
    for (auto && r : buffer)
      out << r << "\n";
    reads_out += buffer.size();
    good_bases_out += buffer.size()*buffer[0].seq.length();
    reads_with_duplicates += (buffer.size() < orig_buffer_size);
  }
  else {
    const size_t selected = rand() % buffer.size();
    out << buffer[selected] << "\n";
    if (hist.size() <= buffer.size())
      hist.resize(buffer.size() + 1);
    hist[buffer.size()]++;
    good_bases_out += buffer[selected].seq.length();
    ++reads_out;
    reads_with_duplicates += (buffer.size() > 1);
  }
  buffer.clear();
}

template<class T>
static void
process_outer_buffer(const bool USE_SEQUENCE,
                     const bool ALL_C,
                     size_t &reads_out,
                     size_t &good_bases_out,
                     size_t &reads_with_duplicates,
                     vector<size_t> &hist,
                     vector<sam_rec> &outer_buffer,
                     vector<sam_rec> &inner_buffer,
                     T &out) {
  sort(begin(outer_buffer), end(outer_buffer), precedes_by_end_and_strand);
  auto it(begin(outer_buffer));

  // give inner buffer the first element before processing
  inner_buffer.push_back(*it);
  ++it;

  for (; it != end(outer_buffer); ++it) {
    if (!equivalent_end_and_strand(inner_buffer.front(), *it)) {
      process_inner_buffer(USE_SEQUENCE, ALL_C, reads_out, good_bases_out,
                           reads_with_duplicates, hist, inner_buffer, out);
    }
    inner_buffer.push_back(*it);
  }
  process_inner_buffer(USE_SEQUENCE, ALL_C, reads_out, good_bases_out,
                       reads_with_duplicates, hist, inner_buffer, out);
  outer_buffer.clear();
}

template <class T>
static void
duplicate_remover(const bool VERBOSE,
                  const bool USE_SEQUENCE,
                  const bool ALL_C,
                  const bool DISABLE_SORT_TEST,
                  const string &infile,
                  const string &statfile,
                  const string &histfile,
                  T &out) {

  // histogram is tabulated whether or not user requests it
  vector<size_t> hist;

  SAMReader in(infile);

  sam_rec aln;
  if (!(in >> aln))
    throw runtime_error("error reading file: " + infile);

  size_t reads_in = 1;
  size_t reads_out = 0;
  size_t good_bases_in = aln.seq.length();
  size_t good_bases_out = 0;
  size_t reads_with_duplicates = 0;

  // header of input = header of output
  out << in.get_header();

  vector<sam_rec> outer_buffer(1, aln), inner_buffer;
  unordered_set<string> chroms_seen;
  string cur_chrom = "";
  while (in >> aln) {
    ++reads_in;
    good_bases_in += aln.seq.length();
    if (!DISABLE_SORT_TEST) {
      if (precedes_by_start(aln, outer_buffer.front()))
        throw runtime_error("input not properly sorted within chrom:\n" +
                            toa(outer_buffer.front()) + "\n" + toa(aln));

      if (aln.rname != cur_chrom) {
        if (chroms_seen.find(aln.rname) != end(chroms_seen))
          throw runtime_error("input not grouped by chromosomes: " +
                toa(aln));

        cur_chrom = aln.rname;
      }
    }
    if (!equivalent_chrom_and_start(outer_buffer.front(), aln))
      process_outer_buffer(USE_SEQUENCE, ALL_C, reads_out, good_bases_out,
                           reads_with_duplicates, hist, outer_buffer,
                           inner_buffer, out);
    outer_buffer.push_back(aln);
  }

  process_outer_buffer(USE_SEQUENCE, ALL_C, reads_out, good_bases_out,
                       reads_with_duplicates, hist, outer_buffer,
                       inner_buffer, out);

  if (!statfile.empty()) {

    const size_t reads_removed = reads_in - reads_out;
    const double non_dup_fraction =
      static_cast<double>(reads_out - reads_with_duplicates)/reads_in;
    const double duplication_rate =
      static_cast<double>(reads_removed + reads_with_duplicates)/
      reads_with_duplicates;

    ofstream out_stat(statfile);
    out_stat << "total_reads: " << reads_in << endl
             << "total_bases: " << good_bases_in << endl
             << "unique_reads: " << reads_out << endl
             << "unique_read_bases: " << good_bases_out << endl
             << "non_duplicate_fraction: " << non_dup_fraction << endl
             << "duplicate_reads: " << reads_with_duplicates << endl
             << "reads_removed: " << reads_removed << endl
             << "duplication_rate: "
             << duplication_rate << endl;
  }
  if (!histfile.empty()) {
    ofstream out_hist(histfile);
    for (size_t i = 0; i < hist.size(); ++i)
      if (hist[i] > 0)
        out_hist << i << '\t' << hist[i] << '\n';
  }
}

int main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    bool USE_SEQUENCE = false;
    bool ALL_C = false;
    bool DISABLE_SORT_TEST = false;

    size_t the_seed = 408;
    string statfile;
    string histfile;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "program to remove "
                           "duplicate reads from sorted mapped reads",
                           "<in-file> <out-file>", 2);
    opt_parse.add_opt("stats", 'S', "statistics output file", false, statfile);
    opt_parse.add_opt("hist", '\0', "histogram output file for library"
                      " complexity analysis", false, histfile);
    opt_parse.add_opt("seq", 's', "use sequence info", false, USE_SEQUENCE);
    opt_parse.add_opt("all-cytosines", 'A', "use all cytosines (default: CpG)",
                      false, ALL_C);
    opt_parse.add_opt("disable", 'D', "disable sort test",
                      false, DISABLE_SORT_TEST);
    opt_parse.add_opt("seed", 's', "specify random seed", false, the_seed);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.set_show_defaults();
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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }

    const string infile(leftover_args.front());
    const string outfile(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    srand(the_seed);

    std::ofstream out(outfile);
    if (!out)
      throw runtime_error("failed to open output file: " + outfile);

    duplicate_remover(VERBOSE, USE_SEQUENCE, ALL_C, DISABLE_SORT_TEST,
                      infile, statfile, histfile, out);
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
