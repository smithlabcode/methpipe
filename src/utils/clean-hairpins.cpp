/* clean-hairpins: a program for identifying and removing hairping
 * reads from paired-end WGBS or RRBS reads.
 *
 * Copyright (C) 2018 Andrew D. Smith, Liz Ji and Jenny Qu
 *
 * Authors: Liz Ji and Andrew D. Smith
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
#include <fstream>
#include <iterator>
#include <algorithm>
#include <stdexcept>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::runtime_error;

// store each read from one end
struct FASTQRecord {
  string name;
  string seq;
  string score;

  string tostring() const {
    std::ostringstream s;
    s << '@'
      << name << '\n'
      << seq << '\n'
      << '+' << '\n'
      << score;
    return s.str();
  }
};

// see if two reads from two ends match to each other (they should
// have the same name)
static bool
mates(const size_t to_ignore_at_end, // in case names have #0/1 name ends
      const FASTQRecord &a, const FASTQRecord &b) {
  assert(to_ignore_at_end < a.name.length());
  return equal(begin(a.name), end(a.name) - to_ignore_at_end, begin(b.name));
}

std::ostream&
operator<<(std::ostream& s, const FASTQRecord &r) {
  return s << r.tostring();
}

// Read 4 lines one time from fastq and fill in the FASTQRecord structure
std::istream&
operator>>(std::istream& s, FASTQRecord &r) {
  if (getline(s, r.name)) {

    if (r.name.empty() || r.name[0] != '@')
      throw std::runtime_error("bad name line: " + r.name);

    r.name = r.name.substr(1);

    if (!getline(s, r.seq))
      throw runtime_error("failed to read expected seq line");

    string tmp;
    if (!getline(s, tmp))
      throw runtime_error("failed to read expected + line");

    if (!getline(s, r.score))
      throw runtime_error("failed to read expected score line");
  }
  return s;
}

////////////////////////////////////////////////////////////////////////////////

static bool
similar_letters_bisulfite_tc_and_ag(const char a, const char b) {
  return (a == b) || (a == 'T' && b == 'C') || (a == 'G' && b == 'A');
}

// compare two reads to detect the overlapped region
static size_t
similarity_both_bisulfite_convsersions(const string &s1, const string &s2) {

  const size_t lim = min(s1.length(), s2.length());

  size_t count = 0;
  for (size_t i = 0; i < lim; ++i)
    if (similar_letters_bisulfite_tc_and_ag(s1[i], s2[i]))
      count++;

  return count;
}


static double
check_hairpins(const size_t name_suffix_len,
               const size_t reads_to_check, const double cutoff,
               const string &reads_file1, const string &reads_file2) {

  // Input: paired-end reads with end1 and end2
  std::ifstream in1(reads_file1);
  if (!in1)
    throw runtime_error("cannot open input file: " + reads_file1);

  std::ifstream in2(reads_file2);
  if (!in2)
    throw runtime_error("cannot open input file: " + reads_file2);

  size_t n_reads = 0;
  size_t n_bad_reads = 0;

  FASTQRecord end_one, end_two;
  while (n_reads < reads_to_check && in1 >> end_one && in2 >> end_two) {
    ++n_reads;
    // two reads should be in paired-ends
    if (!mates(name_suffix_len, end_one, end_two))
      throw runtime_error("expected mates, got:\n" +
                          end_one.tostring() + "\nand:\n" +
                          end_two.tostring());

    // See if inverted duplicates emerge
    const size_t sim =
      similarity_both_bisulfite_convsersions(end_one.seq, end_two.seq);

    const double min_len = min(end_one.seq.length(), end_two.seq.length());
    const double percent_match = sim/min_len;

    if (percent_match > cutoff)
      ++n_bad_reads;
  }

  return static_cast<double>(n_bad_reads)/n_reads;
}



int
main(int argc, const char **argv) {

  try {

    string outfile;
    string stat_outfile;

    double cutoff = 0.95;
    size_t name_suffix_len = 0;

    size_t reads_to_check = 1000000;
    double max_hairpin_rate = 0.1;

    bool VERBOSE = false;
    bool check_first = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "fix and stat invdup/hairping reads",
                           "<end1-fastq> <end2-fastq>");
    opt_parse.add_opt("output", 'o',
                      "output filename", true, outfile);
    opt_parse.add_opt("stat", 's',
                      "stats output filename", true, stat_outfile);
    opt_parse.add_opt("hairpin", 'h',
                      "max hairpin rate", false, max_hairpin_rate);
    opt_parse.add_opt("check", '\0',
                      "check for hairpin contamination", false, check_first);
    opt_parse.add_opt("nreads", 'n',
                      "number of reads in initial check",
                      false, reads_to_check);
    opt_parse.add_opt("cutoff", 'c',
                      "cutoff for calling an invdup (default: 0.95)",
                      false, cutoff);
    opt_parse.add_opt("ignore", 'i', "length of read name suffix "
                      "to ignore when matching", false, name_suffix_len);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested() || leftover_args.size() != 2) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    const string reads_file1(leftover_args.front());
    const string reads_file2(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    double hairpin_fraction = 0.0;
    if (check_first) {
      if (VERBOSE)
        cerr << "testing for hairpin contamination" << endl;

      hairpin_fraction =
        check_hairpins(name_suffix_len,
                       reads_to_check, cutoff, reads_file1, reads_file2);

      if (VERBOSE) {
        if (hairpin_fraction > max_hairpin_rate)
          cerr << "hairpin problem detected" << endl;
        else
          cerr << "no hairpin problem detected" << endl;
      }
    }

    if (!check_first || hairpin_fraction > max_hairpin_rate) {

      // Input: paired-end reads with end1 and end2
      std::ifstream in1(reads_file1);
      if (!in1)
        throw runtime_error("cannot open input file: " + reads_file1);

      std::ifstream in2(reads_file2);
      if (!in2)
        throw runtime_error("cannot open input file: " + reads_file2);

      // output read2 with hairpins removed
      std::ofstream out(outfile);
      if (!out)
        throw runtime_error("cannot open output file: " + outfile);

      size_t n_reads = 0;
      size_t n_bad_reads = 0;
      double sum_percent_match_good = 0.0, sum_percent_match_bad = 0.0;

      FASTQRecord end_one, end_two;
      while (in1 >> end_one && in2 >> end_two) {
        ++n_reads;

        // two reads should be in paired-ends
        if (!mates(name_suffix_len, end_one, end_two))
          throw runtime_error("expected mates, got:\n" +
                              end_one.tostring() + "\nand:\n" +
                              end_two.tostring());

        // See if inverted duplicates emerge
        const size_t sim =
          similarity_both_bisulfite_convsersions(end_one.seq, end_two.seq);

        const double min_len = min(end_one.seq.length(), end_two.seq.length());
        const double percent_match = sim/min_len;

        if (percent_match > cutoff) {

          sum_percent_match_bad += percent_match;
          ++n_bad_reads;

          end_two.seq.clear();
          end_two.score.clear();
        }
        else
          sum_percent_match_good += percent_match;

        out << end_two << '\n';
      }

      std::ofstream stat_out(stat_outfile);
      stat_out << "total_reads_pairs: " << n_reads << endl
               << "hairpin_read_pairs: " << n_bad_reads << endl
               << "hairpin_cutoff: " << cutoff << endl
               << "mean_percent_match_non_hairpin: "
               << sum_percent_match_good/(n_reads - n_bad_reads) << endl
               << "mean_percent_match_hairpin: "
               << sum_percent_match_bad/n_bad_reads << endl;
    }
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
