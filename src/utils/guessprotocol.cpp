/* guessprotocol: a program for guessing whether a wgbs protocol is
 * original, pbat or random pbat
 *
 * Copyright (C) 2019
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
};

// see if two reads from two ends match to each other (they should
// have the same name)
static bool
mates(const size_t to_ignore_at_end, // in case names have #0/1 name ends
      const FASTQRecord &a, const FASTQRecord &b) {
  assert(to_ignore_at_end < a.name.length());
  return equal(begin(a.name), end(a.name) - to_ignore_at_end, begin(b.name));
}

// Read 4 lines one time from fastq and fill in the FASTQRecord structure
static std::istream&
operator>>(std::istream& s, FASTQRecord &r) {
  if (getline(s, r.name)) {

    if (r.name.empty() || r.name[0] != '@')
      throw std::runtime_error("bad name line: " + r.name);

    r.name = r.name.substr(1, r.name.find_first_of(' '));

    if (!getline(s, r.seq))
      throw runtime_error("failed to read expected seq line");

    string tmp;
    if (!getline(s, tmp))
      throw runtime_error("failed to read expected + line");

    if (!getline(s, tmp))
      throw runtime_error("failed to read expected score line");
  }
  return s;
}


static string
guess_protocol(const double fraction_t_rich) {
  if (fraction_t_rich >= 0.8) {
    return "original";
  }
  if (fraction_t_rich <= 0.2) {
    return "pbat";
  }
  if (fraction_t_rich >= 0.4 && fraction_t_rich <= 0.6) {
    return "random";
  }
  return "inconclusive";
}

int
main_guessprotocol(int argc, const char **argv) {

  try {

    size_t reads_to_check = 1000000;
    size_t name_suffix_len = 0;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "guess whether protocol is ordinary, pbat or random",
                           "<end1-fastq> [<end2-fastq>]");
    opt_parse.add_opt("nreads", 'n', "number of reads in initial check",
                      false, reads_to_check);
    opt_parse.add_opt("ignore", 'i', "length of read name suffix "
                      "to ignore when matching", false, name_suffix_len);
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
    if (opt_parse.about_requested() || leftover_args.size() > 2) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    const vector<string> reads_files(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    if (reads_files.size() == 2) {
      // Input: paired-end reads with end1 and end2
      std::ifstream in1(reads_files.front());
      if (!in1)
        throw runtime_error("cannot open input file: " + reads_files.front());

      std::ifstream in2(reads_files.back());
      if (!in2)
        throw runtime_error("cannot open input file: " + reads_files.back());

      size_t n_pairs = 0;
      size_t t_rich_pairs = 0;

      FASTQRecord end_one, end_two;
      while (in1 >> end_one && in2 >> end_two && n_pairs < reads_to_check) {
        ++n_pairs;

        // two reads should be in paired-ends
        if (!mates(name_suffix_len, end_one, end_two))
          throw runtime_error("expected mates, got: " +
                              end_one.name + " and " + end_two.name);

        const double end_one_a =
          count(begin(end_one.seq), end(end_one.seq), 'A') +
          count(begin(end_one.seq), end(end_one.seq), 'C');
        const double end_one_t =
          count(begin(end_one.seq), end(end_one.seq), 'T') +
          count(begin(end_one.seq), end(end_one.seq), 'G');

        const double end_two_a =
          count(begin(end_two.seq), end(end_two.seq), 'A') +
          count(begin(end_two.seq), end(end_two.seq), 'C');
        const double end_two_t =
          count(begin(end_two.seq), end(end_two.seq), 'T') +
          count(begin(end_two.seq), end(end_two.seq), 'G');

        const double t_rich_count = (end_one_t + end_two_a);
        const double pbat_count = (end_one_a + end_two_t);

        t_rich_pairs += (t_rich_count > pbat_count);
      }
      const double fraction_t_rich = static_cast<double>(t_rich_pairs)/n_pairs;
      cout << guess_protocol(fraction_t_rich) << '\t'
           << "fraction_t_rich=" << fraction_t_rich  << '\t'
           << "t_rich_pairs=" << t_rich_pairs << '\t'
           << "pairs_examined=" << n_pairs << endl;
    }
    else { // if (reads_files.size() == 1)
      // Input: single-end reads
      std::ifstream in(reads_files.front());
      if (!in)
        throw runtime_error("cannot open input file: " + reads_files.front());

      size_t n_reads = 0;
      size_t t_rich_reads = 0;

      FASTQRecord r;
      while (in >> r && n_reads < reads_to_check) {
        ++n_reads;
        const double a = (count(begin(r.seq), end(r.seq), 'A') +
                          count(begin(r.seq), end(r.seq), 'C'));
        const double t = (count(begin(r.seq), end(r.seq), 'T') +
                          count(begin(r.seq), end(r.seq), 'G'));
        t_rich_reads += (t > a);
      }
      const double fraction_t_rich = static_cast<double>(t_rich_reads)/n_reads;
      cout << guess_protocol(fraction_t_rich) << '\t'
           << "fraction_t_rich=" << fraction_t_rich  << '\t'
           << "t_rich_reads=" << t_rich_reads << '\t'
           << "reads_examined=" << n_reads << endl;
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
