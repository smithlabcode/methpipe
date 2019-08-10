/* methdiff: Computes probability that individual CpGs have higher
 *           methylation in file A than in file B, where files A and B
 *           are specified on the command line.
 *
 * Copyright (C) 2011-2019 Andrew D Smith
 *
 * Author: Andrew D. Smith
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#include <cmath>
#include <fstream>
#include <utility>
#include <stdexcept>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "zlib_wrapper.hpp"

#include "MethpipeSite.hpp"

#include <gsl/gsl_sf_gamma.h>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::pair;
using std::runtime_error;
using std::min;

using std::ostream_iterator;
using std::ofstream;

static inline double
log_sum_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}


static double
log_hyper_g_greater(size_t meth_a, size_t unmeth_a,
                    size_t meth_b, size_t unmeth_b, size_t k) {
  return  gsl_sf_lnchoose(meth_b + unmeth_b - 1, k) +
    gsl_sf_lnchoose(meth_a + unmeth_a - 1, meth_a + meth_b - 1 - k) -
    gsl_sf_lnchoose(meth_a + unmeth_a + meth_b + unmeth_b - 2,
                    meth_a + meth_b - 1);
}


static double
test_greater_population(const size_t meth_a, const size_t unmeth_a,
                        const size_t meth_b, const size_t unmeth_b) {
  double p = 0;
  for (size_t k = (meth_b > unmeth_a) ? meth_b - unmeth_a : 0; k < meth_b; ++k)
    p = log_sum_log(p, log_hyper_g_greater(meth_a, unmeth_a,
                                           meth_b, unmeth_b, k));
  return exp(p);
}


template <class T> T&
write_methdiff_site(T &out,
                    const MSite &a, const MSite &b, const double diffscore) {

  MSite c(a);
  c.n_reads = a.n_meth();
  c.meth = diffscore;

  std::ostringstream oss;
  oss << c;

  oss << '\t' << a.n_unmeth() // a.n_meth() already output with 'c'
      << '\t' << b.n_meth()
      << '\t' << b.n_unmeth() << endl;

  out << oss.str();
  return out;
}


static bool
same_chrom_and_pos(const MSite &a, const MSite &b) {
  return a.chrom == b.chrom && a.pos == b.pos;
}


template <class T>
static void
process_sites(const bool VERBOSE, igzfstream &in_a, igzfstream &in_b,
              const bool allow_uncovered, const double pseudocount, T &out) {

  MSite a, b;
  string prev_chrom;
  in_b >> b; // load first site for "b"
  while (in_a >> a) {

    if (VERBOSE && a.chrom != prev_chrom)
      cerr << "[processing: " << a.chrom << "]" << endl;

    while (in_b && b < a) in_b >> b; // find appropriate "b" site

    if (same_chrom_and_pos(a, b)) {

      if (allow_uncovered || min(a.n_reads, b.n_reads) > 0) {

        const size_t meth_a = a.n_meth() + pseudocount;
        const size_t unmeth_a = a.n_unmeth() + pseudocount;
        const size_t meth_b = b.n_meth() + pseudocount;
        const size_t unmeth_b = b.n_unmeth() + pseudocount;

        const double diffscore = test_greater_population(meth_b, unmeth_b,
                                                         meth_a, unmeth_a);

        write_methdiff_site(out, a, b, diffscore);
      }
    }
    swap(prev_chrom, a.chrom);
  }
}


int
main(int argc, const char **argv) {

  try {

    string outfile;
    size_t pseudocount = 1;

    // run mode flags
    bool allow_uncovered = true;
    bool VERBOSE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "compute probability that site "
                           "has higher methylation in file A than B",
                           "<methcounts-A> <methcounts-B>");
    opt_parse.add_opt("pseudo", 'p', "pseudocount (default: 1)",
                      false, pseudocount);
    opt_parse.add_opt("nonzero-only", 'A',
                      "process only sites with coveage in both samples",
                      false, allow_uncovered);
    opt_parse.add_opt("out", 'o', "output file", true, outfile);
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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string cpgs_file_a = leftover_args[0];
    const string cpgs_file_b = leftover_args[1];
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "[opening methcounts file: " << cpgs_file_a << "]" << endl;
    igzfstream in_a(cpgs_file_a);
    if (!in_a)
      throw runtime_error("cannot open file: " + cpgs_file_a);

    if (VERBOSE)
      cerr << "[opening methcounts file: " << cpgs_file_b << "]" << endl;
    igzfstream in_b(cpgs_file_b);
    if (!in_b)
      throw runtime_error("cannot open file: " + cpgs_file_b);

    if (outfile.empty() || !has_gz_ext(outfile)) {
      std::ofstream of;
      if (!outfile.empty()) of.open(outfile);
      std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

      process_sites(VERBOSE, in_a, in_b, allow_uncovered, pseudocount, out);
    }
    else {
      ogzfstream out(outfile);
      process_sites(VERBOSE, in_a, in_b, allow_uncovered, pseudocount, out);
    }
  }
  catch (runtime_error &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
