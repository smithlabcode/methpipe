/*    Copyright (C) 2013 University of Southern California and
 *                       Egor Dolzhenko
 *                       Andrew D Smith
 *
 *    Authors: Andrew D. Smith and Egor Dolzhenko
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
 */

#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>

#include "GenomicRegion.hpp"

#include "merge.hpp"

using std::vector;
using std::string;
using std::istringstream;
using std::ostringstream;
using std::stringstream;
using std::ostream;
using std::istream;

// Attemps to find the next significant CpG site. Returns true if one was found
// and flase otherwise.
static bool
read_next_significant_cpg(istream &cpg_stream, GenomicRegion &cpg,
                          double cutoff, bool &skipped_any, bool &sig_raw,
                          size_t &test_cov, size_t &test_meth,
                          size_t &rest_cov, size_t &rest_meth) {
  GenomicRegion region;
  skipped_any = false;
  sig_raw = false;
  string cpg_encoding;

  while (getline(cpg_stream, cpg_encoding)) {
    string record, chrom, name, sign;
    size_t position;
    double raw_pval, adjusted_pval, corrected_pval;

    std::istringstream iss(cpg_encoding);
    iss.exceptions(std::ios::failbit);
    iss >> chrom >> position >> sign >> name >> raw_pval
        >> adjusted_pval >> corrected_pval
        >> test_cov >> test_meth >> rest_cov >> rest_meth;

    if (0 <= corrected_pval && corrected_pval < cutoff) {
      cpg.set_chrom(chrom);
      cpg.set_start(position);
      cpg.set_end(position + 1);
      sig_raw = (0 <= raw_pval && raw_pval < cutoff);
      return true;
    }
    skipped_any = true;
  }

  return false;
}

void
merge(istream &cpg_stream, ostream &dmr_stream, double cutoff) {
  bool skipped_last_cpg, sig_raw;
  GenomicRegion dmr;
  dmr.set_name("dmr");
  size_t dmr_test_cov = 0; size_t dmr_test_meth = 0;
  size_t dmr_rest_cov = 0; size_t dmr_rest_meth = 0;
  size_t test_cov = 0; size_t test_meth = 0;
  size_t rest_cov = 0; size_t rest_meth = 0;

  // Find the first significant CpG, or terminate the function if none exist.
  if (!read_next_significant_cpg(cpg_stream, dmr, cutoff, skipped_last_cpg,
                            sig_raw, test_cov, test_meth, rest_cov, rest_meth))
    return;

  dmr.set_score(sig_raw);
  dmr_test_cov += test_cov;
  dmr_test_meth += test_meth;
  dmr_rest_cov += rest_cov;
  dmr_rest_meth += rest_meth;

  GenomicRegion cpg;
  cpg.set_name("dmr");

  while(read_next_significant_cpg(cpg_stream, cpg, cutoff, skipped_last_cpg,
                          sig_raw, test_cov, test_meth, rest_cov, rest_meth)) {

    if (skipped_last_cpg || cpg.get_chrom() != dmr.get_chrom()) {
      if (dmr.get_score() != 0)
        dmr_stream << dmr.get_chrom() << '\t'
                   << dmr.get_start() << '\t'
                   << dmr.get_end()   << '\t'
                   << dmr.get_name()  << '\t'
                   << dmr.get_score() << '\t'
                   << double(dmr_test_meth)/dmr_test_cov -
                      double(dmr_rest_meth)/dmr_rest_cov << std::endl;
      dmr = cpg;
      dmr.set_score(sig_raw);
      dmr_test_cov = test_cov;
      dmr_test_meth = test_meth;
      dmr_rest_cov = rest_cov;
      dmr_rest_meth = rest_meth;
    } else {
      dmr.set_end(cpg.get_end());
      dmr.set_score(dmr.get_score() + sig_raw);
      dmr_test_cov += test_cov;
      dmr_test_meth += test_meth;
      dmr_rest_cov += rest_cov;
      dmr_rest_meth += rest_meth;
    }
  }
  if(dmr.get_score() != 0) {
      std::cerr << dmr.get_score() << std::endl;
      dmr_stream << dmr << std::endl;
  }
}
