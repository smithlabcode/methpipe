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

using std::vector; using std::string;
using std::istringstream; using std::ostringstream;
using std::stringstream; using std::ostream;
using std::istream;

static bool
read_next_significant_cpg(istream &cpg_stream, GenomicRegion &cpg, double cutoff,
                          bool &skipped_any) {

  GenomicRegion region;
  skipped_any = false;
  string cpg_encoding;

  while (getline(cpg_stream, cpg_encoding)) {
    string record, chrom, name, sign;
    size_t position;
    double raw_pval, adjusted_pval, corrected_pval;

    std::istringstream iss(cpg_encoding);
    iss.exceptions(std::ios::failbit);
    iss >> chrom >> position >> sign >> name >> raw_pval
        >> adjusted_pval >> corrected_pval;

    if (0 <= corrected_pval && corrected_pval < cutoff) {
      cpg.set_chrom(chrom);
      cpg.set_start(position);
      cpg.set_end(position + 1);
      return true;
    }
    skipped_any = true;
  }

  return false;
}

void
merge(istream &cpg_stream, ostream &dmr_stream, double cutoff) {

  bool skipped_last_cpg;
  GenomicRegion dmr;
  dmr.set_name("dmr");

  if (!read_next_significant_cpg(cpg_stream, dmr, cutoff, skipped_last_cpg))
    return;

  dmr.set_score(1);

  GenomicRegion cpg;
  cpg.set_name("dmr");
  while(read_next_significant_cpg(cpg_stream, cpg, cutoff, skipped_last_cpg)) {

    if (skipped_last_cpg || cpg.get_chrom() != dmr.get_chrom()) {
        dmr_stream << dmr << std::endl;

      dmr = cpg;
      dmr.set_score(1);

    } else {
      dmr.set_end(cpg.get_end());
      dmr.set_score(dmr.get_score() + 1);
    }
  }

  dmr_stream << dmr << std::endl;
}
