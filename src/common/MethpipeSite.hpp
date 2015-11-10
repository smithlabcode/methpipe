/*
  Copyright (C) 2015 University of Southern California
  Authors: Andrew D. Smith

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with This program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef METHPIPE_SITE_HPP
#define METHPIPE_SITE_HPP

#include <string>
#include <tr1/cmath>

struct MSite {
  std::string chrom;
  size_t pos;
  char strand;
  std::string context;
  double meth;
  size_t n_reads;

  size_t n_meth() const {return std::tr1::round(meth*n_reads);}

  //////////////////////////////////////////////////////////////
  /// FUNCTIONS BELOW ARE FOR MANIPULATING SYMMETRIC CPG SITES
  //////////////////////////////////////////////////////////////
  void add(const MSite &other) {
    if (!is_mutated() && other.is_mutated())
      context += 'x';
    // ADS: order matters below as n_reads update invalidates n_meth()
    // function until meth has been updated
    const size_t total_c_reads = n_meth() + other.n_meth();
    n_reads += other.n_reads;
    meth = static_cast<double>(total_c_reads)/n_reads;
  }

  // ADS: function below has redundant check for is_cpg, which is
  // expensive and might be ok to remove
  bool is_mate_of(const MSite &first) {
    return (first.pos + 1 == pos && first.is_cpg() && is_cpg() &&
            first.strand == '+' && strand == '-');
  }

  ////////////////////////////////////////////////////////////////////////
  /////  Functions below test the type of site. These are CpG, CHH, and
  /////  CHG divided into two kinds: CCG and CXG, the former including a
  /////  CpG within. Also included is a function that tests if a site
  /////  has a mutation.
  /////  WARNING: None of these functions test for the length of their
  /////  argument std::string, which could cause problems.
  ////////////////////////////////////////////////////////////////////////
  bool is_cpg() const {
    return (context[0] == 'C' && context[1] == 'p' && context[2] == 'G');
  }
  bool is_chh() const {
    return (context[0] == 'C' && context[1] == 'H' && context[2] == 'H');
  }
  bool is_ccg() const {
    return (context[0] == 'C' && context[1] == 'C' && context[2] == 'G');
  }
  bool is_cxg() const {
    return (context[0] == 'C' && context[1] == 'X' && context[2] == 'G');
  }
  bool is_mutated() const {
    return context[3] == 'x';
  }

  std::string tostring() const;

};

std::istream &
operator>>(std::istream &in, MSite &s);

#endif
