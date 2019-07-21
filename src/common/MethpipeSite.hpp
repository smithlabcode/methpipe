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
#include <cmath>

struct MSite {

  MSite() {}
  MSite(const std::string &_chrom,
        const size_t _pos,
        const char _strand,
        const std::string &_context,
        const double _meth,
        const size_t _n_reads) :
    chrom(_chrom), pos(_pos), strand(_strand),
    context(_context), meth(_meth), n_reads(_n_reads) {}
  explicit MSite(const std::string &line);

  std::string chrom;
  size_t pos;
  char strand;
  std::string context;
  double meth;
  size_t n_reads;

  bool operator<(const MSite &other) const {
    int r = chrom.compare(other.chrom);
    return (r < 0 ||
            (r == 0 &&
             (pos < other.pos ||
              (pos == other.pos && strand < other.strand))));
  }

  size_t n_meth() const {return std::round(meth*n_reads);}
  size_t n_unmeth() const {return n_reads - n_meth();}

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
    meth = static_cast<double>(total_c_reads)/std::max(1ul, n_reads);
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
  ////////////////////////////////////////////////////////////////////////
  bool is_cpg() const {
    return context.length() >= 3 &&
      (context[0] == 'C' && context[1] == 'p' && context[2] == 'G');
  }
  bool is_chh() const {
    return context.length() >= 3 &&
      (context[0] == 'C' && context[1] == 'H' && context[2] == 'H');
  }
  bool is_ccg() const {
    return context.length() >= 3 &&
      (context[0] == 'C' && context[1] == 'C' && context[2] == 'G');
  }
  bool is_cxg() const {
    return context.length() >= 3 &&
      (context[0] == 'C' && context[1] == 'X' && context[2] == 'G');
  }
  bool is_mutated() const {
    return context.length() == 4 && context[3] == 'x';
  }

  void set_mutated() {
    if (!is_mutated())
      context += 'x';
  }
  void set_unmutated() {
    if (is_mutated())
      context.resize(context.length() - 1);
  }

  std::string tostring() const;
};

template <class T> T &
operator>>(T &in, MSite &s) {
  std::string line;
  if (getline(in, line))
    s = MSite(line);
  return in;
}

template <class T> T &
operator<<(T &out, const MSite &s) {
  out << s.tostring(); // seems to be an issue returning this directly
  return out;
}

size_t
distance(const MSite &a, const MSite &b);

#endif
