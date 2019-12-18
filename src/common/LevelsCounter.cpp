/* Copyright (C) 2018 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#include "LevelsCounter.hpp"
#include "bsutils.hpp"

#include <string>
#include <sstream>
#include <stdexcept>

using std::string;
using std::to_string;
using std::runtime_error;

void
LevelsCounter::update(const MSite &s) {
  if (s.is_mutated()) {
    ++mutations;
  }
  else if (s.n_reads > 0) {
    ++sites_covered;
    max_depth = std::max(max_depth, s.n_reads);
    total_c += s.n_meth();
    total_t += s.n_reads - s.n_meth();
    mean_agg += s.meth;
    double lower = 0.0, upper = 0.0;
    wilson_ci_for_binomial(alpha, s.n_reads, s.meth, lower, upper);
    called_meth += (lower > 0.5);
    called_unmeth += (upper < 0.5);
  }
  ++total_sites;
}

string
LevelsCounter::tostring() const {
  static const string indent = string(2, ' ');
  const bool good = (sites_covered != 0);
  std::ostringstream oss;
  // directly counted values
  oss << context + ":\n"
      << indent << "total_sites: " << total_sites << '\n'
      << indent << "sites_covered: " << sites_covered << '\n'
      << indent << "total_c: " << total_c << '\n'
      << indent << "total_t: " << total_t << '\n'
      << indent << "max_depth: " << max_depth << '\n'
      << indent << "mutations: " << mutations << '\n'
      << indent << "called_meth: " << called_meth << '\n'
      << indent << "called_unmeth: " << called_unmeth << '\n'
      << indent << "mean_agg: " << mean_agg << '\n';

  // derived values
  oss << indent << "coverage: " << coverage() << '\n'
      << indent << "sites_covered_fraction: "
      << static_cast<double>(sites_covered)/total_sites << '\n'
      << indent << "mean_depth: "
      << static_cast<double>(coverage())/total_sites << '\n'
      << indent << "mean_depth_covered: "
      << static_cast<double>(coverage())/sites_covered << '\n'
      << indent << "mean_meth: "
      << (good ? to_string(mean_meth()) : "NA")  << '\n'
      << indent << "mean_meth_weighted: "
      << (good ? to_string(mean_meth_weighted()) : "NA") << '\n'
      << indent << "fractional_meth: "
      << (good ? to_string(fractional_meth()) : "NA");
  return oss.str();
}

double LevelsCounter::alpha = 0.95;

std::ostream &
operator<<(std::ostream &out, const LevelsCounter &cs) {
  return out << cs.tostring();
}

static void
check_label(const string &observed, const string expected) {
  if (observed != expected)
    throw runtime_error("bad levels format [" + observed + "," + expected + "]");
}

std::istream &
operator>>(std::istream &in, LevelsCounter &cs) {
  in >> cs.context; // get the context
  cs.context = cs.context.substr(0, cs.context.find_first_of(":"));

  string label;
  in >> label >> cs.total_sites; // the total sites
  check_label(label, "total_sites:");

  in >> label >> cs.sites_covered; // the sites covered
  check_label(label, "sites_covered:");

  in >> label >> cs.total_c; // the total c
  check_label(label, "total_c:");

  in >> label >> cs.total_t; // the total t
  check_label(label, "total_t:");

  in >> label >> cs.max_depth; // the max depth
  check_label(label, "max_depth:");

  in >> label >> cs.mutations; // the number of mutations
  check_label(label, "mutations:");

  in >> label >> cs.called_meth; // the number of sites called methylated
  check_label(label, "called_meth:");

  in >> label >> cs.called_unmeth; // the number of sites called unmethylated
  check_label(label, "called_unmeth:");

  in >> label >> cs.mean_agg; // the mean aggregate
  check_label(label, "mean_agg:");

  return in;
}
