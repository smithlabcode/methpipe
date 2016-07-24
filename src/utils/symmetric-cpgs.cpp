/*
 *    symmetric-cpgs: extract the CpG sites from a methcounts output
 *    file and produce a new one with the CpGs treated unstranded.
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
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
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cmath>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::round;

struct SiteInfo {
  bool is_cpg() const {return context.compare(0, 3, "CpG") == 0;}
  void add_meth_info(const SiteInfo &other) {
    const size_t other_meth_count =
      static_cast<size_t>(std::round(other.meth*other.total));
    const size_t meth_count =
      static_cast<size_t>(std::round(meth*total));
    total += other.total;
    meth = (total == 0) ? 0.0 :
      static_cast<double>(meth_count + other_meth_count)/total;
  }
  string chrom;
  string context;
  size_t pos;
  size_t total;
  double meth;
  char strand;
};


static std::istream &
operator>>(std::istream &in, SiteInfo &si) {
  string line;
  if (getline(in, line)) {
    std::istringstream iss(line);
    if (!(iss >> si.chrom >> si.pos >> si.strand
          >> si.context >> si.meth >> si.total))
      in.setstate(std::ios_base::badbit);
  }
  return in;
}


static std::ostream &
operator<<(std::ostream &out, const SiteInfo &si) {
  return out << si.chrom << '\t' << si.pos << '\t'
             << si.strand << '\t' << si.context << '\t'
             << si.meth << '\t' << si.total;
}


static bool
not_mutated(const SiteInfo &si) {
  const size_t len = si.context.length();
  //assert(len > 0);
  // when dealing with the first site, the previous one is empty
  if (si.context.empty()) return true;
  return si.context[len-1] != 'x';
}


static bool
found_symmetric(const SiteInfo &first, const SiteInfo &second) {
  return (first.is_cpg() &&
          second.is_cpg() &&
          (first.strand == '+') &&
          (second.strand == '-') &&
          (first.pos + 1 == second.pos));
}


int
main(int argc, const char **argv) {

  try {

    string outfile;
    bool VERBOSE;
    bool include_mutated = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "get CpG sites and make methylation levels symmetric",
                           "<methcounts-file>");
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("muts", 'm', "include mutated CpG sites",
                      false, include_mutated);
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
    const string filename(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    std::ifstream in(filename.c_str());
    if (!in)
      throw SMITHLABException("could not open file: " + filename);

    SiteInfo prev, si;
    while (in >> si) {
      if (found_symmetric(prev, si)) {
        prev.add_meth_info(si);
        if (not_mutated(si) && not_mutated(prev)) {
          out << prev << '\n';
        }
        else if (include_mutated) {
          prev.context = "CpGx";
          out << prev << '\n';
        }
        prev = SiteInfo();
      }
      else {
        if (prev.is_cpg() &&
            (not_mutated(prev) || include_mutated)) {
          if (prev.strand == '-') {
            prev.strand = '+';
            --prev.pos;
          }
          out << prev << '\n';
        }
        prev = si;
      }
    }

    if (prev.is_cpg() &&
        (not_mutated(prev) || include_mutated)) {
      if (prev.strand == '-') {
        prev.strand = '+';
        --prev.pos;
      }
      out << prev << '\n';
    }

  }
  catch (const SMITHLABException &e)  {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
