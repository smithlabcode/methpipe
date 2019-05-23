/* symmetric-cpgs: extract the CpG sites from a methcounts output
 * file and produce a new one with the CpGs treated unstranded.
 *
 * Copyright (C) 2014 University of Southern California and
 *                    Andrew D. Smith
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
#include <stdexcept>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "MethpipeSite.hpp"

using std::string;
using std::cout;
using std::cerr;
using std::endl;


inline bool
found_symmetric(const MSite &prev_cpg, const MSite &curr_cpg) {
  // assumes check for CpG already done
  return (prev_cpg.strand == '+' &&
          curr_cpg.strand == '-' &&
          prev_cpg.pos + 1 == curr_cpg.pos);
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
    std::vector<string> leftover_args;
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
      throw std::runtime_error("could not open file: " + filename);

    MSite prev_site, curr_site;
    bool prev_is_good_cpg = false;
    if (in >> prev_site)
      if (prev_site.is_cpg() && (include_mutated || !prev_site.is_mutated()))
        prev_is_good_cpg = true;

    while (in >> curr_site) {
      if (curr_site.is_cpg() && (include_mutated || !curr_site.is_mutated())) {
        if (prev_is_good_cpg && found_symmetric(prev_site, curr_site)) {
          prev_site.add(curr_site);
          out << prev_site << '\n';
        }
        prev_is_good_cpg = true;
      }
      else prev_is_good_cpg = false;
      std::swap(prev_site, curr_site);
    }
  }
  catch (const std::runtime_error &e)  {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
