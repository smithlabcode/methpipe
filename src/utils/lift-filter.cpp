/* lift-filter: process lift results
 *
 * Copyright (C) 2014 University of Southern California and
 *                    Andrew D. Smith
 *
 * Authors: Jenny Qu
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <algorithm>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "MethpipeSite.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::runtime_error;

static bool
same_chrom_pos_strand(const MSite &a, const MSite &b) {
  return a.pos == b.pos && a.chrom == b.chrom && a.strand == b.strand;
}

int
main_lift_filter(int argc, const char **argv) {
  try{
    string pfile;
    bool VERBOSE = false;
    bool UNIQUE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "Process duplicated sites from fast-liftover output",
                           "<methcount file>");
    opt_parse.add_opt("output", 'o', "Output processed methcount", true, pfile);
    opt_parse.add_opt("unique", 'u', "keep unique sites", false, UNIQUE);
    opt_parse.add_opt("verbose", 'v', "print more information", false, VERBOSE);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
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
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string mfile(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ifstream in(mfile);
    if (!in)
      throw runtime_error("cannot open input file: " + mfile);

    std::ofstream out(pfile);
    //if (!of)
    //  throw runtime_error("cannot open output file: " + pfile);
    //std::ostream out(of.rdbuf());

    // read first site
    MSite curr_site;
    if (!(in >> curr_site))
      throw runtime_error("failed reading: " + mfile);

    MSite next_site;
    bool site_is_unique = true;
    while (in >> next_site) {
      if (same_chrom_pos_strand(curr_site, next_site)) {
        site_is_unique = false;
        curr_site.add(next_site);
      }
      else {
        if (!UNIQUE || site_is_unique)
          out << curr_site << endl;
        site_is_unique = true;
        curr_site = next_site;
      }
    }
    if (!UNIQUE || site_is_unique)
      out << curr_site << endl;

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
