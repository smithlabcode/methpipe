/*
 *    Copyright (C) 2013 University of Southern California and
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <tr1/cmath>

#include "OptionParser.hpp"
#include "smithlab_os.hpp"

using std::vector;
using std::string;
using std::ifstream;
using std::istream;
using std::ostream;
using std::endl;
using std::cout;
using std::cerr;
using std::tr1::round;

struct cpg_site {
  std::string chrom;
  size_t position;
  size_t total;
  size_t meth;
};

static std::istream &
operator>>(std::istream &in, cpg_site &s) {
  // Taking entire line in case extra columns added
  string line;
  getline(in, line);
  if (!line.empty()) {
    // Now parsing the line.
    char sign;
    string name;
    double meth_level;
    std::istringstream enc_stream(line);
    if (!(enc_stream >> s.chrom >> s.position >> sign
          >> name >> meth_level >> s.total))
      throw (SMITHLABException("bad line: \"" + line + "\""));

    if (meth_level < 0.0 || meth_level > 1.0)
      throw (SMITHLABException("methylation level outside [0, 1]: " + line));

    // Same rounding method as used in other methpipe programs.
    s.meth = static_cast<size_t>(std::tr1::round(meth_level*s.total));
  }
  return in;
}

int
main(int argc, const char **argv) {

  try {
    bool VERBOSE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "make proportion table from methcounts format files",
                           "<methcount-files>");
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
    if (leftover_args.size() <= 1) {
      cerr << "must specify two or more methylomes" << endl
           << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    vector<string> names(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "[COMBINING " << names.size() << " METHYLOMES]." << endl;

    // Output the header line of the table.
    transform(names.begin(), names.end(),
              std::ostream_iterator<string>(cout, "\t"),
              std::ptr_fun(&strip_path));
    cout << endl;

    // Open the files for each methylome.
    vector<ifstream*> methylomes(names.size());

    for (size_t i = 0; i < names.size(); ++i)
      methylomes[i] = new ifstream(names[i].c_str());

    bool all_streams_good = true;
    while (all_streams_good) {
      cpg_site first_s;
      for (size_t i = 0; i < methylomes.size(); ++i) {
        cpg_site s;
        // Pull the site information from the current file.
        if (*(methylomes[i]) >> s && !s.chrom.empty()) {

          if (i == 0) {
            first_s = s;
            cout << s.chrom << ':' << s.position;
          }
          else
            assert(first_s.chrom == s.chrom && first_s.position == s.position);

          cout << '\t' << s.total << '\t' << s.meth
               << (i == methylomes.size() - 1 ? "\n" : "");
        } else
          all_streams_good = false;
      }
    }

    // Close all the open files.
    for (size_t i = 0; i < methylomes.size(); ++i) {
      methylomes[i]->close();
      delete methylomes[i];
    }
  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
