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

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>

#include "OptionParser.hpp"
#include "smithlab_os.hpp"

#include "merge.hpp"

using std::vector; 
using std::string;
using std::ifstream; 
using std::cout;
using std::istream;
using std::endl;
using std::cerr;

int 
main(int argc, const char **argv) {
  try {
  bool VERBOSE = false;

  /****************** COMMAND LINE OPTIONS ********************/
  OptionParser opt_parse(strip_path(argv[0]), 
                          "Make proportion table from "
                          "methylomes in MethPipe format.",
                          "<methpipe-methlomes>");
  
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
      cerr << "Must specify two or more methylomes to combine." << endl 
           << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
  }
  
  //const string table_filename(leftover_args.front());

  /****************** END COMMAND LINE OPTIONS *****************/
  vector<string> names;
  vector<istream*> methylomes;

  for(vector<string>::const_iterator it = leftover_args.begin(); 
      it != leftover_args.end(); it++) {
    names.push_back(strip_path(*it));
    methylomes.push_back(new ifstream(it->c_str()));
  }  

//  for (size_t ind = 1; ind < argc; ++ind) {
//    names.push_back(argv[ind]);
//    methylomes.push_back(new ifstream(argv[ind]));
//  }
  
  if (!methylomes.empty())
    merge_methylomes(names, methylomes, cout);

  for(size_t ind = 0; ind < methylomes.size(); ++ind)
   delete methylomes[ind];
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
