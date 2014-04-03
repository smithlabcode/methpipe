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

// smithlab headers
#include "OptionParser.hpp"
#include "smithlab_os.hpp"

// headers for local code
#include "table_splitter.hpp"

// using block
using std::string;
using std::vector;
using std::istringstream;
using std::cerr;
using std::endl;
using std::ofstream;
using std::stringstream;

static string to_string(size_t num) {
  stringstream ss;
  ss << num;
  return ss.str();
}

int
main(int argc, const char **argv) {

try {
  bool VERBOSE = false;
  size_t rows_per_table = 0;
  
  /****************** COMMAND LINE OPTIONS ********************/
  OptionParser opt_parse(strip_path(argv[0]), "Split a proportion table.",
                          "<proportion-table>");
  
  opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);

  opt_parse.add_opt("num_rows", 'r', "number of rows per table", 
                    true, rows_per_table);

  vector<string> leftover_args;

  std::cerr << "argc = " << argc << std::endl; 

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
  
  const string table_filename(leftover_args.front());

  /****************** END COMMAND LINE OPTIONS *****************/

  // Open the count table.
  std::ifstream table_file(table_filename.c_str());
  if (!table_file)
    throw SMITHLABException("could not open file: " + table_filename);
  

  TableSplitter splitter(table_file, rows_per_table);

  const string base_name = "chunk_table_";

  ofstream table;
  size_t table_num = 1;

  do {
  
    string table_name = base_name + to_string(table_num++) + ".txt";

    if (table.is_open())
      table.close();
    
    table.open (table_name.c_str());
  } while ( splitter.get_table(table) );

  table.close();

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