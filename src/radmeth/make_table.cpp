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

// for the rounding we do... Meng and Qiang both tested it
#include <tr1/cmath>
using std::tr1::round; 

struct cpg_site {
  std::string chrom;
  size_t position;
  size_t total;
  size_t meth;
};


static std::istream &
operator>>(std::istream &in, cpg_site &s) {
  
  // taking whole line in case extra columns added
  string line;
  getline(in, line);
  if (!line.empty()) {
    // now parsing the line
    char sign;
    string name;
    double meth_level;
    std::istringstream enc_stream(line);
    if (!(enc_stream >> s.chrom >> s.position >> sign
          >> name >> meth_level >> s.total))
      throw (SMITHLABException("bad line: \"" + line + "\""));
    
    if (meth_level < 0.0 || meth_level > 1.0)
      throw (SMITHLABException("methylation level outside [0, 1]" + line));
    
    // same rounding method as used in other methpipe programs
    s.meth = static_cast<size_t>(std::tr1::round(meth_level*s.total));
  }
  return in;
}



static bool
get_cpg_info_all_files(const vector<istream*> &methylomes,
                       string &chrom, size_t &position, 
                       vector<double> &totals, vector<double> &meths) {
  
  // these need to be empty for the logic outside this function 
  chrom.clear();
  totals.clear();
  meths.clear();
  
  bool all_streams_good_and_insync = true;
  
  for (vector<istream*>::const_iterator i(methylomes.begin());
       i != methylomes.end(); ++i) {

    // pull the site information from the current file
    cpg_site s;
    *(*i) >> s; // ADS: I hate "*" and double hate this...
    
    if (chrom.empty()) {
      chrom = s.chrom;
      position = s.position;
    }
    
    // make sure the files are still good and the sites consistent 
    all_streams_good_and_insync =
      (all_streams_good_and_insync && (*i)->good() && 
       (s.chrom == chrom && s.position == position));
    
    // add the actual methylation data to the current "row"
    totals.push_back(s.total);
    meths.push_back(s.meth);
  }
  
  return all_streams_good_and_insync;
}



int 
main(int argc, const char **argv) {

  try {
    
    bool VERBOSE = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), 
                           "make proportion table from methcounts format files",
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
      cerr << "must specify two or more methylomes" << endl 
           << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    vector<string> names(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/
    
    if (VERBOSE)
      cerr << "FILES TO PROCESS: " << names.size() << endl;

    // output the header line of the table 
    transform(names.begin(), names.end(), 
              std::ostream_iterator<string>(cout, "\t"),
              std::ptr_fun(&strip_path));
    cout << endl;
    
    // open the files for each methylome
    vector<istream*> methylomes;
    for (size_t i = 0; i < names.size(); ++i)
      methylomes.push_back(new ifstream(names[i].c_str()));
    
    // now build the big table row-by-row
    vector<double> totals, meths;
    string chrom;
    size_t position;
    while (get_cpg_info_all_files(methylomes, chrom, position, totals, meths)) {

      cout << chrom << '\t' << position << '\t' << position + 1;
      
      for (size_t i = 0; i < totals.size(); ++i)
        cout << '\t' << totals[i] << '\t' << meths[i];
      cout << '\n';
      
      // make sure these are ready for next row
      chrom.clear();
      totals.clear();
      meths.clear();
    }
    
    // close all the open files
    for (size_t i = 0; i < methylomes.size(); ++i)
      delete methylomes[i];
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
