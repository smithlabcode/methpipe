/*    levels: a program to compute the three different formulas for
 *    methylation levels described in the paper:
 *
 *        'Leveling' the playing field for analyses of single-base
 *         resolution DNA methylomes
 *         Schultz, Schmitz & Ecker (TIG 2012)
 *
 *    Copyright (C) 2013 University of Southern California and
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
 */

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include <tr1/cmath>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "MethpipeFiles.hpp"

#include "bsutils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;



static void
parse_cpg_line(const string &buffer, size_t &n_meth, size_t &n_unmeth) {
  
  std::istringstream is(buffer);
  string name, dummy;
  double meth_freq = 0.0;
  is >> dummy >> dummy >> dummy >> name >> meth_freq;
  const size_t total = atoi(name.substr(name.find_first_of(":") + 1).c_str());
  n_meth = std::tr1::round(meth_freq*total);
  n_unmeth = std::tr1::round((1.0 - meth_freq)*total);
  assert(n_meth + n_unmeth == total);
}



static bool
get_meth_unmeth(const bool IS_METHPIPE_FILE, std::ifstream &in,
		size_t &n_meth, size_t &n_unmeth) {
  
  if (IS_METHPIPE_FILE) {
    string dummy;
    size_t dummy_pos = 0;
    double meth = 0.0;
    size_t coverage = 0;
    if (!methpipe::read_site(in, dummy, dummy_pos, dummy, 
			     dummy, meth, coverage))
      return false;
    else {
      n_meth = std::tr1::round(meth*coverage);
      n_unmeth = std::tr1::round((1.0 - meth)*coverage);
      assert(n_meth + n_unmeth == coverage);
    }
  }
  else {
    string buffer;
    if (!getline(in, buffer))
      return false;
    parse_cpg_line(buffer, n_meth, n_unmeth);
  }
  return true;
}



int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    bool IS_METHPIPE_FILE = true;
    string outfile;
    double alpha = 0.95;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "compute methylation levels",
			   "<cpgs-file>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("alpha", 'a', "alpha for confidence interval", 
		      false, alpha);
    opt_parse.add_opt("bed", 'b', "file in bed format", 
		      false, IS_METHPIPE_FILE);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (opt_parse.help_requested()) {
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
    const string cpgs_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    std::ifstream in(cpgs_file.c_str());
    if (!in)
      throw SMITHLABException("bad input file: " + cpgs_file);

    size_t mapped_sites = 0, total_sites = 0;
    size_t total_c = 0, total_t = 0;
    size_t called_meth = 0, called_unmeth = 0;
    
    vector<double> meth_levels;
    string buffer;

    size_t n_meth = 0, n_unmeth = 0;
    while (get_meth_unmeth(IS_METHPIPE_FILE,
			   in, n_meth, n_unmeth)) {
      
      if (n_meth + n_unmeth > 0) {
	// get info for weighted mean methylation
	total_c += n_meth;
	total_t += n_unmeth;
	
	// get info for mean methylation
	const size_t N = n_meth + n_unmeth;
	const double level = static_cast<double>(n_meth)/N;
	meth_levels.push_back(level);
	
	// get info for binomial test
	double lower = 0.0, upper = 0.0;
	wilson_ci_for_binomial(alpha, N, level, lower, upper);
	called_meth += (lower > 0.5);
	called_unmeth += (upper < 0.5);
	
	++mapped_sites;
      }
      ++total_sites;
    }
    
    const double weighted_mean_meth = 
      static_cast<double>(total_c)/(total_c + total_t);
    
    const double fractional_meth = 
      static_cast<double>(called_meth)/(called_meth + called_unmeth); 
    
    const double mean_meth = 
      std::accumulate(meth_levels.begin(), meth_levels.end(), 0.0)/
      meth_levels.size();
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    
    out << "mean_meth" << '\t' << mean_meth << endl
	<< "w_mean_meth" << '\t' << weighted_mean_meth << endl
	<< "frac_meth" << '\t' << fractional_meth << '\t'
	<< "(called sites: " << called_meth + called_unmeth << ")" << endl
	<< "frac_mapped" << '\t'
	<< static_cast<double>(mapped_sites)/total_sites << '\t'
	<< "(total sites: " << total_sites << ")" << endl;
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
