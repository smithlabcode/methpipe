/*    cgicaller: a program for calling CGI methylation status in a defined
 *    region from bisulfite capture sequencing with Solexa reads
 *
 *    Copyright (C) 2009 University of Southern California and
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
#include <popt.h>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <iomanip>

#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"
#include "RNG.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::pair;
using std::make_pair;


void
separate_cpgs(const vector<GenomicRegion> &cpgs, 
	      const vector<GenomicRegion> &regions,
	      vector<vector<GenomicRegion> > &sep_cpgs) {
  
  unordered_map<string, size_t> region_id_lookup;
  for (size_t i = 0; i < regions.size(); ++i)
    region_id_lookup[regions[i].get_name()] = i;
  
  sep_cpgs.resize(regions.size());
  for (size_t i = 0; i < cpgs.size(); ++i) {
    const string name(cpgs[i].get_name());
    const size_t region_name_end = name.find_first_of(":");
    const string region_name(name.substr(0, region_name_end));
    const unordered_map<string, size_t>::const_iterator j =
      region_id_lookup.find(region_name);
    if (j == region_id_lookup.end())
      throw RMAPException("ERROR: region " + 
			  region_name + " not found (did "
			  "input come from CpG caller?)");
    sep_cpgs[j->second].push_back(cpgs[i]);
  }
}

struct CPGInfo {
  CPGInfo(const GenomicRegion &r) {
    const string name(r.get_name());
    vector<string> parts = rmap::split(name, ":");
    
    assert(parts.size() >= 3);
    unmeth_count = atoi(parts[1].c_str());
    meth_count = atoi(parts[2].c_str());
  }
  bool has_values() const {return meth_count + unmeth_count > 0;}
  size_t meth_count;
  size_t unmeth_count;
};


void
simulate_island_meth_status(Runif &runif, const gsl_rng *rng, 
			    const size_t n_samples, const double alpha, 
			    const vector<pair<double, double> > &beta_params, 
			    double &lower, double &upper, double &meth_freq) {
  
  static const double BETA_DISTRO_PARAM_1_PRIOR = 0.5;
  static const double BETA_DISTRO_PARAM_2_PRIOR = 0.5;
  
  const size_t n_cpgs = beta_params.size();
  vector<double> outcomes(n_samples);
  for (size_t i = 0; i < n_samples; ++i) {
    size_t sum = 0;
    for (size_t j = 0; j < n_cpgs; ++j) {
      const double p_sample = gsl_ran_beta(rng, beta_params[j].first + 
					   BETA_DISTRO_PARAM_1_PRIOR,
					   beta_params[j].second + 
					   BETA_DISTRO_PARAM_2_PRIOR);
      if (runif.runif(0.0, 1.0) < p_sample)
	++sum;
    }
    outcomes[i] = static_cast<double>(sum)/n_cpgs;
  }
  sort(outcomes.begin(), outcomes.end());

  // RETURN VALUES (PASSED THROUGH REFERENCE PARAMETERS)
  lower = outcomes[static_cast<size_t>(alpha*n_samples)];
  upper = outcomes[static_cast<size_t>((1 - alpha)*n_samples)];
  meth_freq = accumulate(outcomes.begin(), outcomes.end(), 0.0)/n_samples;
}


void
call_island(Runif &runif, const gsl_rng *rng, const size_t n_samples,
	    const double required_cpgs_with_values,
	    const double critical_value, const double alpha, 
	    const double ci_width, const vector<GenomicRegion> &cpgs, 
	    GenomicRegion &r) {
  
  static const size_t NO_CALL = 3ul;
  static const size_t PARTIAL_METH = 2ul;
  static const size_t METHYLATED = 1ul;
  static const size_t UNMETHYLATED = 0ul;

  vector<pair<double, double> > beta_params;
  for (size_t i = 0; i < cpgs.size(); ++i) {
    const CPGInfo cgi(cpgs[i]);
    if (cgi.has_values())
      beta_params.push_back(make_pair(cgi.meth_count, cgi.unmeth_count));
  }
  
  size_t meth_state = NO_CALL;
  bool confident_call = false;
  bool confident_comparison = false;
  
  double lower = 0;
  double upper = 0;
  double meth_freq = 0;
  
  if (beta_params.size() > required_cpgs_with_values*cpgs.size()) {
    
    confident_comparison = true;
    
    simulate_island_meth_status(runif, rng, n_samples, 
				alpha, beta_params, 
				lower, upper, meth_freq);
    
    const double tail_size = (1.0 - critical_value);
    confident_call = ((upper - lower) < ci_width);
    assert(upper <= 1.0 && lower >= 0.0);
    
    if (confident_call)
      meth_state = PARTIAL_METH;
    
    if (upper < tail_size)
      meth_state = UNMETHYLATED;
    else if (lower > critical_value)
      meth_state = METHYLATED;
    
  }

  const string annotated_name = r.get_name() + ":" +
    toa(meth_freq) + ":" + toa(upper) + ":" + toa(lower) + ":" + 
    toa(confident_comparison) + ":" + toa(confident_call);

  r.set_name(annotated_name);
  r.set_score(meth_state);
}

int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    
    string regions_file;
    string outfile;
    
    double alpha = 0.1;
    double crit = 0.75;
    double interval_width = 0.25;
    
    size_t n_samples = 10000ul;
    double required_cpgs_with_values = 0.9;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("cgicaller", "A program for calling CGI methylation "
			   "status from a Solexa bisulfite capture experiment");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("crit", 'C', "critical value", 
		      false , crit);
    opt_parse.add_opt("width", 'w', "interval width for confident call", 
		      false , interval_width);
    opt_parse.add_opt("alpha", 'a', "alpha value", false , alpha);
    opt_parse.add_opt("samples", 's', "random samples to take", false, n_samples);
    opt_parse.add_opt("reqcpg", 'q', "required cpgs with values", false, 
		      required_cpgs_with_values);
    opt_parse.add_opt("regions", 'r', "file of target regions", 
		      true , regions_file);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
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
    const string cpgs_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    srand(time(0) + getpid());
    gsl_rng_set(rng, rand());

    Runif runif(rand());
    
    if (VERBOSE)
      cerr << "reading cpg locations" << endl;
    vector<GenomicRegion> cpgs;
    ReadBEDFile(cpgs_file, cpgs);
    assert(check_sorted(cpgs));
    
    if (VERBOSE)
      cerr << "reading target regions" << endl;
    vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);
    assert(check_sorted(regions));
    
    if (VERBOSE)
      cerr << "separating regions" << endl;
    vector<vector<GenomicRegion> > sep_cpgs;
    separate_cpgs(cpgs, regions, sep_cpgs);
    
    if (VERBOSE)
      cerr << "calling island states" << endl;

    std::ostream *out = (outfile.empty()) ? &cout : 
      new std::ofstream(outfile.c_str());
    for (size_t i = 0; i < regions.size(); ++i)
      if (!sep_cpgs[i].empty()) {
	call_island(runif, rng, n_samples, required_cpgs_with_values,
		    crit, alpha, interval_width, sep_cpgs[i], regions[i]);
	*out << regions[i] << endl;
      }
    if (out != &cout) delete out;
  }
  catch (const RMAPException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
