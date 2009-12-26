/* methdiff: A program for determining the probability that
 * methylation at each CpG (or any nucleotide) differs 
 * between two conditions.
 *
 * Copyright (C) 2009 University of Southern California
 *                    Andrew D Smith
 * Author: Andrew D. Smith
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"

#include <gsl/gsl_sf_gamma.h>

#include <cmath>
#include <fstream>
#include <algorithm>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;

using std::ostream_iterator;
using std::ofstream;

static inline double
log_sum_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}

static void
get_meth_unmeth(const GenomicRegion &cpg, size_t &meth, size_t &unmeth) {
  const double prob = cpg.get_score();
  const string name(cpg.get_name());
  const size_t n_reads = atoi(name.substr(name.find_first_of(":") + 1).c_str());
  meth = prob*n_reads;
  unmeth = n_reads - meth;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

static double
log_hyper_g_greater(size_t meth_a, size_t unmeth_a, 
		    size_t meth_b, size_t unmeth_b, size_t k) {
  return  gsl_sf_lnchoose(meth_b + unmeth_b - 1, k) + 
    gsl_sf_lnchoose(meth_a + unmeth_a - 1, meth_a + meth_b - 1 - k) -
    gsl_sf_lnchoose(meth_a + unmeth_a + meth_b + unmeth_b - 2, 
		    meth_a + meth_b - 1);
}
  

static double
test_greater_population(size_t meth_a, size_t unmeth_a, 
			size_t meth_b, size_t unmeth_b) {
  double p = 0;
  
  for (size_t k = (meth_b > unmeth_a) ? meth_b - unmeth_a : 0; k < meth_b; ++k)
    p = log_sum_log(p, log_hyper_g_greater(meth_a, unmeth_a, meth_b, unmeth_b, k));
  return exp(p);
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// static double
// log_hyper_g(const size_t a, const size_t c, 
// 	    const size_t b, const size_t d) {
//   return  gsl_sf_lnfact(a + b) + gsl_sf_lnfact(c + d) + 
//     gsl_sf_lnfact(a + c) + gsl_sf_lnfact(b + d) -
//     gsl_sf_lnfact(a + b + c + d) - gsl_sf_lnfact(a) - 
//     gsl_sf_lnfact(b) - gsl_sf_lnfact(c) - gsl_sf_lnfact(d);
// }

// static double
// test_similar_population(size_t meth_a, size_t unmeth_a, 
// 			size_t meth_b, size_t unmeth_b) {
//   double p = 0;

// //   size_t a = 0, b = 0, c = 0, d = 0;
// //   if ((meth_a + 1)*(unmeth_b + 1) > unmeth_a*meth_b) {
// //     a = meth_a;   c = meth_b; 
// //     b = unmeth_a; d = unmeth_b;
// //   }
// //   else {
// //     c = meth_a;   a = meth_b;   
// //     d = unmeth_a; b = unmeth_b; 
// //   }
// //   while (b > 0 && c > 0) {
// //     p = log_sum_log(p, log_hyper_g(a, c, b, d));
// //     ++a; --b;
// //     ++d; --c;
// //   }
// //   p = log_sum_log(p, log_hyper_g(a, c, b, d));


//   while (unmeth_a > 0 && meth_b > 0) {
//     p = log_sum_log(p, log_hyper_g(meth_a, meth_b, unmeth_a, unmeth_b));
//     ++meth_a; --unmeth_a;
//     ++unmeth_b; --meth_b;
//   }
//   p = log_sum_log(p, log_hyper_g(meth_a, meth_b, unmeth_a, unmeth_b));
//   return exp(p);
// }


int
main(int argc, const char **argv) {

  try {

    string outfile;
    double pseudocount = 1.0;
    
    // run mode flags
    bool VERBOSE = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("methdiff", 
			   "A program for determining the probability that "
			   "methylation at each CpG (or any nucleotide) differs "
			   "between two conditions.",
			   "<cpgs_file_a> <cpgs_file_b>");
    opt_parse.add_opt("pseudo", 'p', "pseudocount (default: 1)", 
		      false, pseudocount);
    opt_parse.add_opt("out", 'o', "output file (BED format)", 
		      false, outfile);
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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string cpgs_file_a = leftover_args[0];
    const string cpgs_file_b = leftover_args[1];
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "[READING CPGS]";
    vector<GenomicRegion> cpgs_a, cpgs_b;
    ReadBEDFile(cpgs_file_a, cpgs_a);
    if (VERBOSE) cerr << "[READ=" + strip_path(cpgs_file_a) + "]";
    if (!check_sorted(cpgs_a))
      throw RMAPException("CpGs not sorted in file \"" + cpgs_file_a + "\"");
    if (VERBOSE) cerr << "[SORTED]";
    ReadBEDFile(cpgs_file_b, cpgs_b);
    if (VERBOSE) cerr << "[READ=" + strip_path(cpgs_file_b) + "]";
    if (!check_sorted(cpgs_b))
      throw RMAPException("CpGs not sorted in file \"" + cpgs_file_b + "\"");
    if (VERBOSE)
      cerr << "[SORTED]"
	   << "[DONE]" << endl
	   << "CPG COUNT A: " << cpgs_a.size() << endl
	   << "CPG COUNT B: " << cpgs_b.size() << endl;
    
    std::ostream *out = (outfile.empty()) ? &cout : 
      new std::ofstream(outfile.c_str());
    size_t j = 0;
    for (size_t i = 0; i < cpgs_a.size(); ++i) {
      if (VERBOSE && (i == 0 || !cpgs_a[i - 1].same_chrom(cpgs_a[i])))
	cerr << "processing " << cpgs_a[i].get_chrom() << endl;
      
      while (j < cpgs_b.size() && cpgs_b[j] < cpgs_a[i]) ++j;
      
      if (cpgs_a[i].same_chrom(cpgs_b[j]) && 
	  cpgs_a[i].get_start() == cpgs_b[j].get_start()) {
	
	size_t meth_a = 0, unmeth_a = 0;
	get_meth_unmeth(cpgs_a[i], meth_a, unmeth_a);
	
	size_t meth_b = 0, unmeth_b = 0;
	get_meth_unmeth(cpgs_b[j], meth_b, unmeth_b);
	
	meth_a += pseudocount;
	meth_b += pseudocount;
	unmeth_a += pseudocount;
	unmeth_b += pseudocount;
	
 	cpgs_a[i].set_score(test_greater_population(meth_b, unmeth_b, 
 						    meth_a, unmeth_a));
	
	cpgs_a[i].set_name(std::min(cpgs_a[i].get_name(), cpgs_b[j].get_name()));

	*out << cpgs_a[i] << endl;
      }
    }
    if (out != &cout) delete out;
  }
  catch (RMAPException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
