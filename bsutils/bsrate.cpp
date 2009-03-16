/*    bsrate: a program for determining the rate of bisulfite
 *    conversion in a bisulfite sequencing experiment
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

template <class U, class V> static void
separate_regions(const std::vector<U> &regions, 
		 const std::vector<V> &seqs, 
		 std::vector<std::vector<U> > &sep_regions,
		 std::vector<std::vector<V> > &sep_seqs) {
  const size_t n_regions = regions.size();
  assert(n_regions <= seqs.size());
  if (regions.empty()) return;
  
  string prev_chrom = regions.front().get_chrom();
  sep_regions.resize(1);
  sep_seqs.resize(1);
  for (size_t i = 0; i < regions.size(); ++i) {
    std::string current_chrom(regions[i].get_chrom());
    if (current_chrom != prev_chrom) {
      prev_chrom.swap(current_chrom);
      sep_regions.push_back(vector<U>());
      sep_seqs.push_back(vector<V>());
    }
    sep_regions.back().push_back(regions[i]);
    sep_seqs.back().push_back(seqs[i]);
  }
}

static void
relative_sort(const vector<GenomicRegion> &mapped_locations, 
	      const vector<string> &names,
	      vector<size_t> &lookup) {
  
  unordered_map<string, size_t> names_map;
  for (size_t i = 0; i < names.size(); ++i)
    names_map[names[i]] = i;
  
  for (size_t i = 0; i < mapped_locations.size(); ++i) {
    const unordered_map<string, size_t>::const_iterator 
      j(names_map.find(mapped_locations[i].get_name()));
    if (j == names_map.end())
      throw RMAPException("read sequence not found for: " + names[i]);
    lookup.push_back(j->second);
  }
}


// void
// separate_cpgs(const vector<GenomicRegion> &cpgs, 
// 	      const vector<GenomicRegion> &regions,
// 	      vector<vector<GenomicRegion> > &sep_cpgs) {
  
//   unordered_map<string, size_t> region_id_lookup;
//   for (size_t i = 0; i < regions.size(); ++i)
//     region_id_lookup[regions[i].get_name()] = i;
  
//   sep_cpgs.resize(regions.size());
//   for (size_t i = 0; i < cpgs.size(); ++i) {
//     const string name(cpgs[i].get_name());
//     const size_t region_name_end = name.find_first_of(":");
//     const string region_name(name.substr(0, region_name_end));
//     const unordered_map<string, size_t>::const_iterator j =
//       region_id_lookup.find(region_name);
//     if (j == region_id_lookup.end()) 
//       cerr << region_name << endl;
//     assert(j != region_id_lookup.end());
//     sep_cpgs[j->second].push_back(cpgs[i]);
//   }
// }

// struct CPGInfo {
//   CPGInfo(const GenomicRegion &r);
//   bool has_values() const {return meth_count + unmeth_count > 0;}
//   size_t meth_count;
//   size_t unmeth_count;
// };

// CPGInfo::CPGInfo(const GenomicRegion &r) {
//   const string name(r.get_name());
//   vector<string> parts = rmap::split(name, ":");
//   unmeth_count = atoi(parts[1].c_str());
//   meth_count = atoi(parts[2].c_str());
// }

// void
// simulate_island_meth_status(Runif &runif, const gsl_rng *rng, 
// 			    const size_t n_samples, const double alpha, 
// 			    const vector<pair<double, double> > &beta_params, 
// 			    double &lower, double &upper, double &meth_freq) {
//   const size_t n_cpgs = beta_params.size();
//   vector<double> outcomes(n_samples);
//   for (size_t i = 0; i < n_samples; ++i) {
//     size_t sum = 0;
//     for (size_t j = 0; j < n_cpgs; ++j) {
//       const double p_sample = gsl_ran_beta(rng, beta_params[j].first + 0.5, 
// 					   beta_params[j].second + 0.5);
//       if (runif.runif(0.0, 1.0) < p_sample)
// 	++sum;
//     }
//     outcomes[i] = static_cast<double>(sum)/n_cpgs;
//   }
  
//   sort(outcomes.begin(), outcomes.end());
//   lower = outcomes[alpha*n_samples];
//   upper = outcomes[(1 - alpha)*n_samples];
//   meth_freq = accumulate(outcomes.begin(), outcomes.end(), 0.0)/n_samples;
// }


// void
// call_island(Runif &runif,
// 	    const gsl_rng *rng,
// 	    const size_t n_samples,
// 	    const double required_cpgs_with_values,
// 	    const double meth_unmeth_island_critical_value,
// 	    const double meth_unmeth_island_alpha_value,
// 	    const double confidence_interval_width,
// 	    const vector<GenomicRegion> &cpgs, GenomicRegion &r) {


//   vector<pair<double, double> > beta_params;
//   for (size_t i = 0; i < cpgs.size(); ++i) {
//     const CPGInfo cgi(cpgs[i]);
//     if (cgi.has_values())
//       beta_params.push_back(make_pair(cgi.meth_count, cgi.unmeth_count));
//   }
  
//   size_t meth_state = 3ul;
//   bool confident_call = false;
//   bool confident_comparison = false;

//   double lower = 0;
//   double upper = 0;
//   double meth_freq = 0;
  
//   if (beta_params.size() > required_cpgs_with_values*cpgs.size()) {

//     meth_state = 2ul;
    
//     simulate_island_meth_status(runif, rng, n_samples, 
// 				meth_unmeth_island_alpha_value, beta_params, 
// 				lower, upper, meth_freq);
    
//     const double tail_size = (1.0 - meth_unmeth_island_critical_value);
//     confident_call = ((upper - lower) < confidence_interval_width);
//     assert(upper <= 1.0 && lower >= 0.0);
    
//     if (confident_call)
//       meth_state = 2ul;
    
//     if (upper < tail_size)
//       meth_state = 0ul;
//     else if (lower > meth_unmeth_island_critical_value)
//       meth_state = 1ul;
    
//     confident_comparison = true;
//   }
//   const string annotated_name = r.get_name() + ":" +
//     toa(meth_freq) + ":" + toa(upper) + ":" + toa(lower) + ":" + 
//     toa(confident_comparison) + ":" + toa(confident_call);
//   r.set_name(annotated_name);
//   r.set_score(meth_state);
// }

static bool 
is_cytosine(char c) {return (c == 'c' || c == 'C');}

static bool 
is_guanine(char c)  {return (c == 'g' || c == 'G');}

static bool 
is_thymine(char c)  {return (c == 't' || c == 'T');}

static void
check_conversionfun(const string &chrom, 
		    const vector<GenomicRegion> locations,
		    const vector<string> sequences,
		    size_t &total_pos, size_t &conv_pos, 
		    size_t &total_neg, size_t &conv_neg) {

  const size_t read_width = locations.front().get_width();
  const size_t chrom_size = chrom.length();

  for (size_t i = 0; i < locations.size(); ++i) {
    assert(locations[i].get_end() < chrom_size);
    const size_t offset = locations[i].get_start();
    if (locations[i].pos_strand()) {
      for (size_t j = 0; j < read_width - 1; ++j) 
	if (offset + j + 1 < chrom_size) {
	  if (is_cytosine(chrom[offset + j]) && !is_guanine(chrom[offset + j + 1])) {
	    total_pos += (is_cytosine(sequences[i][j]) || is_thymine(sequences[i][j]));
	    conv_pos += is_thymine(sequences[i][j]);
	  }
	}
    }
    else {
      for (size_t j = 0; j < read_width; ++j)
	if (offset + read_width >= j + 2) {
	  if (is_guanine(chrom[offset + read_width - 1 - j]) && 
	      !is_cytosine(chrom[offset + read_width - 2 - j])) {
	    total_neg += (is_cytosine(sequences[i][j]) || is_thymine(sequences[i][j]));
	    conv_neg += is_thymine(sequences[i][j]);
	  }
	}
    }
  }
}

int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    
    string chrom_dir;
    string mapped_reads_file;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("bsrate", "A program for obtaining the conversion "
			   "rate in a bisulfite sequencing experiment");
    opt_parse.add_opt("chrom", 'c', "chromosome directory (FASTA format)", 
		      true , chrom_dir);
    opt_parse.add_opt("mapped", 'm', "file of mapped reads (BED format)", 
		      true , mapped_reads_file);
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
    const string reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "reading mapped read locations" << endl;
    vector<GenomicRegion> mapped_reads;
    ReadBEDFile(mapped_reads_file, mapped_reads);
    assert(check_sorted(mapped_reads));
    
    if (VERBOSE)
      cerr << "loading read sequences" << endl;
    vector<string> names, sequences;
    read_fasta_file(reads_file.c_str(), names, sequences);
    if (VERBOSE)
      cerr << "read " << names.size() << " sequences" << endl;

    if (VERBOSE)
      cerr << "associating read sequences to mapped locations" << endl;
    vector<size_t> lookup;
    relative_sort(mapped_reads, names, lookup);
    vector<string> seq_swapper(sequences.size());
    for (size_t i = 0; i < lookup.size(); ++i)
      seq_swapper[i].swap(sequences[lookup[i]]);
    sequences.swap(seq_swapper);
    seq_swapper.clear();

    if (VERBOSE)
      cerr << "separating chromosomes" << endl;
    vector<vector<GenomicRegion> > mapped_reads_by_chrom;
    vector<vector<string> > sequences_by_chrom;
    separate_regions(mapped_reads, sequences, 
		     mapped_reads_by_chrom, sequences_by_chrom);
    mapped_reads.clear();
    sequences.clear();

    size_t total_pos = 0, conv_pos = 0;
    size_t total_neg = 0, conv_neg = 0;
    
    for (size_t i = 0; i < mapped_reads_by_chrom.size(); ++i) {
      const string chrom(mapped_reads_by_chrom[i].front().get_chrom());
      const string chrom_file(path_join(chrom_dir, chrom + ".fa"));
      
      vector<string> dummy_chrom_names, chrom_seqs;
      if (VERBOSE)
	cerr << "[PROCESSING] " << chrom_file << endl;
      read_fasta_file(chrom_file.c_str(), dummy_chrom_names, chrom_seqs);
      
      check_conversionfun(chrom_seqs.front(), 
			  mapped_reads_by_chrom[i], sequences_by_chrom[i],
			  total_pos, conv_pos, total_neg, conv_neg);
    }
    
    cout << "\tTOTAL\tCONVERTED\tRATE" << endl;
    cout << "P:\t" << total_pos << "\t" << conv_pos << "\t" << conv_pos/static_cast<double>(total_pos) << endl;
    cout << "N:\t" << total_neg << "\t" << conv_neg << "\t" << conv_neg/static_cast<double>(total_neg) << endl;
    cout << "BOTH:\t" << total_pos + total_neg << "\t" 
	 << conv_pos + conv_neg << "\t" 
	 << (conv_pos + conv_neg)/static_cast<double>(total_pos + total_neg) << endl;
    
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
