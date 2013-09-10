/*    amrfinder: A program for resolving epialleles in a sliding
 *    window along a chromosome.
 *
 *    Copyright (C) 2011-2013 University of Southern California and
 *                            Andrew D. Smith and Fang Fang
 *
 *    Authors: Fang Fang and Andrew D. Smith
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
#include <iostream>
#include <numeric>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "EpireadStats.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::streampos;
using std::tr1::unordered_map;



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


static double
get_fdr_cutoff(const size_t n_tests, const vector<GenomicRegion> &amrs, 
	       const double fdr) {
  if (fdr <= 0) return std::numeric_limits<double>::max();
  else if (fdr > 1) return std::numeric_limits<double>::min();
  
  vector<double> pvals;
  for (size_t i = 0; i < amrs.size(); ++i)
    pvals.push_back(amrs[i].get_score());
  std::sort(pvals.begin(), pvals.end());
  
  const size_t n_pvals = pvals.size();
  size_t i = 0;
  while (i < n_pvals - 1 && pvals[i+1] < fdr*static_cast<double>(i+1)/n_tests)
    ++i;
  return pvals[i];
}


static void
clip_read(const size_t start_pos, const size_t end_pos, epiread &r) {
  if (r.pos < start_pos) {
    r.seq = r.seq.substr(start_pos - r.pos);
    r.pos = start_pos;
  }
  if (r.end() > end_pos)
    r.seq = r.seq.substr(0, end_pos - r.pos);
}


static void
get_current_epireads(const vector<epiread> &epireads, 
		     const size_t cpg_window, const size_t start_pos, 
		     size_t &read_id, 
		     vector<epiread> &current_epireads) {
  
  const size_t end_pos = start_pos + cpg_window;
  for (size_t i = read_id; i < epireads.size() && epireads[i].pos < end_pos; ++i)
    if (epireads[i].end() > start_pos) {
      if (current_epireads.empty())
	read_id = i;
      current_epireads.push_back(epireads[i]);
      clip_read(start_pos, end_pos, current_epireads.back());
    }
}


static void
eliminate_amrs_by_fdr(const double fdr_cutoff, vector<GenomicRegion> &amrs) {
  size_t j = 0;
  for (size_t i = 0; i < amrs.size(); ++i)
    if (amrs[i].get_score() < fdr_cutoff) {
      amrs[j] = amrs[i];
      ++j;
    }
  amrs.erase(amrs.begin() + j, amrs.end());
}


static void
eliminate_amrs_by_size(const size_t min_size, vector<GenomicRegion> &amrs) {
  size_t j = 0;
  for (size_t i = 0; i < amrs.size(); ++i)
    if (amrs[i].get_width() >= min_size) {
      amrs[j] = amrs[i];
      ++j;
    }
  amrs.erase(amrs.begin() + j, amrs.end());
}


static void
collapse_amrs(vector<GenomicRegion> &amrs) {
  size_t j = 0;
  for (size_t i = 1; i < amrs.size(); ++i)
    if (amrs[j].same_chrom(amrs[i]) &&
	// The +1 below is because intervals in terms of CpGs are
	// inclusive
        amrs[j].get_end() + 1 >= amrs[i].get_start()) {
      amrs[j].set_end(amrs[i].get_end());
      amrs[j].set_score(std::min(amrs[i].get_score(), amrs[j].get_score()));
    }
    else {
      ++j;
      amrs[j] = amrs[i];
    }
  ++j;
  amrs.erase(amrs.begin() + j, amrs.end());
}


static string
get_amr_name(const size_t x, const size_t y) {
  static const string name_label("AMR");
  return name_label + toa(x) + ":" + toa(y);
}


static void
add_amr(const string &chrom_name, const size_t start_cpg, 
	const size_t cpg_window, const vector<epiread> &reads, 
	const double score, vector<GenomicRegion> &amrs) {
  const size_t end_cpg = start_cpg + cpg_window - 1;
  const string amr_name(get_amr_name(amrs.size(), reads.size()));
  amrs.push_back(GenomicRegion(chrom_name, start_cpg, end_cpg,
			       amr_name, score, '+'));
}


// ////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////
// /////
// /////   CODE FOR RANDOMIZING THE READS TO GET EXPECTED NUMBER OF
// /////   IDENTIFIED AMRS
// /////

// static void
// set_read_states(vector<vector<char> > &state_counts, vector<epiread> &reads) {
//   for (size_t i = 0; i < reads.size(); ++i) {
//     const size_t offset = reads[i].pos;
//     for (size_t j = 0; j < reads[i].length(); ++j) {
//       reads[i].seq[j] = state_counts[offset + j].back();
//       state_counts[offset + j].pop_back();
//     }
//   }
// }

// static void
// get_state_counts(const vector<epiread> &reads, const size_t total_cpgs, 
// 		 vector<vector<char> > &state_counts) {
  
//   state_counts = vector<vector<char> >(total_cpgs);
//   for (size_t i = 0; i < reads.size(); ++i) {
//     const size_t offset = reads[i].pos;
//     for (size_t j = 0; j < reads[i].length(); ++j)
//       state_counts[offset + j].push_back(reads[i].seq[j]);
//   }
//   for (size_t i = 0; i < state_counts.size(); ++i)
//     random_shuffle(state_counts[i].begin(), state_counts[i].end());
// }

// static void
// randomize_read_states(vector<epiread> &reads) {
//   srand(time(0) + getpid());
//   const size_t total_cpgs = get_n_cpgs(reads);
//   vector<vector<char> > state_counts;
//   get_state_counts(reads, total_cpgs, state_counts);
//   set_read_states(state_counts, reads);
// }


static void
merge_amrs(const size_t gap_limit, vector<GenomicRegion> &amrs) {
  size_t j = 0;
  for (size_t i = 1; i < amrs.size(); ++i)
    // check distance between two amrs is greater than gap limit
    if (amrs[j].same_chrom(amrs[i]) &&
        amrs[j].get_end() + gap_limit>= amrs[i].get_start()) {
      amrs[j].set_end(amrs[i].get_end());
      amrs[j].set_score(std::min(amrs[i].get_score(), amrs[j].get_score()));
    }
    else {
      ++j;
      amrs[j] = amrs[i];
    }
  ++j;
  amrs.erase(amrs.begin() + j, amrs.end());
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


inline static bool
is_cpg(const string &s, const size_t idx) {
  return toupper(s[idx]) == 'C' && toupper(s[idx + 1]) == 'G';
}


static void
collect_cpgs(const string &s, unordered_map<size_t, size_t> &cpgs) {
  const size_t lim = s.length() - 1;
  size_t cpg_count = 0;
  for (size_t i = 0; i < lim; ++i)
    if (is_cpg(s, i)) {
      cpgs[cpg_count] = i;
      ++cpg_count;
    }
}


static void
identify_chromosomes(const string chrom_file, const string fasta_suffix, 
                     unordered_map<string, string> &chrom_files) {
  vector<string> the_files;
  if (isdir(chrom_file.c_str())) {
    read_dir(chrom_file, fasta_suffix, the_files);
    for (size_t i = 0; i < the_files.size(); ++i)
      chrom_files[strip_path_and_suffix(the_files[i])] = the_files[i];
  }
  else chrom_files[strip_path_and_suffix(chrom_file)] = chrom_file;
}


static void
convert_coordinates(const unordered_map<size_t, size_t> &cpgs, 
                    GenomicRegion &region)  {
  const unordered_map<size_t, size_t>::const_iterator 
    start_itr(cpgs.find(region.get_start()));
  const unordered_map<size_t, size_t>::const_iterator
    end_itr(cpgs.find(region.get_end()));
  if (start_itr == cpgs.end() || end_itr == cpgs.end())
    throw SMITHLABException("could not convert:\n" + region.tostring());
  region.set_start(start_itr->second);
  region.set_end(end_itr->second);
}


static void
convert_coordinates(const bool VERBOSE, const string chroms_dir,
                    const string fasta_suffix, vector<GenomicRegion> &amrs) {
  
  unordered_map<string, string> chrom_files;
  identify_chromosomes(chroms_dir, fasta_suffix, chrom_files);
  if (VERBOSE)
    cerr << "CHROMS:\t" << chrom_files.size() << endl;
  
  unordered_map<size_t, size_t> cpgs;
  vector<string> chrom_names, chroms;
  GenomicRegion region;
  GenomicRegion chrom_region("chr0", 0, 0);
  for (size_t i = 0; i < amrs.size(); ++i) {
    
    // get the correct chrom if it has changed
    if (!amrs[i].same_chrom(chrom_region)) {
      const unordered_map<string, string>::const_iterator 
        fn(chrom_files.find(amrs[i].get_chrom()));
      if (fn == chrom_files.end())
        throw SMITHLABException("could not find chrom: " + amrs[i].get_chrom());
      
      chrom_names.clear();
      chroms.clear();
      read_fasta_file(fn->second.c_str(), chrom_names, chroms);
      if (chrom_names.size() > 1)
        throw SMITHLABException("multiple chroms/file: " + fn->second);
      if (VERBOSE)
        cerr << "CONVERTING: " << chrom_names.front() << endl;
      collect_cpgs(chroms.front(), cpgs);
      chrom_region.set_chrom(chrom_names.front());
    }
    
    convert_coordinates(cpgs, amrs[i]);
  }
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


static size_t
total_states(const vector<epiread> &epireads) {
  size_t total = 0;
  for (size_t i = 0; i < epireads.size(); ++i)
    total += epireads[i].length();
  return total;
}

static size_t
process_chrom(const bool VERBOSE, const bool PROGRESS,
	      const size_t min_obs_per_cpg, const size_t window_size,
	      const EpireadStats &epistat, const string &chrom_name, 
	      const vector<epiread> &epireads, vector<GenomicRegion> &amrs) {
  
  const size_t min_obs_per_window = window_size*min_obs_per_cpg;
  
  const size_t chrom_cpgs = get_n_cpgs(epireads);
  if (VERBOSE)
    cerr << "PROCESSING: " << chrom_name << " "
	 << "[reads: " << epireads.size() << "] "
	 << "[cpgs: " << chrom_cpgs << "]" << endl;
  
  const size_t PROGRESS_TIMING_MODULUS = std::max(1ul, epireads.size()/1000);

  size_t windows_tested = 0;
  
  size_t start_idx = 0;
  const size_t lim = chrom_cpgs - window_size + 1;
  for (size_t i = 0; i < lim && start_idx < epireads.size(); ++i) {
    
    if (PROGRESS && i % PROGRESS_TIMING_MODULUS == 0) 
      cerr << '\r' << chrom_name << ' ' << percent(i, chrom_cpgs) << "%\r";
    
    vector<epiread> current_epireads;
    get_current_epireads(epireads, window_size, i, start_idx, current_epireads);
    
    if (total_states(current_epireads) > min_obs_per_window) {
      bool is_significant = false;
      const double score = epistat.test_asm(current_epireads, is_significant);
      if (is_significant)
	add_amr(chrom_name, i, window_size, current_epireads, score, amrs);
      ++windows_tested;
    }
  }
  if (PROGRESS)
    cerr << '\r' << chrom_name << " 100%" << endl;
  return windows_tested;
}


int 
main(int argc, const char **argv) {
  
  try {

    static const string fasta_suffix = "fa";
    
    bool VERBOSE = false;
    bool PROGRESS = false;

    string outfile;
    
    size_t max_itr = 10;
    size_t window_size = 10;
    size_t gap_limit = 1000;
    
    double high_prob = 0.75, low_prob = 0.25;
    double min_obs_per_cpg = 4;
    double critical_value = 0.01;
    
    // bool RANDOMIZE_READS = false;
    bool IGNORE_BALANCE = false;
    bool USE_BIC = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), 
			   "identify regions of allele-specific methylation", 
			   "<chroms-dir> <epireads>");
    opt_parse.add_opt("outfile", 'o', "output file", false, outfile);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_itr);
    opt_parse.add_opt("no-bal", 'u', "no penalty for unbalanced alleles", 
		      false, IGNORE_BALANCE);
    opt_parse.add_opt("window", 'w', "size of sliding window", 
		      false, window_size);
    opt_parse.add_opt("min-cov", 'm', "min coverage per cpg to test windows", 
		      false, min_obs_per_cpg);
    opt_parse.add_opt("gap", 'g', "min allowed gap between amrs (in bp)", 
    		      false, gap_limit);
    opt_parse.add_opt("crit", 'c', "critical p-value cutoff (default: 0.01)", 
		      false, critical_value);
    opt_parse.add_opt("bic", 'b', "use BIC to compare models", false, USE_BIC);
    // opt_parse.add_opt("rand", 'R', "randomize reads (for comparison)", 
    // 		      false, RANDOMIZE_READS);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("progress", 'P', "print progress info", false, PROGRESS);
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
    const string chroms_dir(leftover_args.front());
    const string reads_file(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/
    
    if (VERBOSE)
      cerr << "AMR TESTING OPTIONS: "
	   << "[test=" << (USE_BIC ? "BIC" : "LRT") << "] "
	   << "[balance=" << !IGNORE_BALANCE << "] "
	   << "[iterations=" << max_itr << "]" << endl;
    
    const EpireadStats epistat(low_prob, high_prob, critical_value, max_itr,
			       IGNORE_BALANCE, USE_BIC);
    
    std::ifstream in(reads_file.c_str());
    if (!in)
      throw SMITHLABException("cannot open input file: " + reads_file);
    
    vector<GenomicRegion> amrs;
    size_t windows_tested = 0;
    
    vector<epiread> epireads;
    string prev_chrom, curr_chrom, tmp_states;
    size_t tmp_pos;
    while (in >> curr_chrom >> tmp_pos >> tmp_states) {
      if (!epireads.empty() && curr_chrom != prev_chrom) {
	windows_tested += 
	  process_chrom(VERBOSE, PROGRESS, min_obs_per_cpg, window_size,
			epistat, prev_chrom, epireads, amrs);
	epireads.clear();
      }
      epireads.push_back(epiread(tmp_pos, tmp_states));
      std::swap(prev_chrom, curr_chrom);
    }
    if (!epireads.empty()) {
      windows_tested += 
	process_chrom(VERBOSE, PROGRESS, min_obs_per_cpg, window_size,
		      epistat, prev_chrom, epireads, amrs);
    }
    
    //////////////////////////////////////////////////////////////////
    //////  POSTPROCESSING IDENTIFIED AMRS AND COMPUTING SUMMARY STATS
    if (VERBOSE)
      cerr << endl << "========= POST PROCESSING =========" << endl;
    
    const size_t windows_accepted = amrs.size();
    const double fdr_cutoff = (USE_BIC) ? 0.0 :
      get_fdr_cutoff(windows_tested, amrs, critical_value);
    
    collapse_amrs(amrs);
    const size_t collapsed_amrs = amrs.size();
    convert_coordinates(VERBOSE, chroms_dir, fasta_suffix, amrs);
    merge_amrs(gap_limit, amrs);
    const size_t merged_amrs = amrs.size();
    
    eliminate_amrs_by_fdr(fdr_cutoff, amrs);
    const size_t amrs_passing_fdr = amrs.size();
    
    eliminate_amrs_by_size(gap_limit/2, amrs);
    
    if (VERBOSE)
      cerr << "WINDOWS TESTED: " << windows_tested << endl
	   << "WINDOWS ACCEPTED: " << windows_accepted << endl
     	   << "COLLAPSED WINDOWS: " << collapsed_amrs << endl
	   << "MERGED WINDOWS: " << merged_amrs << endl
	   << "FDR CUTOFF: " << fdr_cutoff << endl
     	   << "WINDOWS PASSING FDR: " << amrs_passing_fdr << endl
     	   << "AMRS (WINDOWS PASSING FDR: " << amrs.size() << endl;
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    copy(amrs.begin(), amrs.end(), 
	 std::ostream_iterator<GenomicRegion>(out, "\n"));
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
