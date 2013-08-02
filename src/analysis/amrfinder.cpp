/*    amrfinder: A program for resolving epialleles in a sliding
 *    window along a chromosome.
 *
 *    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith and Fang Fang
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
#include <iomanip>
#include <numeric>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "EpireadStats.hpp"
#include "EpireadIO.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;


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
get_current_window_reads(const vector<epiread> &all_reads_for_chrom, 
			 const size_t cpg_window, const size_t start_pos, 
			 size_t &read_id, 
			 vector<epiread> &current_window_reads) {
  
  const size_t end_pos = start_pos + cpg_window;
  for (size_t i = read_id; i < all_reads_for_chrom.size() && 
	 all_reads_for_chrom[i].pos < end_pos; ++i)
    if (all_reads_for_chrom[i].end() > start_pos) {
      if (current_window_reads.empty())
	read_id = i;
      current_window_reads.push_back(all_reads_for_chrom[i]);
      clip_read(start_pos, end_pos, current_window_reads.back());
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
collapse_amrs(vector<GenomicRegion> &amrs) {
  size_t j = 0;
  for (size_t i = 1; i < amrs.size(); ++i)
    if (amrs[j].same_chrom(amrs[i]) &&
	// The +1 below is because intervals in terms of CpGs are
	// inclusive
        amrs[j].get_end() + 1>= amrs[i].get_start()) {
      amrs[j].set_end(amrs[i].get_end());
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

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////
/////   CODE FOR RANDOMIZING THE READS TO GET EXPECTED NUMBER OF
/////   IDENTIFIED AMRS
/////

static void
set_read_states(vector<vector<char> > &state_counts, vector<epiread> &reads) {
  for (size_t i = 0; i < reads.size(); ++i) {
    const size_t offset = reads[i].pos;
    for (size_t j = 0; j < reads[i].length(); ++j) {
      reads[i].seq[j] = state_counts[offset + j].back();
      state_counts[offset + j].pop_back();
    }
  }
}

static void
get_state_counts(const vector<epiread> &reads, const size_t total_cpgs, 
		 vector<vector<char> > &state_counts) {
  
  state_counts = vector<vector<char> >(total_cpgs);
  for (size_t i = 0; i < reads.size(); ++i) {
    const size_t offset = reads[i].pos;
    for (size_t j = 0; j < reads[i].length(); ++j)
      state_counts[offset + j].push_back(reads[i].seq[j]);
  }

  for (size_t i = 0; i < state_counts.size(); ++i)
    random_shuffle(state_counts[i].begin(), state_counts[i].end());
}

static void
randomize_read_states(vector<epiread> &reads) {
  srand(time(0) + getpid());
  const size_t total_cpgs = get_n_cpgs(reads);
  vector<vector<char> > state_counts;
  get_state_counts(reads, total_cpgs, state_counts);
  set_read_states(state_counts, reads);
}

static void
merge_amrs(vector<GenomicRegion> &amrs,
			const size_t gap_limit) {
  size_t j = 0;
  for (size_t i = 1; i < amrs.size(); ++i)
    if (amrs[j].same_chrom(amrs[i]) &&
	// check the distance between two amrs are greater than the gap limit
        amrs[j].get_end() + gap_limit>= amrs[i].get_start()) {
      amrs[j].set_end(amrs[i].get_end());
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
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    string outfile;
    string chroms_dir;
    bool EPIREAD_FORMAT = false;
    bool USE_BIC = false;

    size_t max_itr = 10;
    size_t cpg_window = 10;
    size_t gap_limit=1000;
    double high_prob = 0.75, low_prob = 0.25;

    double min_reads_per_cpg = 1ul;
    bool RANDOMIZE_READS = false;
    bool PROGRESS = false;

    double critical_value = 0.01;
 
    bool IGNORE_BALANCED_PARTITION_INFO = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), 
			   "resolve epi-alleles", "<reads-file>");
    opt_parse.add_opt("outfile", 'o', "output file", false, outfile);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_itr);
    opt_parse.add_opt("no-bal", 'g', "ignore balanced partition info", 
		      false, IGNORE_BALANCED_PARTITION_INFO);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("epiread", 'E', "reads in epiread format", 
		      false, EPIREAD_FORMAT);
    opt_parse.add_opt("window", 'w', "the window to slide", false, cpg_window);
    opt_parse.add_opt("min-reads", 'm', "min reads per cpg", 
		      false, min_reads_per_cpg);
    opt_parse.add_opt("gap_limit", 'l', "the minimum limit of the gap size between amrs in bp", 
    		      false, gap_limit);
    opt_parse.add_opt("chrom", 'c', "dir of chroms (.fa extn)", 
		      true, chroms_dir);
    opt_parse.add_opt("crit", 'C', "critical p-value cutoff (default: 0.01)", 
		      false, critical_value);
    opt_parse.add_opt("bic", 'b', "use BIC to compare models", false, USE_BIC);
    opt_parse.add_opt("rand", 'R', "randomize reads", false, RANDOMIZE_READS);
    opt_parse.add_opt("progress", 'P', "write progress info", false, PROGRESS);
    
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
      cerr << opt_parse.option_missing_message() << " directory of chromosomes."<< endl;
      return EXIT_SUCCESS;
    }
    const string reads_file_name(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    const size_t min_reads_per_window = cpg_window*min_reads_per_cpg;
    
    EpireadIO eio(reads_file_name, VERBOSE, EPIREAD_FORMAT, chroms_dir);
    
    vector<GenomicRegion> amrs;
    vector<epiread> all_reads_for_chrom;
    string chrom_name;
    
    size_t total_reads_processed = 0;
    size_t windows_tested = 0;
    size_t total_cpgs = 0;
    while (eio.load_reads_next_chrom(chrom_name, all_reads_for_chrom)) {
      total_reads_processed += all_reads_for_chrom.size();

      const size_t chrom_cpgs = get_n_cpgs(all_reads_for_chrom);
      total_cpgs += chrom_cpgs;
      if (VERBOSE)
	cerr << "PROCESSING: " << chrom_name << endl
	     << "INFORMATIVE READS: " << all_reads_for_chrom.size() << endl
	     << "CHROM CPGS: " << chrom_cpgs << endl;
      
      if (RANDOMIZE_READS)
	randomize_read_states(all_reads_for_chrom);
      
      const size_t PROGRESS_TIMING_MODULUS = 
	std::max(1ul, all_reads_for_chrom.size()/1000);
      
      size_t start_read_idx = 0;
      const size_t lim = chrom_cpgs - cpg_window + 1;
      for (size_t i = 0; i < lim && 
	     start_read_idx < all_reads_for_chrom.size(); ++i) {
	
	if (PROGRESS && i % PROGRESS_TIMING_MODULUS == 0) 
	  cerr << '\r' << chrom_name << ' ' << percent(i, chrom_cpgs) << "%\r";
	
	vector<epiread> current_window_reads;
	get_current_window_reads(all_reads_for_chrom, cpg_window, i, 
				 start_read_idx, current_window_reads);
	
	if (current_window_reads.size() > min_reads_per_window) {
	  const double score = (USE_BIC) ?
	    test_asm_bic(max_itr, low_prob, high_prob, current_window_reads) :
	    ((IGNORE_BALANCED_PARTITION_INFO) ?
	     test_asm_lrt(max_itr, low_prob, high_prob, current_window_reads) :
	     test_asm_lrt2(max_itr, low_prob, high_prob, current_window_reads));
	  if (score < critical_value || (USE_BIC && score < 0.0))
	    add_amr(chrom_name, i, cpg_window, current_window_reads, 
		    score, amrs);
	  ++windows_tested;
	}
      }
      if (PROGRESS)
	cerr << '\r' << chrom_name << " 100%" << endl;
    }
    
    const size_t windows_accepted = amrs.size();
    double fdr_cutoff = 0.0;
    
    if(amrs.empty())
	cerr << "No AMR is found. "<<endl;
    else{
    	if (!USE_BIC) {
      		fdr_cutoff = get_fdr_cutoff(total_cpgs, amrs, critical_value);
      		eliminate_amrs_by_fdr(fdr_cutoff, amrs);
    	}
    
    	if (VERBOSE)
	  cerr << "PROCESSED READS: " << total_reads_processed << endl
	       << "TOTAL CPGS: " << total_cpgs << endl
	       << "TESTED WINDOWS: " << windows_tested << endl
	       << "AMR WINDOWS: " << windows_accepted << endl
	       << "AMR/TESTED: " << (windows_accepted/
				     static_cast<double>(windows_tested)) << endl
	       << "FDR CUTOFF: " << fdr_cutoff << endl
	       << "AMR AFTER FDR: " << amrs.size() << endl;
	
	collapse_amrs(amrs);
	eio.convert_coordinates(amrs);
	merge_amrs(amrs,gap_limit);
	
	if (VERBOSE)
	  cerr << "MERGED AMRS: " << amrs.size() << endl;
	
	std::ofstream of;
	if (!outfile.empty()) of.open(outfile.c_str());
	std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
	copy(amrs.begin(), amrs.end(), 
	     std::ostream_iterator<GenomicRegion>(out, "\n"));
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
