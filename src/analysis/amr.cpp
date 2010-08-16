/*    methcounts: a program for counting the methylated and
 *    unmethylated reads mapping over each CpG or C
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
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>

#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"

#include "bsutils.hpp"
#include "FileIterator.hpp"
#include "MappedRead.hpp"


#include <gsl/gsl_sf_gamma.h>

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;


static inline double
log_sum_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}

static double
log_hyper_g(const size_t a, const size_t c, 
	    const size_t b, const size_t d) {
  return  gsl_sf_lnfact(a + b) + gsl_sf_lnfact(c + d) + 
    gsl_sf_lnfact(a + c) + gsl_sf_lnfact(b + d) -
    gsl_sf_lnfact(a + b + c + d) - gsl_sf_lnfact(a) - 
    gsl_sf_lnfact(b) - gsl_sf_lnfact(c) - gsl_sf_lnfact(d);
}

static double
test_similar_population(size_t meth_a, size_t unmeth_a, 
			size_t meth_b, size_t unmeth_b) {
  double p = 0;
  while (unmeth_a > 0 && meth_b > 0) {
    p = log_sum_log(p, log_hyper_g(meth_a, meth_b, unmeth_a, unmeth_b));
    ++meth_a; --unmeth_a;
    ++unmeth_b; --meth_b;
  }
  p = log_sum_log(p, log_hyper_g(meth_a, meth_b, unmeth_a, unmeth_b));
  return exp(p);
}


static void
add_contribution_cpg(const size_t offset_a, 
		     const size_t offset_b, 
		     const MappedRead &r,
		     size_t &meth_meth, size_t &meth_unmeth,
		     size_t &unmeth_meth, size_t &unmeth_unmeth) {
  if (r.r.pos_strand() && 
      (r.r.get_start() <= offset_a) &&
      (r.r.get_start() <= offset_b)) {
    const size_t position_a = offset_a - r.r.get_start();
    if (!(position_a < r.seq.length())) return;
    assert(position_a < r.seq.length());
    const size_t position_b = offset_b - r.r.get_start();
    if (!(position_b < r.seq.length())) return;
    assert(position_b < r.seq.length());
    if (is_cytosine(r.seq[position_a])) {
      if (is_cytosine(r.seq[position_b])) ++meth_meth;
      else if (is_thymine(r.seq[position_b])) ++meth_unmeth;
    }
    else if (is_thymine(r.seq[position_a])) {
      if (is_thymine(r.seq[position_b])) ++unmeth_unmeth;
      else if (is_cytosine(r.seq[position_b])) ++unmeth_meth;
    }
    
  }
  if (r.r.neg_strand() && 
      (offset_a - r.r.get_start() + 2 <= r.seq.length()) &&
      (offset_b - r.r.get_start() + 2 <= r.seq.length())) {
    const size_t position_a = (r.seq.length() - 1) - 
      ((offset_a + 1) - r.r.get_start());
    if (!(position_a < r.seq.length())) return;
    assert(position_a < r.seq.length());
    const size_t position_b = (r.seq.length() - 1) - 
      ((offset_b + 1) - r.r.get_start());
    if (!(position_b < r.seq.length())) return;
    assert(position_b < r.seq.length());
    if (is_cytosine(r.seq[position_a])) {
      if (is_cytosine(r.seq[position_b])) ++meth_meth;
      else if (is_thymine(r.seq[position_b])) ++meth_unmeth;
    }
    else if (is_thymine(r.seq[position_a])) {
      if (is_thymine(r.seq[position_b])) ++unmeth_unmeth;
      else if (is_cytosine(r.seq[position_b])) ++unmeth_meth;
    }
  }
}


static void
add_contribution_c(const size_t offset, const MappedRead &r,
		   size_t &meth, size_t &unmeth) {
  if (r.r.pos_strand()) {
    const size_t position = offset - r.r.get_start();
    assert(position < r.seq.length());
    if (is_cytosine(r.seq[position])) ++meth;
    if (is_thymine(r.seq[position])) ++unmeth;
  }
}


static void
add_contribution_g(const size_t offset, const MappedRead &r,
		   size_t &meth, size_t &unmeth) {
  if (r.r.neg_strand()) {
    const size_t position = r.seq.length() - (offset - r.r.get_start() + 1);
    assert(position < r.seq.length());
    if (is_cytosine(r.seq[position])) ++meth;
    if (is_thymine(r.seq[position])) ++unmeth;
  }
}


static bool
precedes(const MappedRead &r, const size_t offset) {
  return r.r.get_end() <= offset;
}


static bool
succeeds(const MappedRead &r, const size_t offset) {
  return r.r.get_start() > offset;
}


static void
advance(const size_t first, const size_t last,
	const GenomicRegion &chrom_region, 
	FileIterator<MappedRead> &regions) {
  while (regions.last_is_good() && 
	 chrom_region.same_chrom(regions.get_last()->r) &&
	 !succeeds(*regions.get_last(), last)) {
    regions.increment_last();
  }
  //   if (regions.last_is_good() != reads.last_is_good())
  //     throw RMAPException("read and map files seem out of sync");
  while (regions.first_is_good() && 
	 chrom_region.same_chrom(regions.get_first()->r) &&
	 precedes(*regions.get_first(), first)) {
    regions.increment_first();
  }
  //   if (regions.first_is_good() != reads.first_is_good())
  //     throw RMAPException("read and map files seem out of sync");
}


static void
scan_chromosome_cpg(const string &chrom,
		    const GenomicRegion &chrom_region,
		    const double max_mismatches,
		    FileIterator<MappedRead> &regions,
		    std::ostream &out) {
  const string chrom_name(chrom_region.get_chrom());
  size_t prev = std::numeric_limits<size_t>::max();
  size_t read_width = 101;
  for (size_t i = 0; i < chrom.length() - 1 && regions.first_is_good(); ++i) {
    if (is_cpg(chrom, i)) {
      if (prev != std::numeric_limits<size_t>::max()) {
	if (i - prev <= read_width) {
	  advance(prev, i + 1, chrom_region, regions);
	  size_t meth_meth_count = 1, meth_unmeth_count = 1;
	  size_t unmeth_meth_count = 1, unmeth_unmeth_count = 1;
	  for (vector<MappedRead>::const_iterator j(regions.get_first());
	       j != regions.get_last(); ++j)
	    if (j->r.get_score() <= max_mismatches)
	      add_contribution_cpg(prev, i, *j, 
				   meth_meth_count, meth_unmeth_count,
				   unmeth_meth_count, unmeth_unmeth_count);
	  const double total = meth_meth_count + meth_unmeth_count +
	    unmeth_meth_count + unmeth_unmeth_count;
	  out << chrom_name << "\t" << i << "\t" << i + 1 << "\tCpG:" 
	      << total << "\t" 
	      << test_similar_population(meth_meth_count, meth_unmeth_count,
					 unmeth_meth_count, unmeth_unmeth_count) << "\t"
	    //  (meth_meth_count + unmeth_unmeth_count)/max(1.0, total) << "\t+\t"
	      << meth_meth_count << "\t"
	      << meth_unmeth_count << "\t"
	      << unmeth_meth_count << "\t"
	      << unmeth_unmeth_count << "\n";
	}
      }
      prev = i;
    }
  }
}


static void
scan_chromosome(const string &chrom, const GenomicRegion &chrom_region,
		const double max_mismatches,
		FileIterator<MappedRead> &regions, 
		std::ostream &out) {
  const string chrom_name(chrom_region.get_chrom());
  
  for (size_t i = 1; i < chrom.length() - 1 && regions.first_is_good(); ++i) {
    advance(i, i, chrom_region, regions);
    if (is_cytosine(chrom[i]) && !is_guanine(chrom[i + 1])) {
      size_t meth_count = 0, unmeth_count = 0;
      for (vector<MappedRead>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j)
	if (j->r.get_score() <= max_mismatches)
	  add_contribution_c(i, *j, meth_count, unmeth_count);
      const double total = meth_count + unmeth_count;
      out << chrom_name << "\t" << i << "\t" << i + 1 << "\tC:"
	  << total << "\t" << meth_count/max(1.0, total) << "\t+\n";
    }
    if (is_guanine(chrom[i]) && !is_cytosine(chrom[i - 1])) {
      size_t meth_count = 0, unmeth_count = 0;
      for (vector<MappedRead>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j)
	if (j->r.get_score() <= max_mismatches)
	  add_contribution_g(i, *j, meth_count, unmeth_count);
      const double total = meth_count + unmeth_count;
      out << chrom_name << "\t" << i << "\t" << i + 1 << "\tG:" 
	  << total << "\t" << meth_count/max(1.0, total) << "\t+\n";
    }
  }
}


static void
identify_chromosomes(const bool VERBOSE, const string chrom_file,
		     const string fasta_suffix, vector<string> &chrom_files) {
  if (VERBOSE)
    cerr << "[IDENTIFYING CHROMS] ";
  if (isdir(chrom_file.c_str())) 
    read_dir(chrom_file, fasta_suffix, chrom_files);
  else chrom_files.push_back(chrom_file);
  if (VERBOSE) {
    cerr << "[DONE]" << endl 
	 << "chromosome files found (approx size):" << endl;
    for (vector<string>::const_iterator i = chrom_files.begin();
	 i != chrom_files.end(); ++i)
      cerr << *i << " (" << roundf(get_filesize(*i)/1e06) << "Mbp)" << endl;
    cerr << endl;
  }
}


static void
advance_chromosome(const GenomicRegion &chrom_region, 
		   FileIterator<MappedRead> &regions) {
  while (regions.last_is_good() && 
	 (regions.get_last()->r < chrom_region)) {
    assert(regions.last_is_good());
    regions.increment_last();
  }
  while (regions.first_is_good() && 
	 (regions.get_first()->r < chrom_region)) {
    regions.increment_first();
  }
}


static void
scan_chroms(const bool VERBOSE, const bool PROCESS_NON_CPGS,
	    const double max_mismatches,
	    const string &outfile, const vector<string> &chrom_files, 
	    FileIterator<MappedRead> &regions) {
  std::ostream *out = (outfile.empty()) ? &cout : new std::ofstream(outfile.c_str());
  for (size_t i = 0; i < chrom_files.size(); ++i) {
    const string fn(strip_path_and_suffix(chrom_files[i]));
    if (VERBOSE)
      cerr << "[LOADING CHROM FILE=" << fn << "]";
    vector<string> chrom_names, chroms;
    read_fasta_file(chrom_files[i].c_str(), chrom_names, chroms);
    for (size_t j = 0; j < chroms.size(); ++j) {
      if (VERBOSE) cerr << "[SCANNING=" << chrom_names[j] << "]";
      //TODO: WHAT HAPPENS IF A CHROM IS MISSING??
      const GenomicRegion chrom_region(chrom_names[j], 0, 0);
      advance_chromosome(chrom_region, regions);
      if (PROCESS_NON_CPGS)
	scan_chromosome(chroms[j], chrom_region, max_mismatches,
			regions, *out);
      else scan_chromosome_cpg(chroms[j], chrom_region, max_mismatches,
			       regions, *out);
    }
    if (VERBOSE) cerr << " [DONE]" << endl;
  }
  if (out != &cout) delete out;
}


int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    bool PROCESS_NON_CPGS = false;
    
    string chrom_file;
    string outfile;
    string fasta_suffix = "fa";
    
    size_t BUFFER_SIZE = 100000;
    double max_mismatches = std::numeric_limits<double>::max();
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], "a program for counting the "
			   "methylated and unmethylated reads mapping "
			   "over each CpG or C.",
			   "<fasta-reads>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("chrom", 'c', "FASTA file or dir containing chromosome(s)", 
		      true , chrom_file);
    opt_parse.add_opt("suffix", 's', "suffix of FASTA files "
		      "(assumes -c indicates dir)", 
		      false , fasta_suffix);
    opt_parse.add_opt("non", 'N', "process non-CpG cytosines", 
		      false , PROCESS_NON_CPGS);
    opt_parse.add_opt("buffer", 'B', "buffer size (in records, not bytes)", 
		      false , BUFFER_SIZE);
    opt_parse.add_opt("max", 'M', "max mismatches (can be fractional)", 
		      false , max_mismatches);
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
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "MAX MISMATCHES=" << max_mismatches << endl;

    //     const bool FASTQ = is_fastq(reads_file);
    //     if (VERBOSE)
    //       cerr << "READS FILE FORMAT: " << ((FASTQ) ? "FASTQ" : "FASTA") << endl;
    
    vector<string> chrom_files;
    identify_chromosomes(VERBOSE, chrom_file, fasta_suffix, chrom_files);
    sort(chrom_files.begin(), chrom_files.end());
    
    FileIterator<MappedRead> regions(mapped_reads_file, BUFFER_SIZE);
    scan_chroms(VERBOSE, PROCESS_NON_CPGS, max_mismatches, 
		outfile, chrom_files, regions);
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
