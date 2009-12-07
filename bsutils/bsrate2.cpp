/*    bsrate2: a program for determining the rate of bisulfite
 *    conversion in a bisulfite sequencing experiment, with too many
 *    reads to load into memory at once.
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

// TODO: check that the mapped locations and the read sequences are
// equal in number.


#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"

#include "BSUtils.hpp"
#include "FileIterator.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

static void
add_contribution_c(const size_t offset, const GenomicRegion &r,
		   const string &s, size_t &unconv, size_t &conv) {
  if (r.pos_strand()) {
    const size_t position = offset - r.get_start();
    assert(position < s.length());
    if (is_cytosine(s[position])) ++unconv;
    if (is_thymine(s[position])) ++conv;
  }
}


static void
add_contribution_g(const size_t offset, const GenomicRegion &r,
		   const string &s, size_t &unconv, size_t &conv) {
  if (r.neg_strand()) {
    const size_t position = s.length() - (offset - r.get_start() + 1);
    assert(position < s.length());
    if (is_cytosine(s[position])) ++unconv;
    if (is_thymine(s[position])) ++conv;
  }
}


static void
add_contribution_c(const size_t offset, const GenomicRegion &r,
		   const FASTQRecord &s, size_t &unconv, size_t &conv) {
  if (r.pos_strand()) {
    const size_t position = offset - r.get_start();
    // if (position < 40) return;
    assert(position < s.first.length());
    //     if (s.second[position] > '`') {
    if (is_cytosine(s.first[position])) ++unconv;
    if (is_thymine(s.first[position])) ++conv;
  }
}


static void
add_contribution_g(const size_t offset, const GenomicRegion &r,
		   const FASTQRecord &s, size_t &unconv, size_t &conv) {
  if (r.neg_strand()) {
    const size_t position = s.first.length() - (offset - r.get_start() + 1);
    // if (position < 40) return;
    assert(position < s.first.length());
    //     if (s.second[position] > '`') {
    if (is_cytosine(s.first[position])) ++unconv;
    if (is_thymine(s.first[position])) ++conv;
    //     }
  }
}


static bool
precedes(const GenomicRegion &r, const size_t offset) {
  return r.get_end() <= offset;
}


static bool
succeeds(const GenomicRegion &r, const size_t offset) {
  return r.get_start() > offset;
}


template <class T> void
advance(const size_t first, const size_t last,
	const GenomicRegion &chrom_region, 
	FileIterator<GenomicRegion> &regions,
	FileIterator<T> &reads) {
  while (regions.last_is_good() && 
	 chrom_region.same_chrom(*regions.get_last()) &&
	 !succeeds(*regions.get_last(), last)) {
    regions.increment_last();
    reads.increment_last();
  }
  //   if (regions.last_is_good() != reads.last_is_good())
  //     throw RMAPException("read and map files seem out of sync");
  while (regions.first_is_good() && 
	 chrom_region.same_chrom(*regions.get_first()) &&
	 precedes(*regions.get_first(), first)) {
    regions.increment_first();
    reads.increment_first();
  }
  //   if (regions.first_is_good() != reads.first_is_good())
  //     throw RMAPException("read and map files seem out of sync");
}


template <class T> void
scan_chromosome(const string &chrom, const GenomicRegion &chrom_region,
		const double max_mismatches,
		FileIterator<GenomicRegion> &regions, 
		FileIterator<T> &reads,
 		size_t &unconv_count_pos, size_t &conv_count_pos,
 		size_t &unconv_count_neg, size_t &conv_count_neg) {
  const string chrom_name(chrom_region.get_chrom());
  
  
  for (size_t i = 0; i < chrom.length() - 1 && regions.first_is_good(); ++i) {
    advance(i, i, chrom_region, regions, reads);
    if (is_cytosine(chrom[i]) && !is_guanine(chrom[i + 1])) {
      typename vector<T>::const_iterator k(reads.get_first());
      for (vector<GenomicRegion>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j, ++k) {
	if (j->get_score() <= max_mismatches)
	  add_contribution_c(i, *j, *k, unconv_count_pos, conv_count_pos);
      }
    }
    if (is_guanine(chrom[i]) && !is_cytosine(chrom[i - 1])) {
      typename vector<T>::const_iterator k(reads.get_first());
      for (vector<GenomicRegion>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j, ++k)
	if (j->get_score() <= max_mismatches)
	  add_contribution_g(i, *j, *k, unconv_count_neg, conv_count_neg);
    }
  }
}


template <class T> void
scan_chromosome_cpg(const string &chrom, const GenomicRegion &chrom_region,
		    const double max_mismatches,
		    FileIterator<GenomicRegion> &regions, 
		    FileIterator<T> &reads,
		    size_t &unconv_count_pos, size_t &conv_count_pos,
		    size_t &unconv_count_neg, size_t &conv_count_neg) {
  const string chrom_name(chrom_region.get_chrom());
  
  for (size_t i = 0; i < chrom.length() - 1 && regions.first_is_good(); ++i) {
    advance(i, i, chrom_region, regions, reads);
    if (is_cytosine(chrom[i])) {
      typename vector<T>::const_iterator k(reads.get_first());
      for (vector<GenomicRegion>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j, ++k) {
	if (j->get_score() <= max_mismatches)
	  add_contribution_c(i, *j, *k, unconv_count_pos, conv_count_pos);
      }
    }
    if (is_guanine(chrom[i])) {
      typename vector<T>::const_iterator k(reads.get_first());
      for (vector<GenomicRegion>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j, ++k)
	if (j->get_score() <= max_mismatches)
	  add_contribution_g(i, *j, *k, unconv_count_neg, conv_count_neg);
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


template <class T> void
advance_chromosome(const GenomicRegion &chrom_region, 
		   FileIterator<GenomicRegion> &regions, 
		   FileIterator<T> &reads) {
  while (regions.last_is_good() && 
	 (*regions.get_last() < chrom_region)) {
    assert(regions.last_is_good());
    regions.increment_last();
    reads.increment_last();
  }
  while (regions.first_is_good() && 
	 (*regions.get_first() < chrom_region)) {
    regions.increment_first();
    reads.increment_first();
  }
}


template <class T> 
void
scan_chroms(const bool VERBOSE, const bool PROCESS_CPGS,
	    const double max_mismatches,
	    const string &outfile, const vector<string> &chrom_files, 
	    FileIterator<GenomicRegion> &regions,
	    FileIterator<T> &reads) {

  size_t unconv_count_pos = 0, conv_count_pos = 0;
  size_t unconv_count_neg = 0, conv_count_neg = 0;
  
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
      advance_chromosome(chrom_region, regions, reads);
      if (PROCESS_CPGS)
	scan_chromosome_cpg(chroms[j], chrom_region, max_mismatches,
			    regions, reads, 
			    unconv_count_pos, conv_count_pos,
			    unconv_count_neg, conv_count_neg);
      else scan_chromosome(chroms[j], chrom_region, max_mismatches, 
			   regions, reads, 
			   unconv_count_pos, conv_count_pos,
			   unconv_count_neg, conv_count_neg);
    }
    if (VERBOSE) cerr << " [DONE]" << endl;
  }

  std::ostream *out = (outfile.empty()) ? &cout : new std::ofstream(outfile.c_str());
  const double total_pos = unconv_count_pos + conv_count_pos;
  const double total_neg = unconv_count_neg + conv_count_neg;
  *out << "\tTOTAL\tCONVERTED\tRATE" << endl
       << "POS:\t" << static_cast<size_t>(total_pos) << "\t" 
       << conv_count_pos << "\t" 
       << conv_count_pos/total_pos << endl
       << "NEG:\t" << static_cast<size_t>(total_neg) << "\t" 
       << conv_count_neg << "\t" 
       << conv_count_neg/total_neg << endl
       << "BOTH:\t" << total_pos + total_neg << "\t" 
       << conv_count_pos + conv_count_neg << "\t" 
       << (conv_count_pos + conv_count_neg)/(total_pos + total_neg) << endl;
  if (out != &cout) delete out;
}


int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    bool PROCESS_CPGS = false;
    
    string mapped_file;
    string chrom_file;
    string outfile;
    string fasta_suffix = "fa";
    
    size_t BUFFER_SIZE = 100000;
    double max_mismatches = std::numeric_limits<double>::max();
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("bsrate2", "a program for determining the "
			   "rate of bisulfite conversion in a "
			   "bisulfite sequencing experiment, with "
			   "too many reads to load into memory at once.",
			   "<FASTA_READS>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("chrom", 'c', "FASTA file or dir containing chromosome(s)", 
		      true , chrom_file);
    opt_parse.add_opt("suffix", 's', "suffix of FASTA files "
		      "(assumes -c indicates dir)", 
		      false , fasta_suffix);
    opt_parse.add_opt("mapped", 'm', "file of mapped locations", 
		      true , mapped_file);
    opt_parse.add_opt("all", 'A', "process all Cs", 
		      false , PROCESS_CPGS);
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
    const string reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "MAX MISMATCHES=" << max_mismatches << endl;

    const bool FASTQ = is_fastq(reads_file);
    if (VERBOSE)
      cerr << "READS FILE FORMAT: " << ((FASTQ) ? "FASTQ" : "FASTA") << endl;
    
    vector<string> chrom_files;
    identify_chromosomes(VERBOSE, chrom_file, fasta_suffix, chrom_files);
    sort(chrom_files.begin(), chrom_files.end());
    
    FileIterator<GenomicRegion> regions(mapped_file, BUFFER_SIZE);
    if (FASTQ) {
      FileIterator<FASTQRecord> reads(reads_file, BUFFER_SIZE);
      scan_chroms(VERBOSE, PROCESS_CPGS, max_mismatches,
		  outfile, chrom_files, regions, reads);
    }
    else {
      FileIterator<string> reads(reads_file, BUFFER_SIZE);
      scan_chroms(VERBOSE, PROCESS_CPGS, max_mismatches, 
		  outfile, chrom_files, regions, reads);
    }
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
