/*    bsrate: a program for determining the rate of bisulfite
 *    conversion in a bisulfite sequencing experiment
 *
 *    Copyright (C) 2009-2012 University of Southern California and
 *                            Andrew D. Smith
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

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MappedRead.hpp"

#include "bsutils.hpp"
#include "FileIterator.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::accumulate;

using std::tr1::unordered_map;

static void
count_states_pos(const bool COUNT_CPGS, const string &chrom, const MappedRead &r,
		 vector<size_t> &unconv, vector<size_t> &conv, 
		 vector<size_t> &err) {
  
  const size_t width = r.r.get_width();
  const size_t offset = r.r.get_start();
  
  size_t position = offset;
  for (size_t i = 0; i < width; ++i, ++position) {
    assert(position < chrom.length());
    if (is_cytosine(chrom[position]) && 
	(!is_guanine(chrom[position + 1]) || 
	 position == chrom.length() || COUNT_CPGS)) {
      if (is_cytosine(r.seq[i])) ++unconv[i];
      else if (is_thymine(r.seq[i])) ++conv[i];
      else if (toupper(r.seq[i]) != 'N')
	++err[i];
    }
  }
}


static void
count_states_neg(const bool COUNT_CPGS, const string &chrom, const MappedRead &r,
		 vector<size_t> &unconv, vector<size_t> &conv, 
		 vector<size_t> &err) {
  
  const size_t width = r.r.get_width();
  const size_t offset = r.r.get_start();
  
  size_t position = offset + width - 1;
  for (size_t i = 0; i < width; ++i, --position) {
    assert(position < chrom.length());
    if (is_guanine(chrom[position]) && 
	(!is_cytosine(chrom[position - 1]) || position == 0 || COUNT_CPGS)) {
      if (is_cytosine(r.seq[i])) ++unconv[i];
      else if (is_thymine(r.seq[i])) ++conv[i];
      else if (toupper(r.seq[i]) != 'N')
	++err[i];
    }
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
write_output(const string &outfile,
	     const vector<size_t> &ucvt_count_p, const vector<size_t> &cvt_count_p,
	     const vector<size_t> &ucvt_count_n, const vector<size_t> &cvt_count_n,
	     const vector<size_t> &err_p, const vector<size_t> &err_n) {
  
  // Get some totals first
  const double pos_cvt = accumulate(cvt_count_p.begin(), cvt_count_p.end(), 0);
  const double neg_cvt = accumulate(cvt_count_n.begin(), cvt_count_n.end(), 0); 
  const double total_cvt = pos_cvt + neg_cvt;
  const double pos_ucvt = accumulate(ucvt_count_p.begin(), ucvt_count_p.end(), 0);
  const double neg_ucvt = accumulate(ucvt_count_n.begin(), ucvt_count_n.end(), 0);
  const double total_ucvt = pos_ucvt + neg_ucvt;
  
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  
  out << "OVERALL CONVERSION RATE = "
      << total_cvt/(total_cvt + total_ucvt) << endl
      << "POS CONVERSION RATE = "
      << pos_cvt/(pos_cvt + pos_ucvt) << '\t'
      << std::fixed << static_cast<size_t>(pos_cvt + pos_ucvt) << endl
      << "NEG CONVERSION RATE = "
      << neg_cvt/(neg_cvt + neg_ucvt) << '\t'
      << std::fixed << static_cast<size_t>(neg_cvt + neg_ucvt) << endl;
  
  out << "BASE" << '\t'
      << "PTOT" << '\t'
      << "PCONV" << '\t'
      << "PRATE" << '\t'
      << "NTOT" << '\t'
      << "NCONV" << '\t'
      << "NRATE" << '\t'
      << "BTHTOT" << '\t'
      << "BTHCONV" << '\t'
      << "BTHRATE" << '\t'
      << "ERR" << '\t'
      << "ALL" << '\t'
      << "ERRRATE"  << endl;

  // Figure out how many positions to print in the output
  size_t output_len = ucvt_count_p.size();
  while (output_len > 0 && 
	 (ucvt_count_p[output_len-1] + cvt_count_p[output_len-1] +
	  ucvt_count_n[output_len-1] + cvt_count_n[output_len-1] > 0))
    --output_len;
  
  // Now actually output the results
  static const size_t precision_val = 5;
  for (size_t i = 0; i < output_len; ++i) {
    const double total_p = ucvt_count_p[i] + cvt_count_p[i];
    const double total_n = ucvt_count_n[i] + cvt_count_n[i];
    const double total_valid = total_p + total_n;
    out << (i + 1) << "\t";
    
    out.precision(precision_val);
    out << static_cast<size_t>(total_p) << '\t'
	<< cvt_count_p[i] << '\t' << cvt_count_p[i]/max(1.0, total_p) << '\t';
    
    out.precision(precision_val);
    out << static_cast<size_t>(total_n) << '\t' << cvt_count_n[i] << '\t' 
	<< cvt_count_n[i]/max(1.0, total_n) << '\t';
    
    out.precision(precision_val);
    out << static_cast<size_t>(total_valid) 
	<< '\t' << cvt_count_p[i] + cvt_count_n[i] << '\t'
	<< (cvt_count_p[i] + cvt_count_n[i])/max(1.0, total_valid) << '\t';
    
    out.precision(precision_val);
    const double total = total_valid + err_p[i] + err_n[i];
    out << err_p[i] + err_n[i] << '\t' << static_cast<size_t>(total) << '\t'
	<< (err_p[i] + err_n[i])/max(1.0, total) << endl;
  }
}


typedef unordered_map<string, string> chrom_file_map;
static void
get_chrom(const bool VERBOSE, const MappedRead &mr, 
	  const chrom_file_map& chrom_files, GenomicRegion &chrom_region, 
	  string &chrom) {
  const chrom_file_map::const_iterator fn(chrom_files.find(mr.r.get_chrom()));
  if (fn == chrom_files.end())
    throw SMITHLABException("could not find chrom: " + mr.r.get_chrom());
  chrom.clear();
  vector<string> chrom_names, chroms;
  read_fasta_file(fn->second, chrom_names, chroms);
  if (VERBOSE) {
    if (chrom_names.size() > 1)
      cerr << "WARNING: multiple sequences in " << fn->second << endl;
    cerr << "PROCESSING: " << mr.r.get_chrom() << endl;
  }
  chrom_region.set_chrom(chrom_names.front());
  chrom.swap(chroms.front());
}


int 
main(int argc, const char **argv) {
  
  try {
    
    // ASSUMED MAXIMUM LENGTH OF A FRAGMENT
    static const size_t OUTPUT_SIZE = 10000;
    
    bool VERBOSE = false;
    bool COUNT_CPGS = false;
    
    string chrom_file;
    string outfile;
    string fasta_suffix = "fa";
    
    double max_mismatches = std::numeric_limits<double>::max();
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "Program to compute the "
			   "bs conversion rate from BS-seq "
			   "reads mapped to a genome",
			   "-c <chroms> <mapped-reads>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("chrom", 'c', "file or dir of chroms (FASTA format; .fa suffix)",
		      true , chrom_file);
    //!!!!!! OPTION IS HIDDEN BECAUSE USERS DON'T NEED TO CHANGE IT...
    //     opt_parse.add_opt("suffix", 's', "suffix of FASTA files "
    // 		      "(assumes -c indicates dir)", 
    // 		      false , fasta_suffix);
    opt_parse.add_opt("all", 'N', "count all Cs (including CpGs)", 
		      false , COUNT_CPGS);
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

    if (VERBOSE && max_mismatches != std::numeric_limits<double>::max())
      cerr << "MAX MISMATCHES=" << max_mismatches << endl;
    
    chrom_file_map chrom_files;
    identify_chromosomes(chrom_file, fasta_suffix, chrom_files);
    if (VERBOSE)
      cerr << "N_CHROMS:\t" << chrom_files.size() << endl;
    
    std::ifstream in(mapped_reads_file.c_str());
    if (!in) 
      throw SMITHLABException("cannot open input file " + mapped_reads_file);

    vector<size_t> unconv_count_pos(OUTPUT_SIZE, 0ul);
    vector<size_t> conv_count_pos(OUTPUT_SIZE, 0ul);
    vector<size_t> unconv_count_neg(OUTPUT_SIZE, 0ul);
    vector<size_t> conv_count_neg(OUTPUT_SIZE, 0ul);
    vector<size_t> err_pos(OUTPUT_SIZE, 0ul);
    vector<size_t> err_neg(OUTPUT_SIZE, 0ul);
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    
    string chrom;
    MappedRead mr;
    GenomicRegion chrom_region;
    while (!in.eof() && in >> mr) {
      // get the correct chrom if it has changed
      if (chrom.empty() || !mr.r.same_chrom(chrom_region))
	get_chrom(VERBOSE, mr, chrom_files, chrom_region, chrom);
      
      // do the work for this mapped read
      if (mr.r.pos_strand())
	count_states_pos(COUNT_CPGS, chrom, mr, 
			 unconv_count_pos, conv_count_pos, err_pos);
      else 
	count_states_neg(COUNT_CPGS, chrom, mr, 
			 unconv_count_neg, conv_count_neg, err_neg);
    }
    write_output(outfile, unconv_count_pos, 
		 conv_count_pos, unconv_count_neg,
		 conv_count_neg, err_pos, err_neg);
    
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
