/*    bsrate: a program for determining the rate of bisulfite
 *    conversion in a bisulfite sequencing experiment
 *
 *    Copyright (C) 2014-2017 University of Southern California and
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
#include <fstream>
#include <algorithm>
#include <numeric>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MappedRead.hpp"

#include "bsutils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::accumulate;

using std::unordered_map;


static void
revcomp(MappedRead &mr) {
  revcomp_inplace(mr.seq);
  if (mr.r.pos_strand())
    mr.r.set_strand('-');
  else mr.r.set_strand('+');
}


static void
count_states_pos(const bool INCLUDE_CPGS, const string &chrom,
                 const MappedRead &r,
                 vector<size_t> &unconv, vector<size_t> &conv,
                 vector<size_t> &err, size_t &hanging) {

  const size_t width = r.r.get_width();
  const size_t offset = r.r.get_start();

  size_t position = offset;
  assert(offset < chrom.length()); // at least one bp of read on chr
  for (size_t i = 0; i < width; ++i, ++position) {
    if (position >= chrom.length()) // some overhang
      ++hanging;

    if (is_cytosine(chrom[position]) &&
        (!is_guanine(chrom[position+1]) ||
         position == chrom.length() ||
         INCLUDE_CPGS)) {

      if (is_cytosine(r.seq[i])) ++unconv[i];
      else if (is_thymine(r.seq[i])) ++conv[i];
      else if (toupper(r.seq[i]) != 'N')
        ++err[i];
    }
  }
}


static void
count_states_neg(const bool INCLUDE_CPGS, const string &chrom,
                 const MappedRead &r,
                 vector<size_t> &unconv, vector<size_t> &conv,
                 vector<size_t> &err, size_t &hanging) {

  const size_t width = r.r.get_width();
  const size_t offset = r.r.get_start();

  size_t position = offset + width - 1;
  assert(position < chrom.length()); // at least one bp of read on chr
  for (size_t i = 0; i < width; ++i, --position) {
    if (position >= chrom.length())
      ++hanging;
    if (is_guanine(chrom[position]) &&
        (position == 0 ||
         !is_cytosine(chrom[position-1]) ||
         INCLUDE_CPGS)) {

      if (is_cytosine(r.seq[i])) ++unconv[i];
      else if (is_thymine(r.seq[i])) ++conv[i];
      else if (toupper(r.seq[i]) != 'N')
        ++err[i];
    }
  }
}


static void
write_output(const string &outfile,
             const vector<size_t> &ucvt_count_p,
             const vector<size_t> &cvt_count_p,
             const vector<size_t> &ucvt_count_n,
             const vector<size_t> &cvt_count_n,
             const vector<size_t> &err_p, const vector<size_t> &err_n) {

  // Get some totals first
  const size_t pos_cvt = accumulate(cvt_count_p.begin(),
                                    cvt_count_p.end(), 0ul);
  const size_t neg_cvt = accumulate(cvt_count_n.begin(),
                                    cvt_count_n.end(), 0ul);
  const size_t total_cvt = pos_cvt + neg_cvt;

  const size_t pos_ucvt =
    accumulate(ucvt_count_p.begin(), ucvt_count_p.end(), 0ul);
  const size_t neg_ucvt =
    accumulate(ucvt_count_n.begin(), ucvt_count_n.end(), 0ul);
  const size_t total_ucvt = pos_ucvt + neg_ucvt;

  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out << "OVERALL CONVERSION RATE = "
      << static_cast<double>(total_cvt)/(total_cvt + total_ucvt) << endl
      << "POS CONVERSION RATE = "
      << static_cast<double>(pos_cvt)/(pos_cvt + pos_ucvt) << '\t'
      << std::fixed << static_cast<size_t>(pos_cvt + pos_ucvt) << endl
      << "NEG CONVERSION RATE = "
      << static_cast<double>(neg_cvt)/(neg_cvt + neg_ucvt) << '\t'
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

  // Figure out how many positions to print in the output, capped at 1000
  size_t output_len =
    (ucvt_count_p.size() > 1000) ? 1000 : ucvt_count_p.size();

  while (output_len > 0 &&
         (ucvt_count_p[output_len-1] + cvt_count_p[output_len-1] +
          ucvt_count_n[output_len-1] + cvt_count_n[output_len-1] == 0))
    --output_len;

  // Now actually output the results
  static const size_t precision_val = 5;
  for (size_t i = 0; i < output_len; ++i) {

    const size_t total_p = ucvt_count_p[i] + cvt_count_p[i];
    const size_t total_n = ucvt_count_n[i] + cvt_count_n[i];
    const size_t total_valid = total_p + total_n;
    out << (i + 1) << "\t";

    out.precision(precision_val);
    out << total_p << '\t' << cvt_count_p[i] << '\t'
        << static_cast<double>(cvt_count_p[i])/max(size_t(1ul), total_p)
        << '\t';

    out.precision(precision_val);
    out << total_n << '\t' << cvt_count_n[i] << '\t'
        << static_cast<double>(cvt_count_n[i])/max(size_t(1ul), total_n)
        << '\t';

    const double total_cvt = cvt_count_p[i] + cvt_count_n[i];
    out.precision(precision_val);
    out << static_cast<size_t>(total_valid)
        << '\t' << cvt_count_p[i] + cvt_count_n[i] << '\t'
        << total_cvt/max(1ul, total_valid) << '\t';

    const double total_err = err_p[i] + err_n[i];
    out.precision(precision_val);
    const size_t total = total_valid + err_p[i] + err_n[i];
    out << err_p[i] + err_n[i] << '\t' << static_cast<size_t>(total) << '\t'
        << total_err/max(1ul, total) << endl;
  }
}


typedef unordered_map<string, string> chrom_file_map;
static void
get_chrom(const bool VERBOSE, const MappedRead &mr,
          const chrom_file_map& chrom_files,
          GenomicRegion &chrom_region, string &chrom) {

  const chrom_file_map::const_iterator fn(chrom_files.find(mr.r.get_chrom()));
  if (fn == chrom_files.end())
    throw SMITHLABException("could not find chrom: " + mr.r.get_chrom());

  chrom.clear();
  read_fasta_file(fn->second, mr.r.get_chrom(), chrom);
  if (chrom.empty())
    throw SMITHLABException("could not find chrom: " + mr.r.get_chrom());

  chrom_region.set_chrom(mr.r.get_chrom());
}


int
main(int argc, const char **argv) {

  try {

    // ASSUMED MAXIMUM LENGTH OF A FRAGMENT
    static const size_t OUTPUT_SIZE = 10000;

    bool VERBOSE = false;
    bool INCLUDE_CPGS = false;
    bool A_RICH_READS = false;

    string chrom_file;
    string outfile;
    string fasta_suffix = "fa";

    double max_mismatches = std::numeric_limits<double>::max();

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "Program to compute the "
                           "BS conversion rate from BS-seq "
                           "reads mapped to a genome",
                           "-c <chroms> <mapped-reads>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("chrom", 'c', "file or dir of chroms (FASTA format; "
                      ".fa suffix)", true , chrom_file);
    //!!!!!! OPTION IS HIDDEN BECAUSE USERS DON'T NEED TO CHANGE IT...
    //     opt_parse.add_opt("suffix", 's', "suffix of FASTA files "
    //                "(assumes -c indicates dir)",
    //                false , fasta_suffix);
    opt_parse.add_opt("all", 'N', "count all Cs (including CpGs)",
                      false , INCLUDE_CPGS);
    opt_parse.add_opt("max", 'M', "max mismatches (can be fractional)",
                      false , max_mismatches);
    opt_parse.add_opt("a-rich", 'A', "reads are A-rich", false, A_RICH_READS);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE && max_mismatches != std::numeric_limits<double>::max())
      cerr << "MAX_MISMATCHES=" << max_mismatches << endl;

    chrom_file_map chrom_files;
    identify_and_read_chromosomes(chrom_file, fasta_suffix, chrom_files);
    if (VERBOSE)
      cerr << "N_CHROMS=" << chrom_files.size() << endl;

    std::ifstream in(mapped_reads_file.c_str());
    if (!in)
      throw SMITHLABException("cannot open file: " + mapped_reads_file);

    vector<size_t> unconv_count_pos(OUTPUT_SIZE, 0ul);
    vector<size_t> conv_count_pos(OUTPUT_SIZE, 0ul);
    vector<size_t> unconv_count_neg(OUTPUT_SIZE, 0ul);
    vector<size_t> conv_count_neg(OUTPUT_SIZE, 0ul);
    vector<size_t> err_pos(OUTPUT_SIZE, 0ul);
    vector<size_t> err_neg(OUTPUT_SIZE, 0ul);

    string chrom;
    MappedRead mr;
    GenomicRegion chrom_region; // exists only for faster comparison
    size_t hanging = 0;

    while (in >> mr) {

      if (A_RICH_READS)
        revcomp(mr);

      // get the correct chrom if it has changed
      if (chrom.empty() || !mr.r.same_chrom(chrom_region))
        get_chrom(VERBOSE, mr, chrom_files, chrom_region, chrom);

      // do the work for this mapped read
      if (mr.r.pos_strand())
        count_states_pos(INCLUDE_CPGS, chrom, mr,
                         unconv_count_pos, conv_count_pos, err_pos, hanging);
      else
        count_states_neg(INCLUDE_CPGS, chrom, mr,
                         unconv_count_neg, conv_count_neg, err_neg, hanging);
    }
    write_output(outfile, unconv_count_pos,
                 conv_count_pos, unconv_count_neg,
                 conv_count_neg, err_pos, err_neg);
    if (hanging > 0) // some overhanging reads
      cerr << "Warning: a nonzero number (" << hanging << ") of reads mapped"
           << " to the very end of a chromosome. For high numbers, make"
           << " sure you are using the same assembly you mapped with."
           << endl;
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
