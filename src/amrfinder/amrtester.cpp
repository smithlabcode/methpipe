/*    amrtester: A program for testing whether a genomic region has
 *    allele-specific methylation
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Benjamin E Decato and Andrew D. Smith and Fang Fang
 *
 *    Authors: Andrew D. Smith and Benjamin E Decato and Fang Fang
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

#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <smithlab_os.hpp>
#include <GenomicRegion.hpp>

#include "Epiread.hpp"
#include "EpireadStats.hpp"

using std::streampos;
using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::unordered_map;


static void
backup_to_start_of_current_record(std::ifstream &in) {
  while (in.tellg() > 0 && in.peek() != '\n' && in.peek() != '\r') {
    in.seekg(-1, std::ios_base::cur);
  }
}


static streampos
find_first_epiread_ending_after_position(const string &query_chrom,
                                         const size_t query_pos,
                                         std::ifstream &in) {
  in.seekg(0, std::ios_base::end);
  size_t high_pos = in.tellg();
  size_t eof = in.tellg();
  in.seekg(0, std::ios_base::beg);
  size_t low_pos = 0;

  string chrom, seq;
  size_t start = 0ul;

  // This is just binary search on disk
  while (high_pos > low_pos + 1) {
    const size_t mid_pos = (low_pos + high_pos)/2;

    in.seekg(mid_pos);
    backup_to_start_of_current_record(in);

    // we've hit the end of file without finding an epiread
    if(low_pos == eof-2)
      return -1;

    if (!(in >> chrom >> start >> seq)) {
      throw SMITHLABException("problem loading reads");
    }
    if (chrom < query_chrom ||
        (chrom == query_chrom && start + seq.length() <= query_pos))
      low_pos = mid_pos;
    else
      high_pos = mid_pos;
  }
  return low_pos;
}


static void
load_reads(const string &reads_file_name,
           const GenomicRegion &region, vector<epiread> &the_reads) {

  // open and check the file
  std::ifstream in(reads_file_name.c_str());
  if (!in)
    throw SMITHLABException("cannot open input file " + reads_file_name);

  const string query_chrom(region.get_chrom());
  const size_t query_start = region.get_start();
  const size_t query_end = region.get_end();
  const streampos low_offset =
    find_first_epiread_ending_after_position(query_chrom, query_start, in);

  in.seekg(low_offset, std::ios_base::beg);
  backup_to_start_of_current_record(in);

  string chrom, seq;
  size_t start = 0ul;
  while ((in >> chrom >> start >> seq) &&
         chrom == query_chrom && start < query_end)
    the_reads.push_back(epiread(start, seq));
}

static void
convert_coordinates(const vector<size_t> &cpg_positions,
                    GenomicRegion &region) {

  const size_t start_pos =
    lower_bound(cpg_positions.begin(), cpg_positions.end(),
                region.get_start()) - cpg_positions.begin();

  const size_t end_pos =
    lower_bound(cpg_positions.begin(), cpg_positions.end(),
                region.get_end()) - cpg_positions.begin();

  region.set_start(start_pos);
  region.set_end(end_pos);
}


inline static bool
is_cpg(const string &s, const size_t idx) {
  return toupper(s[idx]) == 'C' && toupper(s[idx + 1]) == 'G';
}


static void
collect_cpgs(const string &s, vector<size_t> &cpgs) {
  const size_t lim = s.length() - 1;
  for (size_t i = 0; i < lim; ++i)
    if (is_cpg(s, i))
      cpgs.push_back(i);
}


static void
get_cpg_positions(const string &chrom_file,
                  vector<size_t> &cpg_positions) {
  vector<string> chrom_names, chrom_seqs;
  read_fasta_file(chrom_file.c_str(), chrom_names, chrom_seqs);
  if (chrom_names.size() > 1)
    throw SMITHLABException("error: more than one seq "
                            "in chrom file" + chrom_file);
  cpg_positions.clear();
  collect_cpgs(chrom_seqs.front(), cpg_positions);
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


static void
clip_reads(const size_t start_pos, const size_t end_pos,
           vector<epiread> &r) {
  size_t j = 0;
  for (size_t i = 0; i < r.size(); ++i) {
    if (start_pos < r[i].pos + r[i].seq.length() &&
        r[i].pos < end_pos) {
      if (r[i].pos < start_pos) {
        assert(start_pos - r[i].pos < r[i].seq.length());
        r[i].seq = r[i].seq.substr(start_pos - r[i].pos);
        r[i].pos = start_pos;
      }
      if (r[i].end() > end_pos)
        r[i].seq = r[i].seq.substr(0, end_pos - r[i].pos);
      r[j] = r[i];
      ++j;
    }
  }
  r.erase(r.begin() + j, r.end());
}


int
main(int argc, const char **argv) {

  try {

    static const string fasta_suffix = "fa";

    bool VERBOSE = false;
    bool PROGRESS = false;
    bool USE_BIC = false;

    string outfile;
    string chroms_dir;

    size_t max_itr = 10;
    double high_prob = 0.75, low_prob = 0.25;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "resolve epi-alleles",
                           "<bed-regions> <mapped-reads>");
    opt_parse.add_opt("output", 'o', "output file", false, outfile);
    opt_parse.add_opt("chrom", 'c', "genome sequence file/directory",
                      true, chroms_dir);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_itr);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("progress", 'P', "print progress info", false, PROGRESS);
    opt_parse.add_opt("bic", 'b', "use BIC to compare models", false, USE_BIC);

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
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    const string regions_file(leftover_args.front());
    const string reads_file_name(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);
    if (!check_sorted(regions))
      throw SMITHLABException("regions not sorted in: " + regions_file);

    unordered_map<string, string> chrom_files;
    identify_chromosomes(chroms_dir, fasta_suffix, chrom_files);

    size_t n_regions  = regions.size();
    if (VERBOSE)
      cerr << "NUMBER OF REGIONS: " << n_regions << endl;

    string curr_chrom;
    vector<size_t> cpg_positions;

    vector<GenomicRegion> amrs;

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    for (size_t i = 0; i < regions.size(); ++i) {
      if (PROGRESS)
        cerr << '\r' << percent(i, n_regions) << "%\r";

      if (regions[i].get_chrom() != curr_chrom) {
        curr_chrom = regions[i].get_chrom();
        const unordered_map<string, string>::const_iterator
          chrom_file(chrom_files.find(curr_chrom));
        if (chrom_file == chrom_files.end())
          throw SMITHLABException("no chrom file for:\n" + toa(regions[i]));
        get_cpg_positions(chrom_file->second, cpg_positions);
      }

      GenomicRegion converted_region(regions[i]);
      convert_coordinates(cpg_positions, converted_region);

      vector<epiread> reads;
      load_reads(reads_file_name, converted_region, reads);

      clip_reads(converted_region.get_start(),
                 converted_region.get_end(), reads);

      if (!reads.empty()) {
        regions[i].set_score((USE_BIC) ?
                             test_asm_bic(max_itr, low_prob, high_prob, reads):
                             test_asm_lrt(max_itr, low_prob, high_prob, reads));
      }
      else regions[i].set_score(1.0);

      regions[i].set_name(regions[i].get_name() + ":" + toa(reads.size()));
      out << regions[i] << endl;
    }
    if (PROGRESS) cerr << "\r100%" << endl;
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
