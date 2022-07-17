/* amrtester: A program for testing whether a genomic region has
 * allele-specific methylation
 *
 * Copyright (C) 2014-2022 University of Southern California and
 *                         Benjamin E Decato and Andrew D. Smith and Fang Fang
 *
 * Authors: Andrew D. Smith and Benjamin E Decato and Fang Fang
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>

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
using std::runtime_error;
using std::begin;
using std::end;


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
      throw runtime_error("problem loading reads");
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
    throw runtime_error("cannot open input file " + reads_file_name);

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
    lower_bound(begin(cpg_positions), end(cpg_positions),
                region.get_start()) - begin(cpg_positions);

  const size_t end_pos =
    lower_bound(begin(cpg_positions), end(cpg_positions),
                region.get_end()) - begin(cpg_positions);

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
  r.erase(begin(r) + j, end(r));
}


static void
get_chrom(const string &chrom_name,
          const vector<string> &all_chroms,
          const unordered_map<string, size_t> &chrom_lookup,
          string &chrom) {

  auto the_chrom = chrom_lookup.find(chrom_name);
  if (the_chrom == end(chrom_lookup))
    throw runtime_error("could not find chrom: " + chrom_name);

  chrom = all_chroms[the_chrom->second];
  if (chrom.empty())
    throw runtime_error("problem with chrom: " + chrom_name);
}


int
main(int argc, const char **argv) {

  try {

    static const string fasta_suffix = "fa";

    bool VERBOSE = false;
    bool PROGRESS = false;
    bool USE_BIC = false;

    string outfile;
    string chrom_file;

    size_t max_itr = 10;
    double high_prob = 0.75, low_prob = 0.25;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "resolve epi-alleles",
                           "<bed-regions> <mapped-reads>");
    opt_parse.add_opt("output", 'o', "output file", false, outfile);
    opt_parse.add_opt("chrom", 'c', "genome sequence file/directory",
                      true, chrom_file);
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

    vector<string> chrom_files;
    if (isdir(chrom_file.c_str()))
      read_dir(chrom_file, fasta_suffix, chrom_files);
    else chrom_files.push_back(chrom_file);

    if (VERBOSE)
      cerr << "n_chrom_files:\t" << chrom_files.size() << endl;

    /* first load in all the chromosome sequences and names, and make
       a map from chromosome name to the location of the chromosome
       itself */
    vector<string> all_chroms;
    vector<string> chrom_names;
    unordered_map<string, size_t> chrom_lookup;
    for (auto i(begin(chrom_files)); i != end(chrom_files); ++i) {
      vector<string> tmp_chroms, tmp_names;
      read_fasta_file_short_names(*i, tmp_names, tmp_chroms);
      for (size_t j = 0; j < tmp_chroms.size(); ++j) {
        if (chrom_lookup.find(tmp_names[j]) != end(chrom_lookup))
          throw runtime_error("repeated chromosome or name: " + tmp_names[j]);
        chrom_names.push_back(tmp_names[j]);
        chrom_lookup[chrom_names.back()] = all_chroms.size();
        all_chroms.push_back("");
        all_chroms.back().swap(tmp_chroms[j]);
      }
    }

    vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);
    if (!check_sorted(regions))
      throw runtime_error("regions not sorted in: " + regions_file);

    size_t n_regions  = regions.size();
    if (VERBOSE)
      cerr << "NUMBER OF REGIONS: " << n_regions << endl;

    string chrom_name;
    string chrom;
    vector<size_t> cpg_positions;

    vector<GenomicRegion> amrs;

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    for (size_t i = 0; i < regions.size(); ++i) {

      if (PROGRESS)
        cerr << '\r' << percent(i, n_regions) << "%\r";

      // get the correct chrom if it has changed
      if (regions[i].get_chrom() != chrom_name) {
        chrom_name = regions[i].get_chrom();
        if (VERBOSE)
          cerr << "processing " << chrom_name << endl;

        get_chrom(chrom_name, all_chroms, chrom_lookup, chrom);

        cpg_positions.clear();
        collect_cpgs(chrom, cpg_positions);
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
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
