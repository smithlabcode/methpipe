/* selectsites: program to select sites, specified in a methcounts
 * format file, that are contained in given (bed format) intervals
 *
 * Copyright (C) 2019 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
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
#include "MethpipeSite.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ios_base;
using std::runtime_error;
using std::ifstream;


static void
collapsebed(vector<GenomicRegion> &regions) {
  size_t j = 0;
  for (size_t i = 1; i < regions.size(); ++i) {
    if (regions[j].same_chrom(regions[i]) &&
        regions[i].get_start() <= regions[j].get_end()) {
      regions[j].set_end(std::max(regions[j].get_end(), regions[i].get_end()));
    }
    else {
      regions[++j] = regions[i];
    }
  }
  regions.erase(begin(regions) + j + 1, end(regions));
}

static bool
precedes(const GenomicRegion &r, const MSite &s) {
  return (r.get_chrom() < s.chrom ||
          (r.get_chrom() == s.chrom && r.get_end() <= s.pos));
}


static bool
contains(const GenomicRegion &r, const MSite &s) {
  return (r.get_chrom() == s.chrom &&
          (r.get_start() <= s.pos && s.pos < r.get_end()));
}


static void
process_all_cpgs(const bool VERBOSE,
                 const string &cpgs_file,
                 vector<GenomicRegion> &regions,
                 std::ostream &out) {

  ifstream in(cpgs_file);
  if (!in)
    throw runtime_error("cannot open file: " + cpgs_file);

  MSite the_site;
  size_t i = 0;
  while (in >> the_site) {
    while (i < regions.size() && precedes(regions[i], the_site))
      ++i;

    if (contains(regions[i], the_site))
      out << the_site << '\n';
  }
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///  CODE BELOW HERE IS FOR SEARCHING ON DISK
static void
move_to_start_of_line(ifstream &in) {
  char next;
  while (in.good() && in.get(next) && next != '\n') {
    in.unget();
    in.unget();
  }
  if (in.bad())
    // hope this only happens when hitting the start of the file
    in.clear();
}


static void
find_start_line(const string &chr, const size_t idx, ifstream &cpg_in) {

  cpg_in.seekg(0, ios_base::beg);
  const size_t begin_pos = cpg_in.tellg();
  cpg_in.seekg(0, ios_base::end);
  const size_t end_pos = cpg_in.tellg();

  if (end_pos - begin_pos < 2)
    throw runtime_error("empty meth file");

  size_t step_size = (end_pos - begin_pos)/2;

  cpg_in.seekg(0, ios_base::beg);
  string low_chr;
  size_t low_idx = 0;
  if (!(cpg_in >> low_chr >> low_idx))
    throw runtime_error("failed navigating inside file");

  // MAGIC: need the -2 here to get past the EOF and possibly a '\n'
  cpg_in.seekg(-2, ios_base::end);
  move_to_start_of_line(cpg_in);
  string high_chr;
  size_t high_idx;
  if (!(cpg_in >> high_chr >> high_idx))
    throw runtime_error("failed navigating inside file");

  size_t pos = step_size;
  cpg_in.seekg(pos, ios_base::beg);
  move_to_start_of_line(cpg_in);

  while (step_size > 0) {
    string mid_chr;
    size_t mid_idx = 0;
    if (!(cpg_in >> mid_chr >> mid_idx))
      throw runtime_error("failed navigating inside file");
    step_size /= 2;
    if (chr < mid_chr || (chr == mid_chr && idx <= mid_idx)) {
      std::swap(mid_chr, high_chr);
      std::swap(mid_idx, high_idx);
      pos -= step_size;
    }
    else {
      std::swap(mid_chr, low_chr);
      std::swap(mid_idx, low_idx);
      pos += step_size;
    }
    cpg_in.seekg(pos, ios_base::beg);
    move_to_start_of_line(cpg_in);
  }
}

static void
get_cpgs_in_region(ifstream &cpg_in, const GenomicRegion &region,
                   std::ostream &out) {

  string chrom(region.get_chrom());
  const size_t start_pos = region.get_start();
  const size_t end_pos = region.get_end();
  find_start_line(chrom, start_pos, cpg_in);

  MSite the_site;
  while (cpg_in >> the_site && (the_site.chrom == chrom &&
                                (the_site.pos < end_pos)))
    if (start_pos <= the_site.pos)
      out << the_site << endl;
}


static void
process_with_cpgs_on_disk(const string &cpgs_file,
                          vector<GenomicRegion> &regions,
                          std::ostream &out) {

  ifstream in(cpgs_file);
  if (!in)
    throw runtime_error("cannot open file: " + cpgs_file);

  for (size_t i = 0; i < regions.size() && in; ++i)
    get_cpgs_in_region(in, regions[i], out);
}


///
///  END OF CODE FOR SEARCHING ON DISK
///
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    bool LOAD_ENTIRE_FILE = false;

    string outfile;

    const string description =
      "Select sites inside a set of genomic intervals. "
      "Sites must be specified in methcounts format. "
      "Intervals must be specified in bed format.";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<intervals-bed> <cpgs-bed>", 2);
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("preload", 'p',
                      "preload sites (use if target intervals very large)",
                      false, LOAD_ENTIRE_FILE);
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
    const string regions_file = leftover_args.front();
    const string cpgs_file = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);
    if (!check_sorted(regions))
      throw runtime_error("regions not sorted in file: " + regions_file);

    const size_t n_orig_regions = regions.size();
    collapsebed(regions);
    if (VERBOSE && n_orig_regions != regions.size())
      cerr << "[number of regions merged due to overlap: "
           << n_orig_regions - regions.size() << "]" << endl;

    ifstream in(cpgs_file);
    if (!in)
      throw runtime_error("cannot open file: " + cpgs_file);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    if (!outfile.empty() && !out)
      throw runtime_error("failed to open output file: " + outfile);

    if (LOAD_ENTIRE_FILE)
      process_all_cpgs(VERBOSE, cpgs_file, regions, out);
    else
      process_with_cpgs_on_disk(cpgs_file, regions, out);
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
