/* merge-methcounts: a program for merging methcounts files
 *
 * Copyright (C) 2011-2019 University of Southern California and
 *                         Andrew D. Smith
 *
 * Authors: Benjamin E Decato, Meng Zhou and Andrew D Smith
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
 */

#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <stdexcept>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "MethpipeSite.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::runtime_error;
using std::numeric_limits;

static void
unmutate_site(MSite &s) {
  s.context.resize(s.context.size() - 1);
}

static bool
precedes(const MSite &a, const MSite &b) {
  const int c = a.chrom.compare(b.chrom);
  return (c < 0 || (c == 0 && a.pos < b.pos));
}

static bool
same_location(const MSite &a, const MSite &b) {
  return a.chrom == b.chrom && a.pos == b.pos;
}

static void
set_invalid(MSite &s) {
  s.pos = numeric_limits<size_t>::max();
}

static bool
is_valid(const MSite &s) {
  return s.pos != numeric_limits<size_t>::max();
}

static bool
any_sites_unprocessed(const vector<string> &filenames,
                      const vector<std::ifstream*> &infiles,
                      vector<bool> &outdated, vector<MSite> &sites) {

  const size_t n_files = sites.size();

  bool sites_remain = false;
  for (size_t i = 0; i < n_files; ++i) {
    if (outdated[i]) {
      outdated[i] = false;
      MSite tmp_site;
      if (*infiles[i] >> tmp_site) {
        if (precedes(tmp_site, sites[i]))
          throw runtime_error("error: sites not sorted in " + filenames[i]);
        sites_remain = true;
        sites[i] = tmp_site;
      }
      else set_invalid(sites[i]);
    }
    else if (is_valid(sites[i]))
      sites_remain = true;
  }
  return sites_remain;
}


static size_t
find_minimum_site(const vector<MSite> &sites, const vector<bool> &outdated) {

  const size_t n_files = sites.size();
  size_t ms_id = numeric_limits<size_t>::max();

  for (size_t i = 0; i < n_files && ms_id == numeric_limits<size_t>::max(); ++i)
    if (is_valid(sites[i]) && !outdated[i])
      ms_id = i;

  if (ms_id == numeric_limits<size_t>::max())
    throw runtime_error("failed in find_minimum_site");

  for (size_t i = 0; i < n_files; ++i)
    if (!outdated[i] && is_valid(sites[i]) && precedes(sites[i], sites[ms_id]))
      ms_id = i;

  return ms_id;
}


static size_t
collect_sites_to_print(const vector<MSite> &sites, const vector<bool> &outdated,
                       vector<bool> &to_print) {

  const size_t n_files = sites.size();

  const size_t min_site_idx = find_minimum_site(sites, outdated);

  for (size_t i = 0; i < n_files; ++i)
    // condition below covers "is_valid(sites[i])"
    if (same_location(sites[min_site_idx], sites[i]))
      to_print[i] = true;

  return min_site_idx;
}


static void
write_line_for_tabular(const bool write_fractional,
                       std::ostream &out,
                       const vector<bool> &to_print,
                       const vector<MSite> &sites,
                       MSite min_site) {

  const size_t n_files = sites.size();

  if (min_site.is_mutated())
    unmutate_site(min_site);

  out << min_site.chrom << ':'
      << min_site.pos << ':'
      << min_site.strand << ':'
      << min_site.context;

  if (write_fractional) {
    for (size_t i = 0; i < n_files; ++i) {
      if (to_print[i])
        out << '\t' << sites[i].meth;
      else
        out << '\t' << 0;
    }
  }
  else {
    for (size_t i = 0; i < n_files; ++i) {
      if (to_print[i])
        out << '\t' << sites[i].n_reads << '\t' << sites[i].n_meth();
      else
        out << '\t' << 0 << '\t'<< 0;
    }
  }
  out << '\n';
}

static void
write_line_for_merged_counts(std::ostream &out,
                             const vector<bool> &to_print,
                             const vector<MSite> &sites,
                             MSite min_site) {

  const size_t n_files = sites.size();

  if (min_site.is_mutated())
    unmutate_site(min_site);

  size_t meth_sum = 0;
  min_site.n_reads = 0;
  for (size_t i = 0; i < n_files; ++i)
    if (to_print[i]) {
      meth_sum += sites[i].n_meth();
      min_site.n_reads += sites[i].n_reads;
    }
  min_site.meth = static_cast<double>(meth_sum)/min_site.n_reads;

  out << min_site << '\n';
}

static string
remove_extension(const std::string &filename){
  const size_t last_dot = filename.find_last_of(".");
  if (last_dot == std::string::npos) return filename;
  else return filename.substr(0, last_dot);
}

int
main(int argc, const char **argv) {

  try {

    string outfile;
    bool VERBOSE;
    bool TABULAR = false;
    bool FRAC = false;

    string header_info;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "merge multiple methcounts files",
                           "<methcounts-files>");
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("header", 'h',"header to print (ignored for tabular)",
                      false, header_info);
    opt_parse.add_opt("verbose", 'v',"print more run info", false, VERBOSE);
    opt_parse.add_opt("tabular", 't', "output as table", false, TABULAR);
    opt_parse.add_opt("fractional", 'f', "output fractions (requires tabular)",
                      false, FRAC);
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
    if (FRAC && !TABULAR) {
      cerr << "fractional output only available for tabular format" << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    vector<string> meth_files(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    const size_t n_files = meth_files.size();

    vector<std::ifstream*> infiles(n_files);
    for (size_t i = 0; i < n_files; ++i) {
      infiles[i] = new std::ifstream(meth_files[i]);
      if (!(*infiles[i]))
        throw runtime_error("cannot open file: " + meth_files[i]);
    }

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    // print header if user specifies or if tabular output format
    if (TABULAR) {
      // tabular format header does not include '#' character
      vector<string> colnames;
      for (auto &&i : meth_files)
        colnames.push_back(strip_path(remove_extension(i)));
      copy(begin(colnames), end(colnames),
           std::ostream_iterator<string>(out, "\t"));
      out << endl;
    }
    else if (!header_info.empty())
      out << "#" << header_info << endl;

    vector<MSite> sites(n_files);
    vector<bool> outdated(n_files, true);
    vector<bool> sites_to_print; // declared here to keep allocation

    while (any_sites_unprocessed(meth_files, infiles, outdated, sites)) {
    	
      sites_to_print.clear();
      sites_to_print.resize(n_files, false);

      // below idx is one index among the sites to print
      const size_t idx = collect_sites_to_print(sites, outdated, sites_to_print);

      // output the appropriate sites' data
      if (TABULAR)
        write_line_for_tabular(FRAC, out, sites_to_print, sites, sites[idx]);
      else
        write_line_for_merged_counts(out, sites_to_print, sites, sites[idx]);

      swap(outdated, sites_to_print);
    }

    for (size_t i = 0; i < n_files; ++i) {
      infiles[i]->close();
      delete infiles[i];
    }
  }
  catch (const runtime_error &e)  {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
