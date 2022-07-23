/* merge-methcounts: a program for merging methcounts files
 *
 * Copyright (C) 2011-2022 University of Southern California and
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
#include <unordered_set>
#include <unordered_map>
#include <queue>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "MethpipeSite.hpp"
#include "zlib_wrapper.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::runtime_error;
using std::numeric_limits;
using std::unordered_set;
using std::unordered_map;

static void
set_invalid(MSite &s) {s.pos = numeric_limits<size_t>::max();}

static bool
is_valid(const MSite &s) {return s.pos != numeric_limits<size_t>::max();}

static bool
any_sites_unprocessed(const vector<string> &filenames,
                      const vector<igzfstream*> &infiles,
                      vector<bool> &outdated, vector<MSite> &sites) {

  const size_t n_files = sites.size();

  bool sites_remain = false;
  for (size_t i = 0; i < n_files; ++i) {
    if (outdated[i]) {
      outdated[i] = false;
      MSite tmp_site;
      if ((*infiles[i]) >> tmp_site) {
        // ADS: chrom order within a file already tested
        if (tmp_site.pos <= sites[i].pos && tmp_site.chrom == sites[i].chrom)
          throw runtime_error("sites not sorted in " + filenames[i]);
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


static bool
site_precedes(const vector<MSite> &sites,
              const unordered_map<string, size_t> &chroms_order,
              const size_t a, const size_t b) {
  if (chroms_order.empty())
    return sites[a] < sites[b];

  const size_t a_chr = chroms_order.find(sites[a].chrom)->second;
  const size_t b_chr = chroms_order.find(sites[b].chrom)->second;

  return ((a_chr < b_chr) ||
          (a_chr == b_chr && sites[a].pos < sites[b].pos));
}


static size_t
find_minimum_site(const vector<MSite> &sites,
                  const unordered_map<string, size_t> &chroms_order,
                  const vector<bool> &outdated) {

  const size_t n_files = sites.size();

  // ms_id is the id of the minimum site, the next one to print
  size_t ms_id = numeric_limits<size_t>::max();

  // make sure there is at least one remaining site to print
  for (size_t i = 0; i < n_files && ms_id == numeric_limits<size_t>::max(); ++i)
    if (is_valid(sites[i]) && !outdated[i])
      ms_id = i;

  if (ms_id == numeric_limits<size_t>::max())
    throw runtime_error("failed in a next site to print");

  // now find the earliest site to print among those that could be
  for (size_t i = 0; i < n_files; ++i)
    if (!outdated[i] && is_valid(sites[i]) &&
        site_precedes(sites, chroms_order, i, ms_id))
      ms_id = i;

  return ms_id;
}


static bool
same_location(const MSite &a, const MSite &b) {
  return a.chrom == b.chrom && a.pos == b.pos;
}

static size_t
collect_sites_to_print(const vector<MSite> &sites,
                       const unordered_map<string, size_t> &chroms_order,
                       const vector<bool> &outdated,
                       vector<bool> &to_print) {

  const size_t n_files = sites.size();

  const size_t min_site_idx = find_minimum_site(sites, chroms_order, outdated);

  for (size_t i = 0; i < n_files; ++i)
    // condition below covers "is_valid(sites[i])"
    if (same_location(sites[min_site_idx], sites[i]))
      to_print[i] = true;

  return min_site_idx;
}


static void
write_line_for_tabular(const bool write_fractional,
                       const size_t min_reads,
                       std::ostream &out,
                       const vector<bool> &to_print,
                       const vector<MSite> &sites,
                       MSite min_site) {

  const size_t n_files = sites.size();

  min_site.set_unmutated();

  // ADS: is this the format we want for the row names?
  out << min_site.chrom << ':'
      << min_site.pos << ':'
      << min_site.strand << ':'
      << min_site.context;

  if (write_fractional) {
    for (size_t i = 0; i < n_files; ++i) {
      const size_t r = sites[i].n_reads;
      if (to_print[i] && r > min_reads)
        out << '\t' << sites[i].meth;
      else
        out << '\t' << "NA";
    }
  }
  else
    for (size_t i = 0; i < n_files; ++i) {
      if (to_print[i])
        out << '\t' << sites[i].n_reads << '\t' << sites[i].n_meth();
      else out << '\t' << 0 << '\t'<< 0;
    }
  out << '\n';
}

static void
write_line_for_merged_counts(std::ostream &out,
                             const vector<bool> &to_print,
                             const vector<MSite> &sites,
                             MSite min_site) {

  const size_t n_files = sites.size();

  min_site.set_unmutated();

  double meth_sum = 0;
  min_site.n_reads = 0;
  for (size_t i = 0; i < n_files; ++i)
    if (to_print[i]) {
      meth_sum += sites[i].n_meth();
      min_site.n_reads += sites[i].n_reads;
    }
  min_site.meth = meth_sum/std::max(1ul, min_site.n_reads);

  out << min_site << '\n';
}


static string
remove_extension(const std::string &filename) {
  const size_t last_dot = filename.find_last_of(".");
  if (last_dot == std::string::npos) return filename;
  else return filename.substr(0, last_dot);
}


static string
remove_suffix(const string &suffix, const std::string &filename) {
  if (filename.substr(filename.size() - suffix.size(), suffix.size()) == suffix)
    return filename.substr(0, filename.size() - suffix.size());
  return filename;
}


static void
get_orders_by_file(const string &filename, vector<string> &chroms_order) {

  igzfstream in(filename);
  if (!in)
    throw runtime_error("bad file: " + filename);
  chroms_order.clear();

  unordered_set<string> chroms_seen;
  string line;
  string prev_chrom;

  while (getline(in, line)) {
    line.resize(line.find_first_of(" \t"));
    if (line != prev_chrom) {
      if (chroms_seen.find(line) != end(chroms_seen))
        throw runtime_error("chroms out of order in: " + filename);
      chroms_seen.insert(line);
      chroms_order.push_back(line);
      std::swap(line, prev_chrom);
    }
  }
}


static void
get_chroms_order(const vector<string> &filenames,
                unordered_map<string, size_t> &chroms_order) {

  // get order of chroms in each file
  vector<vector<string>> orders(filenames.size());
  for (size_t i = 0; i < filenames.size(); ++i)
    get_orders_by_file(filenames[i], orders[i]);

  // get the union of chrom sets
  unordered_set<string> the_union;
  for (size_t i = 0; i < orders.size(); ++i)
    for (size_t j = 0; j < orders[i].size(); ++j)
      the_union.insert(orders[i][j]);

  // get an adjacency list and in-degree for each node
  unordered_map<string, vector<string>> adj_lists;
  unordered_map<string, size_t> in_degree;
  for (auto &&i : the_union) {
    in_degree[i] = 0;
    adj_lists[i] = vector<string>();
  }
  for (auto &&i : orders) {
    auto j = begin(i);
    for (auto k = j + 1; k != end(i); ++j, ++k) {
      adj_lists[*j].push_back(*k);
      ++in_degree[*k];
    }
  }

  std::queue<string> q; // invariant: nodes with no incoming edge
  for (auto &&i : the_union)
    if (in_degree[i] == 0)
      q.push(i);

  while (!q.empty()) {
    const string u = q.front();
    q.pop();

    // iterate over the edges (u, v)
    for (auto &&v : adj_lists[u]) {
      --in_degree[v]; // degree should not appear here as 0
      if (in_degree[v] == 0) // this should only happen once per v
        q.push(v);
    }
    adj_lists[u].clear(); // delete node; already had in_degree 0

    chroms_order.insert(make_pair(u, chroms_order.size()));
  }

  // finally, make sure we found a consistent order
  for (auto &&i : adj_lists)
    if (!i.second.empty())
      throw runtime_error("inconsistent order of chroms between files");
}



/*
  This utility does two things, and they are grouped together here
  because of how they are done, not because the uses are related. (1)
  merge-methcounts can take a set of methcounts output files and
  combine them into one. There are several reasons a user might want
  to do this. An example is when technical replicates are performed,
  and analyzed separately to understand technical variance (e.g.,
  between sequencing runs or library preps). After examining the
  technical variation, subsequent analyses might be best conducted on
  all the data together. So all the methcounts files can be combined
  into one using merge-methcounts. In this case, the coverage at any
  site is the sum of the coverages in the original methcounts files,
  and the methylation level at any site is the weighted mean. (2)
  merge-methcounts can take a set of methcounts output files, and
  create a table that contains all the same information. The table
  format is helpful if subsequent analyses are best done using a data
  table, for example in R. When producing a tabular format,
  merge-methcounts allows the user to select whether the desired
  output is in counts or fractions.
 */
int
main(int argc, const char **argv) {

  try {

    static const string description = "merge multiple methcounts files";

    string outfile;
    bool VERBOSE;
    bool write_tabular_format = false;
    bool write_fractional = false;
    bool ignore_chroms_order = false;

    string header_info;
    string header_prefix;
    string column_name_suffix = "RM";
    string suffix_to_remove;

    size_t min_reads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<methcounts-files>");
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("header", 'h',"header to print (ignored for tabular)",
                      false, header_info);
    opt_parse.add_opt("tabular", 't', "output as table",
                      false, write_tabular_format);
    opt_parse.add_opt("prefix", 'p', "prefix header with character",
                      false, header_prefix);
    opt_parse.add_opt("remove", '\0', "Suffix to remove from filenames when "
                      "making column names for tabular format. If not "
                      "specified, suffix including from final dot is removed.",
                      false, suffix_to_remove);
    opt_parse.add_opt("suff", 's',
                      "column name suffixes, one for total reads and one for "
                      "methylated reads, to be separated from sample name "
                      "with underscore",
                      false, column_name_suffix);
    opt_parse.add_opt("fractional", 'f', "output fractions (requires tabular)",
                      false, write_fractional);
    opt_parse.add_opt("reads", 'r', "min reads (for fractional)",
                      false, min_reads);
    opt_parse.add_opt("ignore", '\0',"Do not attempt to determine chromosome. "
                      "Lexicographic order will be used.",
                      false, ignore_chroms_order);
    opt_parse.add_opt("verbose", 'v',"print more run info", false, VERBOSE);
    opt_parse.set_show_defaults();
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
    if (write_fractional && !write_tabular_format) {
      cerr << "fractional output only available for tabular format" << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (column_name_suffix.size() != 2) {
      cerr << "column name suffix must be 2 letters" << endl;
      return EXIT_SUCCESS;
    }
    vector<string> meth_files(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    unordered_map<string, size_t> chroms_order;
    if (!ignore_chroms_order) {
      if (VERBOSE)
        cerr << "resolving chromosome order" << endl;
      get_chroms_order(meth_files, chroms_order);
      if (VERBOSE) {
        cerr << "chromosome order" << endl;
        vector<string> v(chroms_order.size());
        for (auto &&i : chroms_order)
          v[i.second] = i.first;
        for (auto &&i : v)
          cerr << i << endl;
      }
    }

    const size_t n_files = meth_files.size();

    vector<igzfstream*> infiles(n_files);
    for (size_t i = 0; i < n_files; ++i) {
      infiles[i] = new igzfstream(meth_files[i]);
      if (!(*infiles[i]))
        throw runtime_error("cannot open file: " + meth_files[i]);
    }

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile);
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    // print header if user specifies or if tabular output format
    if (write_tabular_format) {

      vector<string> colnames;
      for (auto &&i : meth_files)
        colnames.push_back(strip_path(i));

      for (auto &&i : colnames)
        i = suffix_to_remove.empty() ?
          remove_extension(i) : remove_suffix(suffix_to_remove, i);

      if (!write_fractional) {
        vector<string> tmp;
        for (auto &&i : colnames) {
          tmp.push_back(i + "_" + column_name_suffix[0]);
          tmp.push_back(i + "_" + column_name_suffix[1]);
        }
        swap(tmp, colnames);
      }

      copy(begin(colnames), end(colnames),
           std::ostream_iterator<string>(out, "\t"));
      out << endl;
    }
    else if (!header_info.empty())
      out << "#" << header_info << endl;

    vector<MSite> sites(n_files);
    vector<bool> outdated(n_files, true);
    vector<bool> sites_to_print; // declared here to keep allocation
    vector<unordered_set<string> > chroms_seen(n_files);

    while (any_sites_unprocessed(meth_files, infiles, outdated, sites)) {

      sites_to_print.clear();
      sites_to_print.resize(n_files, false);

      // below idx is one index among the sites to print
      const size_t idx =
        collect_sites_to_print(sites, chroms_order, outdated, sites_to_print);

      // output the appropriate sites' data
      if (write_tabular_format)
        write_line_for_tabular(write_fractional, min_reads, out,
                               sites_to_print, sites, sites[idx]);
      else
        write_line_for_merged_counts(out, sites_to_print, sites, sites[idx]);

      swap(outdated, sites_to_print);
    }

    for (size_t i = 0; i < n_files; ++i)
      delete infiles[i];
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
