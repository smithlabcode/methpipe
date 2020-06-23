/* multimethstat: program to summarize methylation values according to
 * genomic intervals specified in a BED format file. The column
 * headings must be the names associated with intervals in a separate
 * BED file. The methylation table must have columns corresponding to
 * sites in the genome
 *
 * Aug 30 2019 (ADS): This program is written for data frames with
 * columns corresponding to features (probes) and intended for EPIC
 * data. It should be modified for the transpose and also for the
 * format accepted by RADMeth.
 *
 * Copyright (C) 2019 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
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

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <iterator>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"

#include <sys/stat.h>

using std::string;
using std::vector;
using std::ifstream;
using std::cerr;
using std::cout;
using std::endl;
using std::unordered_map;
using std::istringstream;
using std::numeric_limits;
using std::unordered_set;

using std::begin;
using std::end;

using std::pair;
using std::make_pair;

using std::greater;
using std::runtime_error;

static void
parse_table_row(const string &row, vector<double> &values) {

  if (row.empty() || !isalnum(row[0]))
    throw runtime_error("rownames must begin with alphanumeric");

  values.clear(); // keep the reserve, but use push_back

  const char* a = &row[row.find_first_of(" \t")];
  char* b = 0;
  double val = strtod(a, &b);

  while (a != b) {
    values.push_back(val);
    a = b;
    val = strtod(a, &b);
  }
}

static void
parse_header(const string &header, vector<string> &colnames) {
  std::istringstream h(header);
  colnames = {std::istream_iterator<string>{h}, {}};
}

// precedes and follows functions are not like > and !<= since they
// must exclude overlap of intervals
static bool
follows(const GenomicRegion &a, const GenomicRegion &b) {
  return (b.get_chrom() < a.get_chrom() ||
          (b.same_chrom(a) && b.get_end() < a.get_start()));
}
static bool
precedes(const GenomicRegion &a, const GenomicRegion &b) {
  return (a.get_chrom() < b.get_chrom() ||
          (a.same_chrom(b) && a.get_end() <= b.get_start()));
}

static void
build_probe_to_frag(const vector<GenomicRegion> &probes,
                    const unordered_map<string, size_t> &probe_names,
                    const vector<string> &colnames,
                    const vector<GenomicRegion> &frags,
                    vector<size_t> &probe_to_frag,
                    vector<double> &n_probes_per_frag) {

  const size_t n_probes = probes.size();
  probe_to_frag = vector<size_t>(n_probes, numeric_limits<size_t>::max());

  const size_t n_frags = frags.size();

  size_t curr_frag = 0;
  for (size_t i = 0; i < n_probes; ++i) {
    while (curr_frag < n_frags && precedes(frags[curr_frag], probes[i]))
      ++curr_frag;
    if (curr_frag < n_frags && !follows(frags[curr_frag], probes[i]))
      probe_to_frag[i] = curr_frag;
  }

  const size_t n_colnames = colnames.size();
  vector<size_t> tmp(n_colnames, numeric_limits<size_t>::max());

  for (size_t i = 0; i < n_colnames; ++i) {
    auto idx(probe_names.find(colnames[i]));
    if (idx != end(probe_names))
      tmp[i] = probe_to_frag[idx->second];
    /* ADS: the line below is commented out because there are "probes"
     * in some data sets that have no associated genomic location in
     * various files kicking around for mapping probes to genome
     * coordinates
     */
    // else throw runtime_error("cannot find column: " + colnames[i]);
  }
  probe_to_frag.swap(tmp);

  n_probes_per_frag = vector<double>(n_frags, 0);
  for (auto &&i : probe_to_frag)
    if (i != numeric_limits<size_t>::max())
      ++n_probes_per_frag[i];
}


static bool
all_names_unique(const vector<GenomicRegion> &regions) {
  unordered_set<string> names_seen;
  for (size_t i = 0; i < regions.size(); ++i) {
    auto idx = names_seen.find(regions[i].get_name());
    if (idx == end(names_seen))
      names_seen.insert(regions[i].get_name());
    else return false;
  }
  return true;
}


static string
wrong_n_vals(const size_t lines_read,
             const vector<double> &probe_values,
             const vector<string> &colnames) {
  std::ostringstream oss;
  oss << "wrong number of values "
      << "(row=" << (lines_read + 1) << ", "
      << "n_vals=" << probe_values.size() << ", "
      << "expect=" << colnames.size() << ")";
  return oss.str();
}


struct end_point {
  end_point(const string c, const size_t s, const bool isf) :
    chr(c), start(s), is_first(isf) {}
  bool operator<(const end_point &other) const {
    return (chr < other.chr ||
            (chr == other.chr &&
             (start < other.start ||
              (start == other.start &&
               is_first < other.is_first))));
  }
  string chr;
  size_t start;
  bool is_first;
};


static void
get_frags(const vector<GenomicRegion> &features,
          vector<GenomicRegion> &frags) {

  vector<end_point> end_points;
  for (auto &&f : features) {
    end_points.push_back(end_point(f.get_chrom(), f.get_start(), true));
    end_points.push_back(end_point(f.get_chrom(), f.get_end(), false));
  }
  sort(begin(end_points), end(end_points));

  size_t count = 0;
  GenomicRegion region;
  for (size_t i = 0; i < end_points.size() - 1; ++i) {
    if (end_points[i].is_first) count++;
    else count--;
    if (count > 0) {
      region.set_chrom(end_points[i].chr);
      region.set_start(end_points[i].start);
      region.set_end(end_points[i + 1].start);
      if (region.get_width() > 0)
        frags.push_back(region);
    }
  }
}


static void
get_frag_to_feature(const vector<GenomicRegion> &features,
                    const vector<GenomicRegion> &frags,
                    vector<vector<size_t> > &frag_to_feature) {

  const size_t n_frags = frags.size();
  frag_to_feature.resize(n_frags);

  size_t j = 0;
  for (size_t i = 0; i < features.size(); ++i) {

    while (j < n_frags && precedes(frags[j], features[i])) ++j;

    for (size_t k = j; k < n_frags && !precedes(features[i], frags[k]); ++k)
      frag_to_feature[k].push_back(i);
  }
}


int
main(int argc, const char **argv) {

  try {

    static const string description =
      "program to summarize methylation values according to genomic intervals"
      " specified in a BED format file. The column headings must be the names"
      " associated with intervals in a separate BED file. The methylation"
      " table must have columns corresponding to sites in the genome";

    string outfile;
    bool VERBOSE = false;
    bool report_empty_intervals = false;
    bool name_by_interval = false;
    bool report_progress = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<intervals-bed> <probes-bed> <meth-data-frame>");
    opt_parse.add_opt("outfile", 'o', "output file", false, outfile);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("progress", '\0', "report progress",
                      false, report_progress);
    opt_parse.add_opt("name-by-interval", '\0',
                      "name features by interval in output",
                      false, name_by_interval);
    opt_parse.add_opt("empty", 'e', "report empty intervals",
                      false, report_empty_intervals);
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
    if (leftover_args.size() != 3) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string features_file(leftover_args.front());
    const string probe_locations_file(leftover_args[1]);
    const string table_file(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    // read the features
    if (VERBOSE)
      cerr << "[reading intervals]" << endl;
    ifstream feature_in(features_file);
    if (!feature_in)
      throw runtime_error("could not open file: " + features_file);

    vector<GenomicRegion> features;
    GenomicRegion feature;
    while (feature_in >> feature)
      features.push_back(feature);

    if (!check_sorted(features))
      throw runtime_error("features not sorted in: " + features_file);

    if (!name_by_interval && !all_names_unique(features))
      throw runtime_error("duplicate feature names in: " + features_file);

    const size_t n_features = features.size();

    if (VERBOSE)
      cerr << "[n regions to summarize: " << n_features << "]" << endl;

    // read the probe locations
    if (VERBOSE)
      cerr << "[reading probe locations]" << endl;
    ifstream probe_in(probe_locations_file);
    if (!probe_in)
      throw runtime_error("could not open file: " + probe_locations_file);

    vector<GenomicRegion> probes;
    GenomicRegion probe;
    while (probe_in >> probe)
      probes.push_back(probe);

    if (!check_sorted(probes))
      throw runtime_error("probes not sorted in: " + probe_locations_file);

    unordered_map<string, size_t> probe_names;
    for (size_t i = 0; i < probes.size(); ++i)
      probe_names[probes[i].get_name()] = i;

    const size_t n_probes = probes.size();
    if (VERBOSE)
      cerr << "[n probes with locations: " << n_probes << "]" << endl;

    // read and process the probe names from the table header
    if (VERBOSE)
      cerr << "[checking probe names in table header]" << endl;
    ifstream in(table_file);
    if (!in)
      throw runtime_error("could not open file: " + table_file);

    string header;
    getline(in, header);
    vector<string> colnames;
    parse_header(header, colnames);

    if (VERBOSE)
      cerr << "[n columns in data matrix: " << colnames.size() << "]" << endl;

    // vector<GenomicRegion> features,
    vector<GenomicRegion> frags;
    get_frags(features, frags);

    vector<vector<size_t> > frag_to_feature;
    get_frag_to_feature(features, frags, frag_to_feature);

    // for each feature, get an index for the corresponding column
    vector<size_t> probe_to_frag;
    vector<double> n_probes_per_frag;
    build_probe_to_frag(probes, probe_names, colnames, frags,
                        probe_to_frag, n_probes_per_frag);

    vector<double> n_probes_per_feature(features.size(), 0.0);
    for (size_t i = 0; i < frag_to_feature.size(); ++i)
      for (size_t j = 0; j < frag_to_feature[i].size(); ++j)
        n_probes_per_feature[frag_to_feature[i][j]] += n_probes_per_frag[i];

    if (VERBOSE)
       cerr << "[processing table by sample]" << endl;
    vector<vector<double> > feature_values;

    // get file size for reporting progress reading file
    struct stat st;
    stat(table_file.c_str(), &st);

    ProgressBar progress(st.st_size);
    if (report_progress)
      progress.report(cerr, 0);

    // read and parse all the lines in the methylation table
    string line;
    size_t lines_read = 0;
    vector<string> sample_names; // name for each row
    vector<double> probe_values(colnames.size());
    while (getline(in, line)) {

      sample_names.push_back(line.substr(0, line.find_first_of(" \t")));

      parse_table_row(line, probe_values);

      if (probe_values.size() != colnames.size())
        throw runtime_error(wrong_n_vals(lines_read, probe_values, colnames));

      vector<double> tmp(n_features, 0.0);
      for (size_t i = 0; i < probe_values.size(); ++i) {
        // get the frag that the i-th probe maps to
        const size_t curr_frag = probe_to_frag[i];
        // if it has a valid feature mapping, then add it's value
        if (curr_frag != numeric_limits<size_t>::max())
          for (auto &&j : frag_to_feature[curr_frag])
            tmp[j] += probe_values[i];
      }

      ++lines_read;
      if (report_progress && progress.time_to_report(in.tellg()))
        progress.report(cerr, in.tellg());

      feature_values.push_back(vector<double>());
      feature_values.back().swap(tmp);
    }
    if (VERBOSE)
      cerr << "[samples read: " << lines_read << "]" << endl;

    // open output stream
    if (VERBOSE)
      cerr << "[writing output]" << endl;
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    // write the header line for the output file
    bool write_tab = false;
    for (size_t i = 0; i < n_features; ++i)
      if (report_empty_intervals || n_probes_per_feature[i] > 0) {
        if (write_tab)
          out << '\t';
        out << (name_by_interval ?
                assemble_region_name(features[i]) :
                features[i].get_name());
        write_tab = true;
      }
    out << endl;

    // write each row of the output file
    for (size_t i = 0; i < lines_read; ++i) {
      out << sample_names[i];
      for (size_t j = 0; j < n_features; ++j)
        if (n_probes_per_feature[j] > 0)
          out << '\t' << feature_values[i][j]/n_probes_per_feature[j];
        else if (report_empty_intervals)
          out << '\t' << "NA";
      out << endl;
    }
  }
  catch (const std::runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
