/*    merge-methcounts: a program for merging methcounts files
 *
 *    Copyright (C) 2011-2014 University of Southern California and
 *                            Andrew D. Smith
 *
 *    Authors: Benjamin E Decato, Meng Zhou and Andrew D Smith
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

#include <cmath>

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
#include "MethpipeFiles.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::accumulate;
using std::round;


struct Site {
  string chrom;
  size_t pos;
  string strand;
  string seq;
  double meth;
  size_t coverage;

  Site() {}
  Site(const string &chr, const size_t &position,
       const string &str, const string &sequ,
       const double &met, const size_t &cov) :
    chrom(chr), pos(position), strand(str), seq(sequ),
    meth(met), coverage(cov) {}
};


static std::istream &
read_site(std::istream &in, string &chrom,
          size_t &pos, string &strand, string &seq,
          double &meth, size_t &coverage) {
  return methpipe::read_site(in, chrom, pos, strand,
                             seq, meth, coverage);
}

struct SiteLocationLessThan {
  bool operator()(const Site &a, const Site &b) {
    return a.chrom < b.chrom ||
      (a.chrom == b.chrom && a.pos < b.pos);
  }
};


struct SiteLocationEqual {
  bool operator()(const Site &a, const Site &b) {
    return a.chrom == b.chrom && a.pos == b.pos;
  }
};

//checks if any files still have data to read in
static bool
any_files_are_good(vector<std::ifstream*> infiles){
  for(size_t i=0; i< infiles.size(); ++i)
    if(!(*infiles[i]).eof())return true;
  return false;
}

//updates outdated sites
static bool
load_sites(vector<std::ifstream*> &infiles,
           vector<bool> &outdated, vector<Site> &sites) {
  bool sites_loaded = false;
  for (size_t i=0; i<sites.size(); ++i){
    if (outdated[i]){
      if(read_site(*infiles[i],
                   sites[i].chrom, sites[i].pos, sites[i].strand,
                   sites[i].seq, sites[i].meth, sites[i].coverage)){
        outdated[i]=false;
        sites_loaded = true;
        if ((*infiles[i]).fail()) sites_loaded= false;
      }
    }
  }
  return sites_loaded;
}

//finds first not outdated site
static size_t
find_first_site(vector<bool> &outdated) {
  size_t first_site_pos = std::numeric_limits<size_t>::max();
  for (size_t i = 0; i < outdated.size(); ++i) {
    if (!outdated[i]) {
      first_site_pos = i;
    }
  }
  return first_site_pos;
}

static void
find_minimum_site_location( vector<Site> &sites,
                           vector<bool> &outdated, Site &min_site){
  SiteLocationLessThan comparator; // bad name
  size_t index;
  index = find_first_site(outdated);
  min_site = sites[index];

  for (size_t i=0; i< sites.size(); ++i){
    if(!outdated[i]){
      if (comparator(sites[i], min_site))
        min_site = sites[i];
    }
  }
}

static void
collect_equivalent_locations(Site &min_site,
               vector<Site> &sites,vector<bool> &sites_to_print){
  SiteLocationEqual comparator;
  for(size_t i=0; i < sites.size(); ++i){
    if (comparator(sites[i],min_site))
       sites_to_print[i] = true;
  }
}

static string
format_line_for_tabular(Site &min_site, vector<bool> &to_print,
                        vector<Site> &sites){
  std::ostringstream oss;

  if (*min_site.seq.rbegin() == 'x'){
    min_site. seq = min_site.seq.substr(0,min_site.seq.size()-1);
  }

  oss<< min_site.chrom << ':' << min_site.pos << ':' << min_site.strand
     << ':' << min_site.seq << '\t';
  for (size_t i = 0; i < sites.size(); ++i){
    if (to_print[i]){
      size_t total_meth = round((sites[i].meth)*(sites[i].coverage));
      oss<< sites[i].coverage << '\t' << total_meth  << '\t';
    }
    else oss<< 0 << '\t' << 0 << '\t';
  }
  return oss.str();
}

static string
format_line_for_merged_counts(Site &min_site, vector<bool> &to_print,
                         vector<Site> &sites){
  size_t meth_sum=0;
  size_t cov_sum=0;
  std::ostringstream oss;

  if (*min_site.seq.rbegin() == 'x'){
    min_site. seq = min_site.seq.substr(0,min_site.seq.size()-1);
  }

  oss<< min_site.chrom << '\t'<< min_site.pos << '\t'<< min_site.strand
     << '\t'<< min_site.seq << '\t';

  for(size_t i = 0; i < sites.size(); ++i){
    if (to_print[i]){
      meth_sum += round(sites[i].meth*sites[i].coverage);
      cov_sum += sites[i].coverage;
    }
  }

  double percent;
  if(cov_sum != 0 ) percent =(double)meth_sum/(double)cov_sum;
  else percent = 0;

  oss << percent << '\t' << cov_sum;
  return oss.str();
}

static string
remove_extension(const std::string &filename){
  size_t last_dot = filename.find_last_of(".");
  if (last_dot == std::string::npos) return filename;
  else return filename.substr(0, last_dot);
}

int
main(int argc, const char **argv) {

  try {

    string outfile;
    bool VERBOSE;
    bool TABULAR = false;

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
    vector<string> methcounts_files(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    vector<std::ifstream*> infiles(methcounts_files.size());
    for (size_t i = 0; i < methcounts_files.size(); ++i)
      infiles[i] = new std::ifstream(methcounts_files[i].c_str());

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    for (size_t i = 0; i < methcounts_files.size(); i++){
      methcounts_files[i] = remove_extension(methcounts_files[i]);
    }

    // Print the header if the user specifies or if the output is to
    // be in tabular format
    if (!TABULAR && !header_info.empty())
      out << "#" << header_info << endl;

    if (TABULAR) {
      // tabular format does not include the '#' character
      transform(methcounts_files.begin(), methcounts_files.end(),
                std::ostream_iterator<string>(out, "\t"),
                std::ptr_fun(&strip_path));
      out << endl;
    }

    vector<Site> sites;
    vector<bool> outdated(infiles.size(), true);

    for (size_t i = 0; i< infiles.size(); ++i){ // initialize site vector
      Site new_site;
      sites.push_back(new_site);
    }

    while (any_files_are_good(infiles) &&
           load_sites(infiles, outdated, sites)) {
      Site min_site;
      // find minimum site location
      find_minimum_site_location(sites, outdated, min_site);

      // collect equivalent locations to minimum
      vector<bool> sites_to_print(sites.size(), false);
      collect_equivalent_locations(min_site, sites, sites_to_print);

      // output the appropriate sites' data
      out << ((TABULAR) ?
              format_line_for_tabular(min_site, sites_to_print, sites) :
              format_line_for_merged_counts(min_site, sites_to_print, sites))
          << endl;

      for (size_t i = 0; i < outdated.size(); ++i)
        outdated[i] = (outdated[i] || sites_to_print[i]);
    }

    while (any_files_are_good(infiles)) {
      Site min_site;
      find_minimum_site_location(sites, outdated, min_site);

      vector<bool> sites_to_print(sites.size(), false);
      collect_equivalent_locations(min_site, sites, sites_to_print);

      out << ((TABULAR) ?
              format_line_for_tabular(min_site, sites_to_print, sites) :
              format_line_for_merged_counts(min_site, sites_to_print, sites))
          << endl;

      for (size_t i = 0; i < outdated.size(); ++i)
        outdated[i] = (outdated[i] || sites_to_print[i]);
      load_sites(infiles, outdated, sites);
    }

    for (size_t i = 0; i < infiles.size(); ++i) {
      infiles[i]->close();
      delete infiles[i];
    }
  }
  catch (const SMITHLABException &e)  {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
