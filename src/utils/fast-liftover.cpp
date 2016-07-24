/*    fast-liftover: lift over sites using index file
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Jenny Qu, Qiang Song
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


/*
Sample indexfile line:
[T-chr] [T-start] [T-end] [S-chr]:[S-start]:[S-end]:[S-strand] [] [T-strand]
chr21   26608683        26608684        chr1:3007015:3007016:-  0       +
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>


#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "GenomicRegion.hpp"
#include "MethpipeFiles.hpp"


using std::string;
using std::ios_base;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::unordered_map;

struct GenomicSite {
  string chrom;
  size_t pos;
  string strand;
  string name;
  GenomicSite(const string &c = "",
              const size_t p = 0,
              const string &s = ""):
    chrom(c), pos(p), strand(s) {}
};

struct GenomicSiteMeth {
  string chrom;
  size_t pos;
  string strand;
  string name;
  double meth;
  size_t coverage;
  GenomicSiteMeth(const string &c = "",
                  const size_t p = 0,
                  const string &s = "",
                  const string &n = "",
                  const double m = 0,
                  const size_t cov = 0):
    chrom(c), pos(p), strand(s), name(n), meth(m), coverage(cov) {}
};


void
flip_strand(GenomicSite &site){
  if (site.strand == "-"){
    site.pos--;
    site.strand = "+";
  }
  return;
}

typedef unordered_map<string,
                      unordered_map<size_t, GenomicSite> > liftover_index;

static void
read_index_file(const bool SS, const string &indexFile,
                unordered_map<string, liftover_index> &index) {
  std::ifstream in(indexFile.c_str());
  if (!in)
    throw SMITHLABException("problem opening index file");

  size_t toPos, toEnd, toScore;
  string toChrom, fromName, toStrand;
  while (in >> toChrom >> toPos >> toEnd >> fromName >> toScore >> toStrand){
    const size_t dim1 = fromName.find_first_of(":");
    const size_t dim2 = fromName.find(":", dim1+1);
    const size_t dim3 = fromName.find(":", dim2+1);
    const string chrom = fromName.substr(0, dim1);
    const size_t pos = atoi(fromName.substr(dim1 + 1, dim2-dim1).c_str());
    const string strand = fromName.substr(dim3 + 1);
    GenomicSite site(toChrom, toPos, toStrand);
    if (SS) flip_strand(site);
    index[chrom][strand][pos] = site;
  }
}

int
main(int argc, const char **argv) {
  try{
    string indexfile;
    string tofile;
    string fromfile;
    string leftfile;

    bool VERBOSE = false;
    bool SS = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "Fast liftOver-all cytosine-by strand" );
    opt_parse.add_opt("indexfile", 'i', "index file", true, indexfile);
    opt_parse.add_opt("from", 'f', "Original file", true, fromfile);
    opt_parse.add_opt("to", 't', "Output file liftovered", true, tofile);
    opt_parse.add_opt("unmapped", 'u', "(optional) File for unmapped sites",
                      false, leftfile);
    opt_parse.add_opt("plus-strand", 'p', "(optional) Report sites on + strand",
                      false, SS);
    opt_parse.add_opt("verbose", 'v', "(optional) Print more information",
                      false, VERBOSE);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
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
    /****************** END COMMAND LINE OPTIONS *****************/

    unordered_map<string, liftover_index> index;
    if (VERBOSE)
      cerr << "Loading index file " << indexfile << endl;
    read_index_file(SS, indexfile, index);

    std::ifstream from(fromfile.c_str());
    std::ofstream to(tofile.c_str());
    std::ofstream unmapped;
    if (!leftfile.empty()) unmapped.open(leftfile.c_str());

    string chrom;
    size_t pos;
    string strand;
    string seq;
    double meth;
    size_t coverage;

    size_t total = 0;
    size_t good = 0;
    size_t nogood = 0;

    if (VERBOSE)
      cerr << "Lifting " << fromfile << " to " << tofile << endl;

    while (methpipe::read_site(from, chrom, pos, strand,
                               seq, meth, coverage)) {
      ++total;
      GenomicSite loc = index[chrom][strand][pos];
      if (!loc.chrom.empty()){
        methpipe::write_site(to, loc.chrom, loc.pos, loc.strand,
                             seq, meth, coverage);
        ++good;
      }
      else {
        if (unmapped.good())
          methpipe::write_site(unmapped, chrom, pos, strand,
                               seq, meth, coverage);
        ++nogood;
        index[chrom][strand].erase(pos);
      }
    }

    if (VERBOSE)
      cerr << "Total sites: " << total << ";\tMapped: "
           << good << ";\tUnmapped: " << nogood << endl;
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
