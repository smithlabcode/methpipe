/*    lift-filter: process lift results
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Jenny Qu
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

int
main(int argc, const char **argv) {
  try{
    string pfile;
    bool VERBOSE = false;
    bool UNIQUE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "Process duplicated sites from fast-liftover output",
                           "<methcount file>");
    opt_parse.add_opt("output", 'o', "Output processed methcount", true, pfile);
    opt_parse.add_opt("unique", 'u', "keep unique sites", false, UNIQUE);
    opt_parse.add_opt("verbose", 'v', "print more information", false, VERBOSE);

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
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string mfile(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "Loading methcount file " << mfile << endl;

    std::ifstream in(mfile.c_str());
    std::ofstream output(pfile.c_str());

    string chrom;
    size_t pos;
    string strand;
    string seq;
    double meth;
    size_t coverage;

    vector<GenomicSiteMeth>  newmeth;
    size_t i = 0;
    while (methpipe::read_site(in, chrom, pos, strand,
                               seq, meth, coverage)) {
      const GenomicSiteMeth loc(chrom, pos, strand, seq, meth, coverage);
      if (i==0) {
        newmeth.push_back(loc);
        ++i;
      }
      else if(newmeth.back().chrom == chrom &&
              newmeth.back().pos == pos &&
              newmeth.back().strand == strand){
        if (!UNIQUE) {
          newmeth[i].meth = (newmeth[i].meth*newmeth[i].coverage +
                             meth*coverage)/(newmeth[i].coverage + coverage);
          newmeth[i].coverage =  newmeth[i].coverage + coverage;
        }
      }
      else {
        newmeth.push_back(loc);
        ++i;
      }
    }

    if(VERBOSE)
      cerr << "Keeping " << i << " sites" << endl;

    for(size_t j=0; j < newmeth.size(); ++j){
      methpipe::write_site(output, newmeth[j].chrom,
                           newmeth[j].pos, newmeth[j].strand,
                           newmeth[j].name, newmeth[j].meth,
                           newmeth[j].coverage);
    }
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
