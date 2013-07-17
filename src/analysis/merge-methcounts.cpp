/*  
 *    conversion in a bisulfite sequencing experiment
 *
 *    Copyright (C) 2009-2012 University of Southern California and
 *                            Andrew D. Smith
 *
 *    Authors: Song Qiang
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
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MappedRead.hpp"
#include "MethpipeFiles.hpp"


using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::accumulate;


int 
main(int argc, const char **argv) 
{
  
  try {
    string outfile("/dev/stdout");
    bool VERBOSE;
        
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "merge multiple methcounts files",
                           "<methcounts-files>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
                      false, outfile);
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
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const vector<string> methcounts_files(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    vector<std::ifstream*> infiles(methcounts_files.size());
    for (size_t i = 0; i < methcounts_files.size(); ++i)
      infiles[i] = new std::ifstream(methcounts_files[i].c_str());
    const bool new_methcount_fmt =
      methpipe::is_methpipe_file_single(methcounts_files.front());
    std::ofstream outf(outfile.c_str());

    string chrom, strand, seq;
    size_t pos, coverage;
    double meth;
    
    while (new_methcount_fmt
           ? methpipe::read_site(*infiles.front(), chrom, pos, strand,
                                 seq, meth, coverage)
           : methpipe::read_site_old(*infiles.front(), chrom, pos, strand,
                           seq, meth, coverage)) {	

      const string ref_chrom = chrom;
      const size_t ref_pos = pos;
      const string ref_strand = strand;
      
      size_t n_total = coverage;
      double n_meth = coverage * meth;
            
      for (size_t i = 1; i < infiles.size(); ++i) {
        if ((new_methcount_fmt
            ? methpipe::read_site(*infiles[i], chrom, pos, strand,
                                  seq, meth, coverage)
            : methpipe::read_site_old(*infiles[i], chrom, pos, strand,
                            seq, meth, coverage))
            && ref_chrom == chrom
            && ref_pos  == pos
            && ref_strand == strand) {
          n_total += coverage;
          n_meth += coverage * meth;
        } else
          throw SMITHLABException("error reading methcount file: "
                                  + methcounts_files[i]);
      }
      if (new_methcount_fmt)
        methpipe::write_site(outf, chrom, pos, strand, seq,
                             n_total == 0 ? 0 : (n_meth / n_total), n_total);
      else
        methpipe::write_site_old(outf, chrom, pos, strand, seq,
                       n_total == 0 ? 0 : (n_meth / n_total), n_total);
    } 
    for (size_t i = 0; i < infiles.size(); ++i)
    {
        infiles[i]->close();
        delete infiles[i];
    }
    outf.close();
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
