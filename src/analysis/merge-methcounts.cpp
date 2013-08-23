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

struct MethStat {

  MethStat() : 
    total_sites(0), total_covered(0),
    max_cov(0), sum_cov(0), sum_cov_Cs(0) {}
  
  string tostring() const;
  
  void collect(const size_t meth_count, const size_t total) {
    total_sites++;
    if (total > 0) {
      total_covered++;
      max_cov = max(max_cov, total);
      sum_cov += total;
      sum_cov_Cs += meth_count;
    }
  }
  
  size_t total_sites;
  size_t total_covered;
  size_t max_cov;
  size_t sum_cov;
  size_t sum_cov_Cs;
};


string
MethStat::tostring() const {
  std::ostringstream out;
  
  out << "SITES:\t" << total_sites << endl
      << "SITES COVERED:\t" << total_covered << endl
      << "FRACTION:\t" << static_cast<double>(total_covered)/total_sites << endl;
  
  const double overall_cov = 
    static_cast<double>(sum_cov)/max(static_cast<size_t>(1), total_sites);
  const double covered_cov = 
    static_cast<double>(sum_cov)/max(static_cast<size_t>(1), total_covered);
  out << "MAX COVERAGE:\t" << max_cov << endl
      << "MEAN COVERAGE:\t" << overall_cov << endl
      << "MEAN (WHEN > 0):\t" << covered_cov << endl;
  
  const double meth_level = 
    static_cast<double>(sum_cov_Cs)/max(static_cast<size_t>(1), sum_cov);
  out << "MEAN METHYLATION:\t" << meth_level;
  return out.str();
}


std::ostream& 
operator<<(std::ostream& the_stream, const MethStat& ms) {
  return the_stream << ms.tostring();
}

int 
main(int argc, const char **argv) 
{
  
  try {
    string outfile("/dev/stdout");
    string out_stat;
    bool VERBOSE;
        
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "merge multiple methcounts files",
                           "<methcounts-files>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
                      false, outfile);
    opt_parse.add_opt("output_stat", 'S', "Name of output file with statistics",
                      false , out_stat);
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

    MethStat meth_stat_collector;
    
    while (new_methcount_fmt
           ? methpipe::read_site(*infiles.front(), chrom, pos, strand,
                                 seq, meth, coverage)
           : methpipe::read_site_old(*infiles.front(), chrom, pos, strand,
                           seq, meth, coverage)) {	

      const string ref_chrom = chrom;
      const size_t ref_pos = pos;
      const string ref_strand = strand;
      
      size_t n_total = coverage;
      double tmp_meth = coverage * meth + 0.5;
      double n_meth = static_cast<size_t>(tmp_meth);
            
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
          tmp_meth = coverage * meth + 0.5;
          n_meth += static_cast<size_t>(tmp_meth);
        } else
          throw SMITHLABException("error reading methcount file: "
                                  + methcounts_files[i]);
      }

      meth_stat_collector.collect(n_meth, n_total);

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
    if (VERBOSE || !out_stat.empty()) {
      std::ofstream of;
      if (!out_stat.empty()) of.open(out_stat.c_str());
      std::ostream out(out_stat.empty() ? cerr.rdbuf() : of.rdbuf());
      out << meth_stat_collector << endl;
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
