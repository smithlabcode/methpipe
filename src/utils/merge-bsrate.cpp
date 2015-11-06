/*    merge-bsrate: a program for merging multiple bsrate output files
 *                  to create composite bsrate statistics for biological
 *                  replicates or reference methylomes
 *
 *    Copyright (C) 2011-2014 University of Southern California and
 *                            Andrew D Smith
 *
 *    Authors: Benjamin E Decato
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
#include <iterator>
#include <string>
#include <iomanip>
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MappedRead.hpp"

#include "bsutils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::stringstream;

bool readline(std::vector<std::ifstream*>& infiles,
              std::vector<string>& cur_line) {
  for ( size_t i = 0; i < infiles.size(); ++i) {
    if (infiles[i]->eof() )
      return false;
    else
      getline(*infiles[i],cur_line[i]);
    if(!cur_line[i].compare(""))
      return false;
  }
  return true;
}

int 
main(int argc, const char **argv) {
  
  try {    
    bool VERBOSE = false;
    string outfile;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "Program to merge the "
			   "BS conversion rate from two sets of BS-seq "
			   "reads mapped to a genome",
			   "<bsrate file>, ..., <bsrate file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args; // list of mapped-read files to merge
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
    /****************** END COMMAND LINE OPTIONS *****************/
    vector<std::ifstream*> infiles(leftover_args.size());
    for (size_t i = 0; i < leftover_args.size(); ++i) {
      infiles[i] = new std::ifstream(leftover_args[i].c_str());
      if (!infiles[i])
        throw SMITHLABException("cannot open input file " + leftover_args[i]);
    }

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    static const size_t precision_val = 5;
    out.precision(precision_val);
    // produce weighted average for each length
    vector<string> cur_line(infiles.size());
    vector<string> overall_line(infiles.size());
    vector<string> pos_line(infiles.size());
    vector<string> neg_line(infiles.size());
    vector<string> title_line(infiles.size());

    double overall_conversion_rate = 0;
    double pos_conv_rate = 0;
    double neg_conv_rate = 0;
    size_t sum_bth = 0ul;
    size_t sum_pos = 0ul;
    size_t sum_neg = 0ul;
    size_t base = 1;
    for (size_t i = 0; i < infiles.size(); ++i)
    {
      getline(*infiles[i],overall_line[i]);
      getline(*infiles[i],pos_line[i]);
      getline(*infiles[i],neg_line[i]);
      getline(*infiles[i],title_line[i]);
    }

    vector<string> ostrings;
    ostrings.clear();

    while(readline(infiles, cur_line)) {
      // declare all values
      vector<double> p_total(cur_line.size());
      vector<double> n_total(cur_line.size());
      vector<double> bth_total(cur_line.size());

      vector<double> p_conv(cur_line.size());
      vector<double> n_conv(cur_line.size());
      vector<double> bth_conv(cur_line.size());

      vector<double> p_rate(cur_line.size());
      vector<double> n_rate(cur_line.size());
      vector<double> bth_rate(cur_line.size());

      vector<double> err(cur_line.size());
      vector<double> all(cur_line.size());
      vector<double> err_rate(cur_line.size());

      for (size_t j=0; j< cur_line.size(); ++j) {
        //parse the line
        stringstream ss(cur_line[j]);     
        string item;
        vector<string> elems;
        while(getline(ss,item,'\t')){
          elems.push_back(item);
        }
        p_total[j] = strtod(elems[1].c_str(), NULL);
        p_conv[j] = strtod(elems[2].c_str(), NULL);
        p_rate[j] = strtod(elems[3].c_str(), NULL);

        n_total[j] = strtod(elems[4].c_str(), NULL);
        n_conv[j] = strtod(elems[5].c_str(), NULL);
        n_rate[j] = strtod(elems[6].c_str(), NULL);

        bth_total[j] = strtod(elems[7].c_str(), NULL);
        bth_conv[j] = strtod(elems[8].c_str(), NULL);
        bth_rate[j] = strtod(elems[9].c_str(), NULL);

        err[j] = strtod(elems[10].c_str(), NULL);
        all[j] = strtod(elems[11].c_str(), NULL);
        err_rate[j] = strtod(elems[12].c_str(), NULL);
      }
      size_t ptot_out = 0, ntot_out = 0, bthtot_out = 0, pconv_out = 0;
      size_t nconv_out = 0, bthconv_out = 0;
      size_t err_out = 0, all_out = 0; 
      double prate_out = 0, nrate_out = 0, bthrate_out = 0, errrate_out = 0;

      for (size_t k=0; k<p_total.size(); ++k) {
        prate_out += p_rate[k]*(p_total[k]);
        nrate_out += n_rate[k]*(n_total[k]);
        bthrate_out += bth_rate[k]*(bth_total[k]);
        errrate_out += err_rate[k]*(all[k]);
      }
      ptot_out = accumulate(p_total.begin(), p_total.end(), 0.0);
      ntot_out = accumulate(n_total.begin(), n_total.end(), 0.0);
      bthtot_out = accumulate(bth_total.begin(), bth_total.end(), 0.0);
      all_out = accumulate(all.begin(), all.end(), 0.0);

      pconv_out = accumulate(p_conv.begin(), p_conv.end(), 0.0);
      nconv_out = accumulate(n_conv.begin(), n_conv.end(), 0.0);
      bthconv_out = accumulate(bth_conv.begin(), bth_conv.end(), 0.0);
      err_out = accumulate(err.begin(), err.end(), 0.0);

      prate_out /= ptot_out;
      nrate_out /= ntot_out;
      bthrate_out /= bthtot_out;
      errrate_out /= all_out;
      std::ostringstream x;
      x.precision(precision_val); 
      x << base << "\t" << ptot_out << "\t" << pconv_out << "\t";
      x << setw(precision_val) << prate_out << "\t";
      x << ntot_out << "\t" << nconv_out << "\t";
      x << setw(precision_val) << nrate_out << "\t";
      x << bthtot_out << "\t" << bthconv_out << "\t";
      x << std::setw(precision_val) << bthrate_out << "\t";
      x << all_out << "\t" << err_out << "\t";
      x << setw(precision_val) << errrate_out << endl;
      ostrings.push_back(x.str());

      overall_conversion_rate += bthrate_out*bthtot_out;
      pos_conv_rate += prate_out*ptot_out;
      neg_conv_rate += nrate_out*ntot_out;
      sum_bth += bthtot_out;
      sum_pos += ptot_out;
      sum_neg += ntot_out;
      ++base;
    } 

    out << "OVERALL CONVERSION RATE = ";
    out << setw(precision_val) << overall_conversion_rate/sum_bth << endl;
    out << "POS CONVERSION RATE = ";
    out << setw(precision_val) << pos_conv_rate/sum_pos << "\t";
    out << sum_pos << endl << "NEG CONVERSION RATE = ";
    out << setw(precision_val) << neg_conv_rate/sum_neg << "\t";
    out << sum_neg << endl;

    out << "BASE" << '\t'
    << "PTOT" << '\t'
    << "PCONV" << '\t'
    << "PRATE" << '\t'
    << "NTOT" << '\t'
    << "NCONV" << '\t'
    << "NRATE" << '\t'
    << "BTHTOT" << '\t'
    << "BTHCONV" << '\t'
    << "BTHRATE" << '\t'
    << "ERR" << '\t'
    << "ALL" << '\t'
    << "ERRRATE"  << endl;

    for(size_t i = 0; i < ostrings.size(); ++i) {
      out << ostrings[i];
    }

    for (size_t i = 0; i < infiles.size(); ++i) {
      infiles[i]->close();
      delete infiles[i];
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
