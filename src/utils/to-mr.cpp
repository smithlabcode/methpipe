/*    to-mr: a program for converting SAM and BAM format to MappedRead
 *    format.
 *    Currently supported mappers: bsmap, bismark.
 *
 *    Copyright (C) 2009-2012 University of Southern California and
 *                            Andrew D. Smith
 *
 *    Authors: Meng Zhou, Qiang Song, Andrew Smith
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
#include <fstream>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "SAM.hpp"

#ifdef HAVE_BAMTOOLS
#include <tr1/unordered_map>

#include "bamtools_interface.hpp"
#include <api/BamReader.h>
#include <api/BamAlignment.h>

using std::tr1::unordered_map;
using BamTools::BamAlignment;
using BamTools::SamHeader;
using BamTools::RefVector;
using BamTools::BamReader;
using BamTools::RefData;
#endif

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;

int 
main(int argc, const char **argv) {
  
  try {
    string outfile_t, outfile_a;
    string mapper;
    bool bam_format = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
        "Supported mappers: bsmap, bismark, bs_seeker", "[options] input_file");
    opt_parse.add_opt("t_rich_output", 't', "Name of output file for T-rich reads.", 
		      true, outfile_t);
    opt_parse.add_opt("a_rich_output", 'a', "Name of output file for A-rich \
                      reads, if input is pair-end.", 
		      false, outfile_a);
    opt_parse.add_opt("bam", 'b', "Input file format is bam. Must have bamtools installed."
		      , false, bam_format);
    opt_parse.add_opt("mapper", 'm', "Mapper used to generate input file. See \
                      supported mappers below. If you don't know which mapper \
                      was used, simply type 'unknown'.", 
		      true, mapper);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc < 3 || opt_parse.help_requested()) {
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
    const string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ofstream of_t,of_a;
    of_t.open(outfile_t.c_str());
    of_a.open(outfile_a.c_str());
    std::ostream out_t(of_t.rdbuf());
    std::ostream out_a(of_a.rdbuf());

    if (!bam_format) {
      ifstream in;
      in.open(mapped_reads_file.c_str());
      string line;

      //skip header
      do {
        std::getline(in, line);
      }
      while (line.substr(0,1) == "@");

      SAM r1(mapper,line);
      MappedRead mr;
      if(r1.is_pairend() && r1.is_mapped()
         && r1.is_primary() ) {
        if (r1.is_Trich())
            out_t << r1 << endl;
        else
          out_a << r1 << endl;
      }
      else if(r1.is_mapped() && r1.is_primary()) {
        out_t << r1 << endl;
      }

      while (!in.eof()) {
        in >> r1;
        in.peek();
        MappedRead mr;

        if(r1.is_pairend() && r1.is_mapped()
           && r1.is_primary() ) {
          if (r1.is_Trich())
              out_t << r1 << endl;
          else
            out_a << r1 << endl;
        }
        else if(r1.is_mapped() && r1.is_primary()) {
          out_t << r1 << endl;
        }

        // if the read is not mapped or is not primary alignment, do nothing
      }
    }

    #ifdef HAVE_BAMTOOLS
    else if (bam_format) {
      BamReader reader;
      reader.Open(mapped_reads_file);
      
      // Get header and reference
      string header = reader.GetHeaderText();
      RefVector refs = reader.GetReferenceData();
      
      unordered_map<size_t, string> chrom_lookup;
      for (size_t i = 0; i < refs.size(); ++i)
        chrom_lookup[i] = refs[i].RefName;
      
      BamAlignment bam_1;
      while (reader.GetNextAlignment(bam_1)) {
        MappedRead mr;
        if(bam_1.IsPaired() && bam_1.IsMapped()
           && bam_1.IsPrimaryAlignment()){
          BamAlignmentToMappedReadWithMapper(chrom_lookup, bam_1, mr, mapper);
          if(bam_1.IsFirstMate()){
            out_t << mr << endl;
          }
          else{
            out_a << mr << endl;
          }
        }
        else if(bam_1.IsMapped() && bam_1.IsPrimaryAlignment()){
          BamAlignmentToMappedReadWithMapper(chrom_lookup, bam_1, mr, mapper);
          out_t << mr << endl;
        }

        // if the read is not mapped or is not primary alignment, do nothing

      }
      reader.Close();
    }
    #endif
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
