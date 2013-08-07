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

/********Below are functions for merging pair-end reads********/
static void
fill_overlap(const bool pos_str, const MappedRead &mr, const size_t start, 
	     const size_t end, const size_t offset, string &seq, string &scr) {
  const size_t a = pos_str ? (start - mr.r.get_start()) : (mr.r.get_end() - end);
  const size_t b = pos_str ? (end -  mr.r.get_start()) : (mr.r.get_end() - start);
  copy(mr.seq.begin() + a, mr.seq.begin() + b, seq.begin() + offset);
  copy(mr.scr.begin() + a, mr.scr.begin() + b, scr.begin() + offset);
}

static bool
merge_mates(const bool VERBOSE, const size_t range, const MappedRead &one, 
	    const MappedRead &two, MappedRead &merged) {
  
  const bool pos_str = one.r.pos_strand();
  const size_t overlap_start = std::max(one.r.get_start(), two.r.get_start());
  const size_t overlap_end = std::min(one.r.get_end(), two.r.get_end());

  const size_t one_left = pos_str ? 
    one.r.get_start() : std::max(overlap_end, one.r.get_start());
  const size_t one_right = 
    pos_str ? std::min(overlap_start, one.r.get_end()) : one.r.get_end();
  
  const size_t two_left = pos_str ? 
    std::max(overlap_end, two.r.get_start()) : two.r.get_start();
  const size_t two_right = pos_str ? 
    two.r.get_end() : std::min(overlap_start, two.r.get_end());
  
  const int len = pos_str ? (two_right - one_left) : (one_right - two_left);
  
  if(len <= 0){
    if(VERBOSE){
      cerr << one << endl;
      cerr << two << endl;
      cerr << "len = " << len << endl;
    } 
    return false;
  }
  assert(one_left <= one_right && two_left <= two_right);
  assert(overlap_start >= overlap_end || static_cast<size_t>(len) == 
	 ((one_right - one_left) + (two_right - two_left) + (overlap_end - overlap_start)));
  
  string seq(len, 'N');
  string scr(len, 'B');
  if (len > 0 && len <= static_cast<int>(range)) {
    // lim_one: offset in merged sequence where overlap starts
    const size_t lim_one = one_right - one_left;
    copy(one.seq.begin(), one.seq.begin() + lim_one, seq.begin());
    copy(one.scr.begin(), one.scr.begin() + lim_one, scr.begin());
    
    const size_t lim_two = two_right - two_left;
    copy(two.seq.end() - lim_two, two.seq.end(), seq.end() - lim_two);
    copy(two.scr.end() - lim_two, two.scr.end(), scr.end() - lim_two);
    
    // deal with overlapping part
    if (overlap_start < overlap_end) {
      const size_t one_bads = count(one.seq.begin(), one.seq.end(), 'N');
      const int info_one = one.seq.length() - (one_bads + one.r.get_score());
      
      const size_t two_bads = count(two.seq.begin(), two.seq.end(), 'N');
      const int info_two = two.seq.length() - (two_bads + two.r.get_score());
      
      // use the mate with the most info to fill in the overlap
      if (info_one >= info_two)
	fill_overlap(pos_str, one, overlap_start, overlap_end, lim_one, seq, scr);
      else
	fill_overlap(pos_str, two, overlap_start, overlap_end, lim_one, seq, scr);
    }
  }
  
  merged = one;
  merged.r.set_start(pos_str ? one.r.get_start() : two.r.get_start());
  merged.r.set_end(merged.r.get_start() + len);
  merged.r.set_score(one.r.get_score() + two.r.get_score());
  merged.seq = seq;
  merged.scr = scr;  
  const string name(one.r.get_name());
  merged.r.set_name("FRAG:" + name.substr(0, name.size()-2));

  return true;
}

inline static bool
same_read(string mapper, 
	  const MappedRead &a, const MappedRead &b) {
  size_t suffix_len = 1;
  //if (mapper.compare("bismark") == 0) suffix_len = 0;
  //else if (mapper.compare("bs_seeker") == 0) suffix_len = 0;
  const string sa(a.r.get_name());
  const string sb(b.r.get_name());
  return std::equal(sa.begin(), sa.end() - suffix_len, sb.begin());
}
/********Above are functions for merging pair-end reads********/

int 
main(int argc, const char **argv) {
  
  try {
    string outfile;
    string mapper;
    bool bam_format = false;
    size_t MAX_FRAG_LENGTH = 500;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
        "Supported mappers: bsmap, bismark, bs_seeker", "[options] input_file");
    opt_parse.add_opt("output", 'o', "Name of output file.", 
		      true, outfile);
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

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    if (!bam_format) {
      ifstream in;
      in.open(mapped_reads_file.c_str());
      string line;

      //skip header
      do {
        std::getline(in, line);
      }
      while (line.substr(0,1) == "@");

      SAM r1(mapper,line), r2(mapper);
      do {
        MappedRead mr;
        if (r1.is_pairend() && r1.is_mapped()
           && r1.is_primary() ) {
          if (r1.is_mapping_paired()) {
            in >> r2;
            MappedRead mr_1, mr_2;
            mr_1 = r1.GetMappedRead();
            mr_2 = r2.GetMappedRead();
            if(!same_read(mapper,mr_1, mr_2)) {
              cerr << mr_1 << endl;
              cerr << mr_2 << endl;
              throw SMITHLABException("Reads not sorted by name");
            }
            bool merge_success = merge_mates(false, MAX_FRAG_LENGTH,
                mr_1, mr_2, mr);
            if(merge_success)
              out << mr << endl;
            else {
              /********Non-concordant mates are discarded********
              cerr << mr_1 << endl;
              cerr << mr_2 << endl;
              throw SMITHLABException("Problem merging mates");
              */
            }
          }
          else {
            out << r1 << endl;
          }
        }
        else if(r1.is_mapped() && r1.is_primary()) {
          out << r1 << endl;
        }
        // if the read is not mapped or is not primary alignment, do nothing

        in >> r1;
      }
      while (!in.eof());
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
      
      BamAlignment bam_1, bam_2;
      while (reader.GetNextAlignment(bam_1)) {
        MappedRead mr;
        if(bam_1.IsPaired() && bam_1.IsMapped()
           && bam_1.IsPrimaryAlignment()) {
	  if(bam_1.IsProperPair()) {
	    reader.GetNextAlignment(bam_2);
	    // if the next read is not the mate, they will not have the same name
	    MappedRead mr_1, mr_2;
	    BamAlignmentToMappedReadWithMapper(chrom_lookup, bam_1, mr_1, mapper);
	    BamAlignmentToMappedReadWithMapper(chrom_lookup, bam_2, mr_2, mapper);
	    if(!same_read(mapper,mr_1, mr_2)) {
	      cerr << mr_1 << endl;
	      cerr << mr_2 << endl;
	      throw SMITHLABException("Reads not sorted by name");
	    }
            bool merge_success = merge_mates(false, MAX_FRAG_LENGTH,
                mr_1, mr_2, mr);
	    if(merge_success)
	      out << mr << endl;
	    else {
              /********Non-concordant mates are discarded********
	      cerr << mr_1 << endl;
	      cerr << mr_2 << endl;
	      throw SMITHLABException("Problem merging mates");
              */
	    }
          }
	  else {
	    BamAlignmentToMappedReadWithMapper(chrom_lookup, bam_1, mr, mapper);
	    out << mr << endl;
	  }
        }
        else if(bam_1.IsMapped() && bam_1.IsPrimaryAlignment()){
          BamAlignmentToMappedReadWithMapper(chrom_lookup, bam_1, mr, mapper);
          out << mr << endl;
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
