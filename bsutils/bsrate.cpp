/*    bsrate: a program for determining the rate of bisulfite
 *    conversion in a bisulfite sequencing experiment
 *
 *    Copyright (C) 2009 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
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

#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"

#include "BSUtils.hpp"

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;

template <class U, class V> static void
separate_regions(const vector<U> &regions, const vector<V> &seqs,
		 vector<vector<U> > &sep_regions, vector<vector<V> > &sep_seqs) {
  const size_t n_regions = regions.size();
  assert(n_regions <= seqs.size());
  if (regions.empty()) return;
  
  string prev_chrom = regions.front().get_chrom();
  sep_regions.resize(1);
  sep_seqs.resize(1);
  for (size_t i = 0; i < regions.size(); ++i) {
    std::string current_chrom(regions[i].get_chrom());
    if (current_chrom != prev_chrom) {
      prev_chrom.swap(current_chrom);
      sep_regions.push_back(vector<U>());
      sep_seqs.push_back(vector<V>());
    }
    sep_regions.back().push_back(regions[i]);
    sep_seqs.back().push_back(seqs[i]);
  }
}


static void
check_conversionfun(const string &chrom, 
		    const vector<GenomicRegion> locations,
		    const vector<string> sequences,
		    size_t &total_pos, size_t &conv_pos, 
		    size_t &total_neg, size_t &conv_neg) {

  const size_t read_width = locations.front().get_width();
  const size_t chrom_size = chrom.length();

  for (size_t i = 0; i < locations.size(); ++i) {
    assert(locations[i].get_end() < chrom_size);
    const size_t offset = locations[i].get_start();
    if (locations[i].pos_strand()) {
      for (size_t j = 0; j < read_width - 1; ++j) {
	if (offset + j + 1 < chrom_size) {
	  if (is_cytosine(chrom[offset + j]) && 
	      !is_guanine(chrom[offset + j + 1])) {
	    total_pos += (is_cytosine(sequences[i][j]) || 
			  is_thymine(sequences[i][j]));
	    conv_pos += is_thymine(sequences[i][j]);
	  }
	}
      }
    }
    else {
      for (size_t j = 0; j < read_width; ++j) {
	if (offset + read_width >= j + 2) {
	  if (is_guanine(chrom[offset + read_width - 1 - j]) && 
	      !is_cytosine(chrom[offset + read_width - 2 - j])) {
	    total_neg += (is_cytosine(sequences[i][j]) || 
			  is_thymine(sequences[i][j]));
	    conv_neg += is_thymine(sequences[i][j]);
	  }
	}
      }
    }
  }
}

int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    string chrom_dir;
    string mapped_reads_file;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("bsrate", "A program for obtaining the conversion "
			   "rate in a bisulfite sequencing experiment");
    opt_parse.add_opt("chrom", 'c', "chromosome directory (FASTA format)", 
		      true , chrom_dir);
    opt_parse.add_opt("mapped", 'm', "file of mapped reads (BED format)", 
		      true , mapped_reads_file);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
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
    const string reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/


    if (VERBOSE)
      cerr << "reading mapped read locations" << endl;
    vector<GenomicRegion> mapped_reads;
    ReadBEDFile(mapped_reads_file, mapped_reads);
    assert(check_sorted(mapped_reads));
    

    if (VERBOSE)
      cerr << "loading read sequences" << endl;
    vector<string> names, sequences;
    read_fasta_file(reads_file.c_str(), names, sequences);
    if (VERBOSE)
      cerr << "read " << names.size() << " sequences" << endl;


    if (VERBOSE)
      cerr << "associating read sequences to mapped locations" << endl;
    vector<size_t> lookup;
    relative_sort(mapped_reads, names, lookup);
    vector<string> seq_swapper(sequences.size());
    for (size_t i = 0; i < lookup.size(); ++i)
      seq_swapper[i].swap(sequences[lookup[i]]);
    sequences.swap(seq_swapper);
    seq_swapper.clear();


    if (VERBOSE)
      cerr << "separating chromosomes" << endl;
    vector<vector<GenomicRegion> > mapped_reads_by_chrom;
    vector<vector<string> > sequences_by_chrom;
    separate_regions(mapped_reads, sequences, 
		     mapped_reads_by_chrom, sequences_by_chrom);
    mapped_reads.clear();
    sequences.clear();


    size_t total_pos = 0, conv_pos = 0;
    size_t total_neg = 0, conv_neg = 0;
    
    for (size_t i = 0; i < mapped_reads_by_chrom.size(); ++i) {
      const string chrom(mapped_reads_by_chrom[i].front().get_chrom());
      const string chrom_file(path_join(chrom_dir, chrom + ".fa"));
      
      if (VERBOSE)
	cerr << "[PROCESSING] " << chrom_file << endl;
      vector<string> dummy_chrom_names, chrom_seqs;
      read_fasta_file(chrom_file.c_str(), dummy_chrom_names, chrom_seqs);
      
      check_conversionfun(chrom_seqs.front(),
			  mapped_reads_by_chrom[i], sequences_by_chrom[i],
			  total_pos, conv_pos, total_neg, conv_neg);
    }

    
    cout << "\tTOTAL\tCONVERTED\tRATE" << endl
	 << "POS:\t" << total_pos << "\t" << conv_pos << "\t" 
	 << conv_pos/static_cast<double>(total_pos) << endl
	 << "NEG:\t" << total_neg << "\t" << conv_neg << "\t" 
	 << conv_neg/static_cast<double>(total_neg) << endl
	 << "BOTH:\t" << total_pos + total_neg << "\t" 
	 << conv_pos + conv_neg << "\t" 
	 << (conv_pos + conv_neg)/
      static_cast<double>(total_pos + total_neg) << endl;
  
  }
  catch (const RMAPException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
