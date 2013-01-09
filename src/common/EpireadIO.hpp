/*
 *    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith and Fang Fang
 *
 *    Authors: Andrew D. Smith and Fang Fang
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

#ifndef EPIREAD_IO
#define EPIREAD_IO

#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <string>
#include <vector>
#include <tr1/unordered_map>

#include "Epiread.hpp"

class EpireadIO {
public:
  EpireadIO(const std::string &rfn,
	    const bool v, const bool ef, const std::string &cd);
  void
  load_reads(const GenomicRegion &region, std::vector<epiread> &the_reads);
  bool
  load_reads_next_chrom(std::string &chrom, std::vector<epiread> &the_reads);
  void
  convert_coordinates(std::vector<GenomicRegion> &amrs) const;
  
  static void
  set_max_fragment_size(const size_t mfs) {max_fragment_size = mfs;}
  
private:
  
  static size_t max_fragment_size;
  
  std::string chrom_name_cache;
  std::string chrom_seq_cache;
  std::string reads_file_name;
  std::tr1::unordered_map<size_t, size_t> cpgs_cache;
  
  size_t position;
  bool VERBOSE;
  bool EPIREAD_FORMAT; // to indicate if Epiread or MappedRead format 
  std::tr1::unordered_map<std::string, std::string> chrom_seq_files;
  std::tr1::unordered_map<std::string, std::vector<size_t> > cpg_rev_lookup;
  
  static std::ios_base::streampos
  find_last_position_before(const std::string &query_chrom, 
			    const size_t query_pos, const size_t range_low, 
			    const size_t range_high, std::ifstream &in);
  
  static std::ios_base::streampos
  find_first_position_after(const std::string &query_chrom,
			    const size_t query_pos,const size_t range_low, 
			    const size_t range_high, std::ifstream &in);
  
  static std::ios_base::streampos
  find_last_position_before_mr(const std::string &query_chrom, 
			       const size_t query_pos, const size_t range_low, 
			       const size_t range_high, std::ifstream &in);
  static std::ios_base::streampos
  find_first_position_after_mr(const std::string &query_chrom, 
			       const size_t query_pos, const size_t range_low, 
			       const size_t range_high, std::ifstream &in);
  
};

#endif
