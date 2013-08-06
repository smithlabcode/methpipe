/*    rmapbs-pe: mapping Illumina BS-seq paired-end reads
 *
 *    Copyright (C) 2012 University of Southern California and
 *                       Andrew D. Smith
 *                       Qiang Song
 *
 *    Authors: Andrew D. Smith and Qiang Song
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
#include <cmath>
#include <tr1/unordered_map>
#include <unistd.h>
#include <ios>

#include "FastRead.hpp"
#include "smithlab_os.hpp"
#include "SeedMaker.hpp"
#include "MapResultForPE.hpp"
#include "MappedRead.hpp"
#include "OptionParser.hpp"
#include "load_reads.hpp"
#include "clip_adaptor_from_reads.hpp"

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::min;
using std::numeric_limits;
using std::pair;
using std::make_pair;
using std::ptr_fun;


static string
format_time(const clock_t start, const clock_t end) {
  const size_t tics = end - start;
  return "[" + toa(static_cast<float>(tics)/CLOCKS_PER_SEC) + " SEC]";
}


static void
bisulfite_treatment(bool AG_WC, size_t &x) {
  if (AG_WC)
    x = ((x & 0x5555555555555555ul) | 
	 (((x & 0x5555555555555555ul) << 1) & 
	  (x & 0xAAAAAAAAAAAAAAAAul)));
  else x |= ((x & 0x5555555555555555ul) << 1);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////
////  THIS STUFF DEALS WITH THE HASHING OF KEYS
////
template <typename T> struct SeedHash {
  typedef
  unordered_map<size_t,
		pair<typename vector<T>::const_iterator,
		     typename vector<T>::const_iterator> > type;
};
typedef vector<pair<size_t, unsigned int> > SeedHashSorter;


static void
load_seeds(const bool VERBOSE, vector<size_t> &the_seeds) {
  /* ORIGINAL SEEDS */
  // 0b1111110011000011111100110000111111001100001111110011000000111100
  the_seeds.push_back(18213668567193169980ul);
  // 0b0011111100110000111111001100001111110011000011111100110000001111
  the_seeds.push_back(4553417141798292495ul);
  // 0b0000111111001100001111110011000011111100110000111111001111111111
  the_seeds.push_back(1138354285449573375ul);
  // 0b1100001111110011000011111100110000111111001100001111110000110011
  the_seeds.push_back(14119646626644556851ul);
  // 0b0011000011111100110000111111001100001111110011000011111111110000
  the_seeds.push_back(3529911656661139440ul);
  // 0b1100110000111111001100001111110011000011111100110000111111001100
  the_seeds.push_back(14717535969447448524ul);
  // 0b1111001100001111110011000011111100110000111111001100001111000011
  the_seeds.push_back(17514442047644025795ul);
  
  // for (size_t i = 0; i < the_seeds.size(); ++i)
  //   the_seeds[i] >>= 8ul;
  
  if (VERBOSE) {
    cerr << endl << "SEED STRUCTURES:" << endl;
    for (size_t i = 0; i < the_seeds.size(); ++i)
      cerr << bits2string_masked(rmap_bits::all_ones, the_seeds[i]) << endl;
    cerr << endl;
  }
}


struct wildcard_score {
  wildcard_score() {}
  template <class T, class U> size_t score_tc(T a, U b) const {
    return a->score_tc_wild_n(b);
  }
  template <class T, class U> size_t score_ag(T a, U b) const {
    return a->score_ag_wild_n(b);
  }
};


inline size_t
bisulfite_treatment_tc(char c) {
  switch(c) {
  case 1 : return 3;
  default : return c;
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////
////  WHERE THE ACTUAL MAPPING HAPPENS
////
template <class T, class U> void
map_reads_tc(const U &specialized_score,
             const string &chrom, const size_t chrom_id,
             const size_t profile, const size_t read_width, 
             const size_t max_diffs, const vector<T> &fast_reads, 
             /*const*/ typename SeedHash<T>::type &seed_hash, 
             const bool strand, vector<MultiMapResult> &best_maps) {

  MASK_t bad_bases = rmap_bits::all_ones;
  MASK_t read_word = rmap_bits::all_zeros;
  T fast_read;
  
  const size_t key_diff = read_width - min(read_width, SeedMaker::max_seed_part);
  const size_t chrom_size = chrom.size();
  
  size_t chrom_offset = 0;
  while(chrom_offset < min(chrom_size, key_diff))
    fast_read.shift(chrom[chrom_offset++]);
  
  while (chrom_offset < min(chrom_size, read_width - 1)) {
    const size_t key_base = chrom[chrom_offset - key_diff];
    fast_read.shift(chrom[chrom_offset++]);
    SeedMaker::update_bad_bases(key_base, bad_bases);
    SeedMaker::update_read_word(key_base, read_word);
  }
  
  string::const_iterator key_pos(chrom.begin() + (chrom_offset - key_diff));
  const string::const_iterator chrom_lim(chrom.end());
  for (string::const_iterator chrom_pos(chrom.begin() + chrom_offset); 
       chrom_pos < chrom_lim; ++chrom_pos, ++key_pos) {
    const size_t key_base = bisulfite_treatment_tc(*key_pos);
    fast_read.shift(*chrom_pos);
    SeedMaker::update_bad_bases(key_base, bad_bases);
    SeedMaker::update_read_word(key_base, read_word);
    if ((bad_bases & profile) == 0) {
      typename SeedHash<T>::type::const_iterator 
	bucket(seed_hash.find(read_word & profile));
      if (bucket != seed_hash.end()) {
	pair<typename vector<T>::const_iterator, 
	  typename vector<T>::const_iterator> tmp(bucket->second);
	const typename vector<T>::const_iterator limit(tmp.second);
	for (typename vector<T>::const_iterator to_test(tmp.first); 
	     to_test != limit; ++to_test) {
	  const size_t score = specialized_score.score_tc(to_test, fast_read);
	  if (score <= max_diffs) {
	    const vector<MultiMapResult>::iterator 
	      current(best_maps.begin() + (to_test - fast_reads.begin()));
	    current->add(score, chrom_id, distance(chrom.begin(), chrom_pos)
			 - read_width + 1, strand);
	  }
	}
      }
    }
  }
}


inline size_t
bisulfite_treatment_ag(char c) {
  switch(c) {
  case 2 : return 0;
  default : return c;
  }
}

template <class T, class U> void
map_reads_ag(const U &specialized_score,
             const string &chrom, const size_t chrom_id,
             const size_t profile, const size_t read_width, 
             const size_t max_diffs, const vector<T> &fast_reads, 
             /*const*/ typename SeedHash<T>::type &seed_hash, 
             const bool strand, vector<MultiMapResult> &best_maps) {

  MASK_t bad_bases = rmap_bits::all_ones;
  MASK_t read_word = rmap_bits::all_zeros;
  T fast_read;
  
  const size_t key_diff = read_width - min(read_width, SeedMaker::max_seed_part);
  const size_t chrom_size = chrom.size();
  
  size_t chrom_offset = 0;
  while(chrom_offset < min(chrom_size, key_diff))
    fast_read.shift(chrom[chrom_offset++]);
  
  while (chrom_offset < min(chrom_size, read_width - 1)) {
    const size_t key_base = chrom[chrom_offset - key_diff];
    fast_read.shift(chrom[chrom_offset++]);
    SeedMaker::update_bad_bases(key_base, bad_bases);
    SeedMaker::update_read_word(key_base, read_word);
  }
  
  string::const_iterator key_pos(chrom.begin() + (chrom_offset - key_diff));
  const string::const_iterator chrom_lim(chrom.end());
  for (string::const_iterator chrom_pos(chrom.begin() + chrom_offset); 
       chrom_pos < chrom_lim; ++chrom_pos, ++key_pos) {
    const size_t key_base = bisulfite_treatment_ag(*key_pos);
    fast_read.shift(*chrom_pos);
    SeedMaker::update_bad_bases(key_base, bad_bases);
    SeedMaker::update_read_word(key_base, read_word);
    if ((bad_bases & profile) == 0) {
      typename SeedHash<T>::type::const_iterator 
	bucket(seed_hash.find(read_word & profile));
      if (bucket != seed_hash.end()) {
	pair<typename vector<T>::const_iterator, 
	  typename vector<T>::const_iterator> tmp(bucket->second);
	const typename vector<T>::const_iterator limit(tmp.second);
	for (typename vector<T>::const_iterator to_test(tmp.first); 
	     to_test != limit; ++to_test) {
	  const size_t score = specialized_score.score_ag(to_test, fast_read);
	  if (score <= max_diffs) {
	    const vector<MultiMapResult>::iterator 
	      current(best_maps.begin() + distance(fast_reads.begin(), to_test));
	    current->add(score, chrom_id, distance(chrom.begin(), chrom_pos)
			 - read_width + 1, strand);
	  }
	}
      }
    }
  }
}


static void
treat_cpgs(const bool AG_WILDCARD, string &chrom) {
  const size_t lim = chrom.length() - (!AG_WILDCARD);
  if (AG_WILDCARD) {
    for (size_t i = 1; i < lim; ++i)
      if (chrom[i] == 0 and chrom[i - 1] == 1) chrom[i] = 2;
  }
  else for (size_t i = 0; i < lim; ++i)
	 if (chrom[i] == 3 and chrom[i + 1] == 2) chrom[i] = 1;
}


static void
get_read_matches(const size_t the_seed, const vector<size_t> &read_words,
                 SeedHashSorter &sh_sorter) {
  const size_t lim = read_words.size();
  sh_sorter.resize(read_words.size());
  for (size_t i = 0; i < lim; ++i)
    sh_sorter[i] = make_pair(the_seed & read_words[i], i);
  sort(sh_sorter.begin(), sh_sorter.end());
}


template <class T> void
sort_by_key(const SeedHashSorter &sh, vector<T> &in) {
  vector<T> tmp(in);
  size_t j = 0;
  for (SeedHashSorter::const_iterator i(sh.begin()); i != sh.end(); ++i, ++j)
    in[j] = tmp[i->second];
}

template <class T> void
swap_sort_by_key(const SeedHashSorter &sh, vector<T> &in) {
  vector<T> tmp(in.size()); 
  size_t j = 0;
  for (SeedHashSorter::const_iterator i(sh.begin()); i != sh.end(); ++i, ++j)
    std::swap(tmp[j], in[i->second]);
  std::swap(tmp, in);
}


template <class T> void
sort_by_key(SeedHashSorter &sh_sorter, vector<MultiMapResult> &best_maps,
            vector<size_t> &reads, vector<unsigned int> &read_index, 
            vector<T> &fast_reads) {
  swap_sort_by_key(sh_sorter, best_maps);
  // sort_by_key(sh_sorter, best_maps);
  sort_by_key(sh_sorter, reads);
  sort_by_key(sh_sorter, read_index);
  sort_by_key(sh_sorter, fast_reads);
  size_t j = 0;
  for (SeedHashSorter::iterator i(sh_sorter.begin()); i != sh_sorter.end(); ++i, ++j)
    i->second = j;
}


template <class T> void
build_seed_hash(const SeedHashSorter &sh_sorter, const vector<T> &fast_reads,
                typename SeedHash<T>::type &seed_hash) {
  seed_hash.clear();
  typename vector<T>::const_iterator frb(fast_reads.begin());
  size_t prev_key = 0, prev_idx = 0, curr_idx = 0;
  for (SeedHashSorter::const_iterator shs(sh_sorter.begin()); 
       shs != sh_sorter.end(); ++shs) {
    curr_idx = shs->second;
    if (shs->first != prev_key) {
      seed_hash[prev_key] = make_pair(frb + prev_idx, frb + curr_idx);
      prev_key = shs->first;
      prev_idx = curr_idx;
    }
  }
  seed_hash[prev_key] = make_pair(frb + prev_idx, fast_reads.end());
}


template <class T> static void
resort_reads(const size_t the_seed,
             vector<T> &fast_reads,  vector<size_t> &read_words, 
             vector<unsigned int> &read_index,
             vector<MultiMapResult> &best_maps,
             typename SeedHash<T>::type &seed_hash) {
  seed_hash.clear();
  SeedHashSorter sh_sorter;
  get_read_matches(the_seed, read_words, sh_sorter);
  sort_by_key(sh_sorter, best_maps, read_words, read_index, fast_reads);
  build_seed_hash(sh_sorter, fast_reads, seed_hash);
}


unsigned char
b2i(char c) {
  switch(c) {
  case 'A' : return 0;
  case 'C' : return 1;
  case 'G' : return 2;
  case 'T' : return 3;
  default : return 4;
  }
}


unsigned char
b2i_rc(char c) {
  switch(c) {
  case 'A' : return 3;
  case 'C' : return 2;
  case 'G' : return 1;
  case 'T' : return 0;
  default : return 4;
  }
}


template <class T, class U> void
iterate_over_seeds(const bool VERBOSE, const bool AG_WILDCARD, 
                   const U &specialized_score,
                   const vector<size_t> &the_seeds, 
                   const vector<string> &chrom_files,
                   vector<string> &chrom_names, vector<size_t> &chrom_sizes,
                   vector<T> &fast_reads,  vector<size_t> &read_words, 
                   vector<unsigned int> &read_index,
                   vector<MultiMapResult> &best_maps,
                   const size_t max_mismatches, const size_t read_width) {
  
  if (VERBOSE)
    cerr << "[SCANNING CHROMOSOMES]" << endl;
  
  for (size_t j = 0; j < the_seeds.size() && !fast_reads.empty(); ++j) {
    if (VERBOSE)
      cerr << "[SEED:" << j + 1 << "/" << the_seeds.size() << "] "
	   << "[FORMATTING READS]" << endl;
    
    typename SeedHash<T>::type seed_hash;
    resort_reads(the_seeds[j], fast_reads, read_words, 
		 read_index, best_maps, seed_hash);
    
    size_t prev_chrom_count = 0;
    for (size_t i = 0; i < chrom_files.size() && !fast_reads.empty(); ++i) {
      
      vector<string> tmp_chrom_names, chroms;
      if (VERBOSE)
	cerr << "[SEED:" << j + 1 << "/" << the_seeds.size() << "] "
	     << "[LOADING CHROM] ";
      read_fasta_file(chrom_files[i].c_str(), tmp_chrom_names, chroms);
      
      if (VERBOSE)
	cerr << "[SCANNING=" << tmp_chrom_names.front() << "] ";
      
      const clock_t start(clock());
      for (size_t k = 0; k < chroms.size(); ++k) {
	if (j == 0) {
	  chrom_sizes.push_back(chroms[k].length());
	  chrom_names.push_back(tmp_chrom_names[k]);
	}
	transform(chroms[k].begin(), chroms[k].end(), chroms[k].begin(), ptr_fun(&toupper));
	string tmp_chrom(chroms[k]);
	transform(chroms[k].begin(), chroms[k].end(), tmp_chrom.begin(), ptr_fun(&b2i));
	treat_cpgs(AG_WILDCARD, tmp_chrom);
	if (AG_WILDCARD)
	  map_reads_ag(specialized_score, 
		       tmp_chrom, prev_chrom_count + k, the_seeds[j], read_width, 
		       max_mismatches, fast_reads, seed_hash, true, best_maps);
	else
	  map_reads_tc(specialized_score, 
		       tmp_chrom, prev_chrom_count + k, the_seeds[j], read_width, 
		       max_mismatches, fast_reads, seed_hash, true, best_maps);
	string().swap(tmp_chrom);
	
	std::reverse(chroms[k].begin(), chroms[k].end());
	transform(chroms[k].begin(), chroms[k].end(), chroms[k].begin(), ptr_fun(&b2i_rc));
	treat_cpgs(AG_WILDCARD, chroms[k]);
	if (AG_WILDCARD)
	  map_reads_ag(specialized_score, 
		       chroms[k], prev_chrom_count + k, the_seeds[j], read_width, 
		       max_mismatches, fast_reads, seed_hash, false, best_maps);
	else
	  map_reads_tc(specialized_score, 
		       chroms[k], prev_chrom_count + k, the_seeds[j], read_width,
		       max_mismatches, fast_reads, seed_hash, false, best_maps);
	string().swap(chroms[k]);
      }
      const clock_t end(clock());
      if (VERBOSE) cerr << format_time(start, end) << endl;
      prev_chrom_count += chroms.size();
    }
    for (size_t x = 0; x < best_maps.size(); ++x)
      best_maps[x].collapse();
  }
  vector<T>().swap(fast_reads);
  vector<size_t>().swap(read_words);
}


static void
identify_chromosomes(const bool VERBOSE, const string fasta_suffix, 
		     const string chrom_file, vector<string> &chrom_files) {
  if (VERBOSE)
    cerr << "[IDENTIFYING CHROMS] ";
  if (isdir(chrom_file.c_str())) 
    read_dir(chrom_file, fasta_suffix, chrom_files);
  else chrom_files.push_back(chrom_file);
  if (VERBOSE) {
    cerr << "[DONE]" << endl 
	 << "chromosome files found (approx size):" << endl;
    for (vector<string>::const_iterator i = chrom_files.begin();
	 i != chrom_files.end(); ++i)
      cerr << *i << " (" << roundf(get_filesize(*i)/1e06) << "Mbp)" << endl;
    cerr << endl;
  }
}


static void
load_reads(const bool VERBOSE, const bool AG_WILDCARD,
           const string &adaptor, const string &reads_file, 
	   const size_t read_start_index, const size_t n_reads_to_process, 
	   vector<FastRead> &fast_reads, vector<unsigned int> &read_index, 
	   vector<size_t> &read_words, size_t &read_width) {
  
  //////////////////////////////////////////////////////////////
  // LOAD THE READS (AS SEQUENCES OR PROBABILITIES) FROM DISK
  if (VERBOSE) cerr << "[LOADING READ SEQUENCES] ";
  load_reads_from_fastq_file(reads_file, read_start_index, n_reads_to_process,
			     adaptor, read_width,
			     fast_reads, read_words, read_index);
  for (size_t i = 0; i < read_words.size(); ++i)
    bisulfite_treatment(AG_WILDCARD, read_words[i]);
  if (VERBOSE)
    cerr << "[DONE]" << endl
	 << "TOTAL HQ READS: " << read_index.size() << endl
	 << "READ WIDTH: " << read_width << endl;
}


struct indexed_best_less {
  bool operator()(const pair<unsigned int, MultiMapResult> &a,
		  const pair<unsigned int, MultiMapResult> &b) const {
    return a.first < b.first;
  }
};


static void
invert_bests_list(vector<unsigned int> &read_index, 
                  vector<MultiMapResult> &bests) {
  vector<pair<unsigned int, size_t> > sorter;
  for (size_t i = 0; i < bests.size(); ++i)
    sorter.push_back(make_pair(read_index[i], i));
  sort(sorter.begin(), sorter.end());
  vector<MultiMapResult> tmp(bests.size());
  for (size_t i = 0; i < sorter.size(); ++i) {
    read_index[i] = sorter[i].first;
    std::swap(tmp[i], bests[sorter[i].second]);
  }
  std::swap(tmp, bests);
}


static size_t
count_mapped(const vector<MultiMapResult> &best_maps) {
  size_t total_mapped = 0;
  for (size_t i = 0; i < best_maps.size(); ++i) {
    const size_t maps = best_maps[i].mr.size();
    total_mapped += (maps > 0);
  }
  return total_mapped;
}


static double
get_average_maps(const vector<MultiMapResult> &best_maps) {
  size_t total_mapped = 0;
  double average_maps = 0.0;
  for (size_t i = 0; i < best_maps.size(); ++i) {
    const size_t maps = best_maps[i].mr.size();
    average_maps += maps;
    total_mapped += (maps > 0);
  }
  return average_maps = average_maps/total_mapped;
}



static void
map_reads(const bool VERBOSE, const bool AG_WILDCARD,
	  const string &adaptor_sequence, const vector<string> &chrom_files, 
	  const string &reads_file, const size_t read_start_index, 
	  const size_t n_reads_to_process, const vector<size_t> the_seeds, 
	  size_t max_mismatches, vector<unsigned int> &read_index,
	  vector<MultiMapResult> &best_maps, size_t &read_width) {
  
  //////////////////////////////////////////////////////////////
  // OBTAIN THE READS
  // 
  vector<FastRead> fast_reads;
  vector<size_t> read_words;
  read_width = 0;
  load_reads(VERBOSE, AG_WILDCARD, adaptor_sequence, 
	     reads_file, read_start_index, n_reads_to_process,
	     fast_reads, read_index, read_words, read_width);
  
  best_maps.resize(read_words.size());    
  
  if (max_mismatches == numeric_limits<size_t>::max())
    max_mismatches = static_cast<size_t>(0.10*read_width);
  
  if (VERBOSE)
    cerr << "MAX MISMATCHES: " << max_mismatches << endl;
  
  vector<size_t> chrom_sizes;
  vector<string> chrom_names;
  
  const wildcard_score specialized_score;
  iterate_over_seeds(VERBOSE, AG_WILDCARD, specialized_score,
		     the_seeds, chrom_files, chrom_names, 
		     chrom_sizes, fast_reads, read_words, read_index,
		     best_maps, max_mismatches, read_width);
  
  invert_bests_list(read_index, best_maps);
  
  for (size_t i = 0; i < chrom_names.size(); ++i) {
    const size_t chr_name_end = chrom_names[i].find_first_of(" \t");
    if (chr_name_end != string::npos)
      chrom_names[i].erase(chr_name_end);
  }
  MapResult::chrom_names = chrom_names;
  MapResult::chrom_sizes = chrom_sizes;

  if (VERBOSE)
    cerr << "MAPPED: " << count_mapped(best_maps) << endl
	 << "AVERAGE MAPS: " << get_average_maps(best_maps) << endl;
}
  

static MappedRead
make_mapped_read(const size_t read_len, const MapResult &mr, const string &name,
		 const string &sequence, const string &scores) {
  const size_t chrom_id = mr.chrom;
  size_t start = mr.strand ? mr.site :
    MapResult::chrom_sizes[chrom_id] - mr.site - read_len;
  size_t end = start + read_len;
  const char strand = mr.strand ? '+' : '-';
  const size_t real_read_len = sequence.length();
  if (strand == '+')
    end = start + real_read_len;
  else
    start = end - real_read_len;
  MappedRead r;
  r.r = GenomicRegion(MapResult::chrom_names[chrom_id], start, end, 
		      name, mr.score, strand);
  r.seq = sequence;
  r.scr = scores;
  return r;
}

static string
read_name_from_fastq(const string &buffer) {
  const size_t truncpos = buffer.find_first_of(" \t");
  return string(buffer.begin() + 1, (truncpos != string::npos) ?
		buffer.begin() + truncpos : buffer.end());
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

static void
fill_overlap(const bool pos_str, const MappedRead &mr, const size_t start, 
	     const size_t end, const size_t offset, string &seq, string &scr) {
  const size_t a = pos_str ? (start - mr.r.get_start()) : (mr.r.get_end() - end);
  const size_t b = pos_str ? (end -  mr.r.get_start()) : (mr.r.get_end() - start);
  copy(mr.seq.begin() + a, mr.seq.begin() + b, seq.begin() + offset);
  copy(mr.scr.begin() + a, mr.scr.begin() + b, scr.begin() + offset);
}

static void
merge_mates(const size_t suffix_len, const size_t range,
	    const MappedRead &one, const MappedRead &two, MappedRead &merged) {
  
  const bool pos_str = one.r.pos_strand();
  const size_t overlap_start = max(one.r.get_start(), two.r.get_start());
  const size_t overlap_end = min(one.r.get_end(), two.r.get_end());

  const size_t one_left = pos_str ? 
    one.r.get_start() : max(overlap_end, one.r.get_start());
  const size_t one_right = 
    pos_str ? min(overlap_start, one.r.get_end()) : one.r.get_end();
  
  const size_t two_left = pos_str ? 
    max(overlap_end, two.r.get_start()) : two.r.get_start();
  const size_t two_right = pos_str ? 
    two.r.get_end() : min(overlap_start, two.r.get_end());
  
  const int len = pos_str ? (two_right - one_left) : (one_right - two_left);
  
  assert(len > 0);
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
  merged.r.set_name("FRAG:" + name.substr(0, name.size() - suffix_len));
}


struct MappedReadLess {
  bool operator()(const MappedRead &a, const MappedRead &b) const {
    return a.r < b.r;
  }
};

static int
get_fragment_length(const MappedRead &a, const MappedRead &b) {
  return (a.r.pos_strand()) ? 
    (b.r.get_end() - a.r.get_start()) : (a.r.get_end() - b.r.get_start());
}


static double
merge_score(const int max_fraglen, const MappedRead &a, const MappedRead &b) {
  if (!a.r.same_chrom(a.r) || a.r.get_strand() != b.r.get_strand())
    return numeric_limits<double>::max();
  int frag_size = get_fragment_length(a, b);
  if (frag_size <= 0 || frag_size > max_fraglen)
    return numeric_limits<double>::max();
  return a.r.get_score() + b.r.get_score(); 
}


static bool
within_range(const MappedRead &a, const MappedRead &b, const size_t d) {
  return static_cast<size_t>(abs(get_fragment_length(a, b))) < d;
}


static bool
before_range(const MappedRead &a, const MappedRead &b, const size_t d) {
  return (a.r.get_end() < b.r.get_start()) && !within_range(a, b, d);
}


static double
find_best_mates(const size_t range, 
		const vector<MappedRead> &one, const vector<MappedRead> &two,
		size_t &one_best, size_t &two_best) {
  double best_score = numeric_limits<double>::max();
  bool ambig = false;
  size_t k = 0;
  for (size_t i = 0; i < one.size(); ++i) {
    while (k < two.size() && before_range(two[k], one[i], range))
      ++k;
    for (size_t j = k; j < two.size() && within_range(two[j], one[i], range); ++j) {
      const double score = merge_score(range, one[i], two[j]);
      if (score < best_score) {
	one_best = i;
	two_best = j;
	best_score = score;
	ambig = false;
      }	
      else if (score == best_score) 
	ambig = true;
    }
  }
  return ambig ? numeric_limits<double>::max() : best_score;
}


static size_t
find_best(const vector<MappedRead> &a) {
  size_t best_score = numeric_limits<size_t>::max();
  size_t best_pos = 0;
  bool ambig = false;
  for (size_t i = 0; i < a.size(); ++i) {
    const size_t score = a[i].r.get_score();
    if (score < best_score) {
      best_pos = i;
      best_score = score;
      ambig = false;
    }
    else if (score == best_score) ambig = true;
  }
  return ambig ? numeric_limits<size_t>::max() : best_pos;
}


static void
process_same_read(const size_t range, const size_t suffix_len, 
		  const MappedReadLess mrl, vector<MappedRead> &one, 
		  vector<MappedRead> &two, std::ostream &out) {
  sort(one.begin(), one.end(), mrl);
  sort(two.begin(), two.end(), mrl);
  
  size_t one_best = 0, two_best = 0;
  const double best_score = 
    find_best_mates(range, one, two, one_best, two_best);
  
  if (best_score < numeric_limits<double>::max()) {
    MappedRead merged;
    merge_mates(suffix_len, range, one[one_best], two[two_best], merged);
    out << merged << '\n';
  }
  else {
    one_best = find_best(one);
    if (one_best != numeric_limits<size_t>::max())
      out << one[one_best] << '\n';
    two_best = find_best(two);
    if (two_best != numeric_limits<size_t>::max())
      out << two[two_best] << '\n';
  }
}


static void
process_single_read(const vector<MappedRead> &mr, std::ostream &out) {
  const size_t best = find_best(mr);
  if (best != numeric_limits<size_t>::max())
    out << mr[best] << '\n';
}


static void
get_mapped_reads(const bool SECOND_END,
		 const size_t read_len, const string &adaptor,
		 const size_t target_read, const vector<MapResult> &r,
		 std::ifstream &in, size_t &line_count,
		 vector<MappedRead> &mr) {
  // move to the correct position
  const size_t target_line = 4*target_read;
  assert(line_count <= target_line);
  
  string buffer;
  while (!in.eof() && line_count < target_line) {
    in.ignore(numeric_limits<std::streamsize>::max(), '\n');
    ++line_count;
  }
  assert(!in.eof());
  
  // load the read info
  getline(in, buffer);
  const string name(read_name_from_fastq(buffer));
  getline(in, buffer);
  string sequence(buffer);
  if (!adaptor.empty())
    clip_adaptor_from_read(adaptor, sequence);
  if (SECOND_END)
    revcomp_inplace(sequence);
  in.ignore(numeric_limits<std::streamsize>::max(), '\n');
  getline(in, buffer);
  string scores(buffer);
  if (SECOND_END)
    reverse(scores.begin(), scores.end());
  
  line_count += 4;
  
  // convert to mapped reads
  mr.clear();
  for (size_t i = 0; i < r.size(); ++i) {
    mr.push_back(make_mapped_read(read_len, r[i], name, sequence, scores));
    if (SECOND_END)
      mr.back().r.set_strand(mr.back().r.pos_strand() ? '-' : '+');
  }
}


static void
get_next(const vector<MultiMapResult> &bests, size_t &idx) {
  ++idx;
  while (idx < bests.size() && bests[idx].empty()) ++idx;
}


static void
clip_and_output(const bool VERBOSE, const size_t range, 
		const size_t read_len, const size_t suffix_len, 
		const string &end_one_file, const string &end_two_file, 
		const string &adaptor_one, const string &adaptor_two, 
		vector<unsigned int> &read_index_one, vector<unsigned int> &read_index_two,
		vector<MultiMapResult> &bests_one, vector<MultiMapResult> &bests_two,
		const string &outfile) {
  
  std::ifstream in_one(end_one_file.c_str());
  std::ifstream in_two(end_two_file.c_str());
  
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
  
  size_t one_idx = 0, two_idx = 0;
  while (one_idx < bests_one.size() && bests_one[one_idx].empty()) ++one_idx;
  while (two_idx < bests_two.size() && bests_two[two_idx].empty()) ++two_idx;

  MappedReadLess mrl;
  size_t line_count_one = 0, line_count_two = 0;
  vector<MappedRead> one, two;
  while (one_idx < bests_one.size() && two_idx < bests_two.size()) {
    
    // if they are the same
    if (read_index_one[one_idx] == read_index_two[two_idx]) {
      get_mapped_reads(false, read_len, adaptor_one, read_index_one[one_idx],
		       bests_one[one_idx].mr, in_one, line_count_one, one);
      get_mapped_reads(true, read_len, adaptor_two, read_index_two[two_idx],
		       bests_two[two_idx].mr, in_two, line_count_two, two);
      process_same_read(range, suffix_len, mrl, one, two, out);
      get_next(bests_one, one_idx);
      get_next(bests_two, two_idx);
    }
    else if (read_index_one[one_idx] < read_index_two[two_idx]) {
      get_mapped_reads(false, read_len, adaptor_one, read_index_one[one_idx],
		       bests_one[one_idx].mr, in_one, line_count_one, one);
      process_single_read(one, out);
      get_next(bests_one, one_idx);
    }
    else {
      get_mapped_reads(true, read_len, adaptor_two, read_index_two[two_idx],
		       bests_two[two_idx].mr, in_two, line_count_two, two);
      process_single_read(two, out);
      get_next(bests_two, two_idx);
    }
  }
  
  while (one_idx < bests_one.size()) {
    get_mapped_reads(false, read_len, adaptor_one, read_index_one[one_idx],
		     bests_one[one_idx].mr, in_one, line_count_one, one);
    process_single_read(one, out);
    get_next(bests_one, one_idx);
  }
  
  while (two_idx < bests_two.size()) {
    get_mapped_reads(true, read_len, adaptor_two, read_index_two[two_idx],
		     bests_two[two_idx].mr, in_two, line_count_two, two);
    process_single_read(two, out);
    get_next(bests_two, two_idx);
  }
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

static void
extract_adaptors(const string &adaptor, string & T_adaptor, string & A_adaptor) {
  const size_t sep_idx = adaptor.find_first_of(":");
  if (adaptor.find_last_of(":") != sep_idx)
    throw SMITHLABException("ERROR: adaptor format \"T_adaptor[:A_adaptor]\"");
  if (sep_idx == string::npos)
    T_adaptor = A_adaptor = adaptor;
  else {
    T_adaptor = adaptor.substr(0, sep_idx);
    A_adaptor = adaptor.substr(sep_idx + 1);
  }
}



int 
main(int argc, const char **argv) {
  try {
    
    string chrom_file;
    string outfile;
    string ambiguous_file;
    string fasta_suffix = "fa";
    string adaptor_sequence;
    
    size_t max_mismatches = numeric_limits<size_t>::max();
    size_t max_mappings = 100;
    size_t read_start_index = 0;
    size_t n_reads_to_process = numeric_limits<size_t>::max();
    
    bool VERBOSE = false;
    size_t range = 1000;
    size_t suffix_len = 1;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "map Illumina reads from BS-seq",
			   "<fastq-reads-file>");
    opt_parse.add_opt("output", 'o', "output file name", 
		      true , outfile);
    opt_parse.add_opt("chrom", 'c', "chromosomes in FASTA file or dir", 
		      true , chrom_file);
    opt_parse.add_opt("start", 'T', "index of first read to map", 
		      false , read_start_index);
    opt_parse.add_opt("number", 'N', "number of reads to map", 
		      false , n_reads_to_process);
    opt_parse.add_opt("suffix", 's', "suffix of chrom files "
		      "(assumes dir provided)", false , fasta_suffix);
    opt_parse.add_opt("mismatch", 'm', "maximum allowed mismatches", 
		      false , max_mismatches);
    opt_parse.add_opt("max-map", 'M', "maximum allowed mappings for a read", 
		      false, max_mappings);
    opt_parse.add_opt("clip", 'C', "clip the specified adaptor",
		      false, adaptor_sequence);
    opt_parse.add_opt("fraglen", 'L', "max fragment length", false, range);
    opt_parse.add_opt("suffix-len", '\0', "Suffix length of reads name", 
		      false, suffix_len);
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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string reads_file_one(leftover_args.front());
    const string reads_file_two(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/
    
    string T_adaptor, A_adaptor;
    extract_adaptors(adaptor_sequence, T_adaptor, A_adaptor);
    
    //////////////////////////////////////////////////////////////
    //  DETERMINE WHICH CHROMOSOMES WILL USED IN MAPPING
    //
    vector<string> chrom_files;
    identify_chromosomes(VERBOSE, fasta_suffix, chrom_file, chrom_files);
    
    //////////////////////////////////////////////////////////////
    // INITIALIZE THE STRUCTURES THAT HOLD THE RESULTS
    //
    MultiMapResult::init(max_mappings);
    
    //////////////////////////////////////////////////////////////
    // INITIALIZE THE SEED STRUCTURES
    //
    vector<size_t> the_seeds;
    load_seeds(VERBOSE, the_seeds);

    if (VERBOSE)
      cerr << "[STARTING END ONE]" << endl;
    bool AG_WILDCARD = false;
    vector<unsigned int> read_index_one;
    vector<MultiMapResult> best_maps_one;
    size_t read_width = 0;
    map_reads(VERBOSE, false, T_adaptor, chrom_files,
    	      reads_file_one, read_start_index, n_reads_to_process, the_seeds, 
    	      max_mismatches, read_index_one, best_maps_one, read_width);

    if (VERBOSE)
      cerr << endl << "[STARTING END TWO]" << endl;
    AG_WILDCARD = true;
    vector<unsigned int> read_index_two;
    vector<MultiMapResult> best_maps_two;
    map_reads(VERBOSE, true, A_adaptor, chrom_files,
	      reads_file_two, read_start_index, n_reads_to_process, the_seeds, 
	      max_mismatches, read_index_two, best_maps_two, read_width);
    
    clip_and_output(VERBOSE, range, read_width, suffix_len, reads_file_one, 
		    reads_file_two, T_adaptor, A_adaptor, read_index_one, 
		    read_index_two, best_maps_one, best_maps_two, outfile);
    
  }
  catch (const SMITHLABException &e) {
    cerr << endl << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
