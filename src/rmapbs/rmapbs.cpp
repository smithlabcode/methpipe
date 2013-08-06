/*    rmapbs: a program for mapping bisulfite treated Solexa reads
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

#include "FastRead.hpp"
#include "smithlab_os.hpp"
#include "SeedMaker.hpp"
#include "MapResult.hpp"
#include "OptionParser.hpp"
#include "load_reads.hpp"
#include "clip_adaptor_from_reads.hpp"

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ostream_iterator;
using std::max;
using std::min;
using std::numeric_limits;
using std::pair;
using std::make_pair;
using std::ptr_fun;


static string
format_time(const clock_t start, const clock_t end) {
  return string("[") + toa(static_cast<float>(end - start)/
			   CLOCKS_PER_SEC) + string(" SEC]");
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
////
template <typename T> struct SeedHash {
  typedef
  unordered_map<size_t,
		pair<typename vector<T>::const_iterator,
		     typename vector<T>::const_iterator> > type;
};
typedef vector<pair<size_t, unsigned int> > SeedHashSorter;


static void
load_seeds(const bool VERBOSE, 
           const size_t read_width, vector<size_t> &the_seeds) {
  // 0b0000111111001100001111110011000011111100110000111111001111111111
  the_seeds.push_back(1138354285449573375ul);
  // 0b1111110011000011111100110000111111001100001111110011000000111100
  the_seeds.push_back(18213668567193169980ul);
  // 0b0011111100110000111111001100001111110011000011111100110000001111
  the_seeds.push_back(4553417141798292495ul);
  // 0b1100001111110011000011111100110000111111001100001111110000110011
  the_seeds.push_back(14119646626644556851ul);
  // 0b0011000011111100110000111111001100001111110011000011111111110000
  the_seeds.push_back(3529911656661139440ul);
  // 0b1100110000111111001100001111110011000011111100110000111111001100
  the_seeds.push_back(14717535969447448524ul);
  // 0b1111001100001111110011000011111100110000111111001100001111000011
  the_seeds.push_back(17514442047644025795ul);
  if (VERBOSE) {
    cerr << endl << "SEED STRUCTURES:" << endl;
    for (size_t i = 0; i < the_seeds.size(); ++i)
      cerr << bits2string_masked(rmap_bits::all_ones, the_seeds[i]) << endl;
    cerr << endl;
  }
}

inline size_t
bisulfite_treatment_tc(char c) {
  switch(c) {
  case 1 : return 3;
  default : return c;
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


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////
////  WHERE THE ACTUAL MAPPING HAPPENS
////
////
template <class T, class U> void
map_reads_tc(const U &specialized_score,
             const string &chrom, const size_t chrom_id,
             const size_t profile, const size_t read_width, 
             const size_t max_diffs, const vector<T> &fast_reads, 
             //const
             typename SeedHash<T>::type &seed_hash, 
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
             //const
             typename SeedHash<T>::type &seed_hash, 
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
sort_by_key(SeedHashSorter &sh_sorter, vector<MultiMapResult> &best_maps,
            vector<size_t> &reads, vector<unsigned int> &read_index, 
            vector<T> &fast_reads) {
  sort_by_key(sh_sorter, best_maps);
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


template <class T> void
eliminate_ambigs(const size_t max_mismatches, vector<MultiMapResult> &best_maps, 
                 vector<unsigned int> &read_index, vector<size_t> &reads, 
                 vector<pair<unsigned int, unsigned int> > &ambigs, vector<T> &fast_reads) 
{
  size_t j = 0;
  for (size_t i = 0; i < best_maps.size(); ++i) 
    {
      best_maps[i].collapse();
      if (best_maps[i].ambiguous() && best_maps[i].score <= max_mismatches)
	ambigs.push_back(make_pair(read_index[i], best_maps[i].score));
      else 
        {
	  best_maps[j] = best_maps[i];
	  read_index[j] = read_index[i];
	  reads[j] = reads[i];
	  fast_reads[j] = fast_reads[i];
	  ++j;
        }
    }
  best_maps.erase(best_maps.begin() + j, best_maps.end());
  vector<MultiMapResult>(best_maps).swap(best_maps);
  read_index.erase(read_index.begin() + j, read_index.end());
  vector<unsigned int>(read_index).swap(read_index);
  reads.erase(reads.begin() + j, reads.end());
  vector<size_t>(reads).swap(reads);
  fast_reads.erase(fast_reads.begin() + j, fast_reads.end());
  vector<T>(fast_reads).swap(fast_reads);
}


template <class T, class U> void
iterate_over_seeds(const bool VERBOSE, const bool AG_WILDCARD, 
                   const U &specialized_score,
                   const vector<size_t> &the_seeds, 
                   const vector<string> &chrom_files,
                   vector<pair<unsigned int, unsigned int> > &ambigs, 
                   vector<string> &chrom_names, 
                   vector<size_t> &chrom_sizes,
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
	transform(chroms[k].begin(), chroms[k].end(), chroms[k].begin(), 
		  ptr_fun(&toupper));
	string tmp_chrom(chroms[k]);
	transform(chroms[k].begin(), chroms[k].end(), tmp_chrom.begin(),
		  ptr_fun(&b2i));
	treat_cpgs(AG_WILDCARD, tmp_chrom);
	if (AG_WILDCARD)
	  map_reads_ag(specialized_score, 
		       tmp_chrom, prev_chrom_count + k, the_seeds[j], read_width, 
		       max_mismatches, fast_reads, seed_hash, true, best_maps);
	else map_reads_tc(specialized_score, 
			  tmp_chrom, prev_chrom_count + k, the_seeds[j], read_width, 
			  max_mismatches, fast_reads, seed_hash, true, best_maps);
	string().swap(tmp_chrom);
	
	std::reverse(chroms[k].begin(), chroms[k].end());
	transform(chroms[k].begin(), chroms[k].end(), chroms[k].begin(),
		  ptr_fun(&b2i_rc));
	treat_cpgs(AG_WILDCARD, chroms[k]);
	if (AG_WILDCARD)
	  map_reads_ag(specialized_score, 
		       chroms[k], prev_chrom_count + k, the_seeds[j], read_width, 
		       max_mismatches, fast_reads, seed_hash, false, best_maps);
	else map_reads_tc(specialized_score, 
			  chroms[k], prev_chrom_count + k, the_seeds[j], read_width,
			  max_mismatches, fast_reads, seed_hash, false, best_maps);
	string().swap(chroms[k]);
      }
      const clock_t end(clock());
      if (VERBOSE)
	cerr << format_time(start, end) << endl;
      prev_chrom_count += chroms.size();
    }
    for (size_t x = 0; x < best_maps.size(); ++x)
      best_maps[x].collapse();
  }
  eliminate_ambigs(max_mismatches, best_maps, read_index, 
		   read_words, ambigs, fast_reads);
  if (VERBOSE)
    cerr << "[AMBIG=" << ambigs.size() << "] " << endl;
  vector<T>().swap(fast_reads);
  vector<size_t>().swap(read_words);
}


static void
write_non_uniques(string filename, 
                  const vector<pair<unsigned int, unsigned int> > &ambigs,
                  const vector<string> &read_names) {
  std::ofstream out(filename.c_str());
  for (size_t i = 0; i < ambigs.size(); ++i)
    out << read_names[ambigs[i].first] << "\t" << ambigs[i].second << endl;
  out.close();
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
load_read_names(string filename, vector<string> &names) {
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in)
    throw SMITHLABException("cannot open input file " + string(filename));
  size_t line_count = 0;
  string buffer;
  while (!in.eof()) {
    getline(in, buffer);
    if (in.good()) {
      if (line_count % 4 == 0) {
	const size_t truncpos = buffer.find_first_of(" \t");
	names.push_back(string(buffer.begin() + 1,
			       (truncpos != string::npos) ?
			       buffer.begin() + truncpos : buffer.end()));
      }
      ++line_count;
    }
  }
}


static void
load_reads(const bool VERBOSE, const bool AG_WILDCARD,
           const size_t max_mismatches, const string &adaptor,
           const string &reads_file, const size_t read_start_index, 
	   const size_t n_reads_to_process, vector<FastRead> &fast_reads,
           vector<unsigned int> &read_index, vector<size_t> &read_words,
           size_t &read_width) {
  
  //////////////////////////////////////////////////////////////
  // LOAD THE READS (AS SEQUENCES OR PROBABILITIES) FROM DISK
  if (VERBOSE) cerr << "[LOADING READ SEQUENCES] ";
  load_reads_from_fastq_file(reads_file, read_start_index, n_reads_to_process,
			     adaptor, read_width, fast_reads, 
			     read_words, read_index);
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
  vector<pair<unsigned int, MultiMapResult> > sorter;
  for (size_t i = 0; i < bests.size(); ++i)
    sorter.push_back(make_pair(read_index[i], bests[i]));
  sort(sorter.begin(), sorter.end(), indexed_best_less());
  for (size_t i = 0; i < sorter.size(); ++i) {
    read_index[i] = sorter[i].first;
    bests[i] = sorter[i].second;
  }
}


static void
iterate_over_reads(const bool VERBOSE, const string &adaptor,
                   const string &filename, const string &outfile,
                   const size_t read_len, const vector<unsigned int> &reads_index, 
                   vector<MultiMapResult> &bests, const vector<size_t> &chrom_sizes,
                   const vector<string> &chrom) {
  
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in)
    throw SMITHLABException("cannot open input file " + string(filename));

  std::ofstream out(outfile.c_str());
  if (!out)
    throw SMITHLABException("cannot open output file " + outfile);
  
  size_t read_idx = 0, line_count = 0, curr_idx = 0, n_reads_mapped = 0;
  string name, sequence;
  
  while (!in.eof() && curr_idx < reads_index.size()) {
    
    string buffer;
    getline(in, buffer);
  
    if (read_idx == reads_index[curr_idx]) {
      if (line_count % 4 == 0) {
	const size_t truncpos = buffer.find_first_of(" \t");
	string(buffer.begin() + 1,
	       (truncpos != string::npos) ?
	       buffer.begin() + truncpos : buffer.end()).swap(name);
      }
      else if (line_count % 4 == 1) {
	sequence.swap(buffer);
      }
      else if (line_count % 4 == 3) {
	if (!bests[curr_idx].empty()) {
	  if (!adaptor.empty())
	    clip_adaptor_from_read(adaptor, sequence);
	  
	  bests[curr_idx].sort();
	  for (size_t j = 0; j < bests[curr_idx].mr.size(); ++j) {
	    const size_t chrom_id = bests[curr_idx].mr[j].chrom;
	    size_t start = bests[curr_idx].mr[j].strand ? 
	      bests[curr_idx].mr[j].site : 
	      chrom_sizes[chrom_id] - bests[curr_idx].mr[j].site - read_len;
	    size_t end = start + read_len;
	    const double score = bests[curr_idx].score;
	    const char strand = ((bests[curr_idx].mr[j].strand) ? '+' : '-');
	    size_t real_read_len = sequence.length();
	    if (strand == '+')
	      end = start + real_read_len;
	    else
	      start = end - real_read_len;
	    out << GenomicRegion(chrom[chrom_id], start, end, name, score, strand) 
		<< '\t' << sequence << "\t" << buffer << endl;
	    n_reads_mapped++;
	  }
	  bests[curr_idx].clear();
	}
	++curr_idx;
      }
    }
    if (line_count % 4 == 3)
      ++read_idx;
    ++line_count;
  }
  if (VERBOSE)
    cerr << "[DONE] " << endl
  	 << "TOTAL READS MAPPED: " << n_reads_mapped << endl;
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
    size_t max_mappings = 1;
    size_t read_start_index = 0;
    size_t n_reads_to_process = std::numeric_limits<size_t>::max();
    
    bool VERBOSE = false;
    bool AG_WILDCARD = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "map Illumina BS-seq reads",
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
    opt_parse.add_opt("ambiguous", 'a', "file to write names of ambiguously "
		      "mapped reads", false , ambiguous_file);
    opt_parse.add_opt("max-map", 'M', "maximum allowed mappings for a read", 
		      false, max_mappings);
    opt_parse.add_opt("ag-wild", 'A', "map using A/G bisulfite wildcards", 
		      false, AG_WILDCARD);
    opt_parse.add_opt("clip", 'C', "clip the specified adaptor",
		      false, adaptor_sequence);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    //////////////////////////////////////////////////////////////
    //  DETERMINE WHICH CHROMOSOMES WILL USED IN MAPPING
    //
    vector<string> chrom_files;
    identify_chromosomes(VERBOSE, fasta_suffix, chrom_file, chrom_files);
    
    //////////////////////////////////////////////////////////////
    // OBTAIN THE READS
    // 
    vector<FastRead> fast_reads;
    vector<unsigned int> read_index;
    vector<size_t> read_words;
    size_t read_width = 0;
    load_reads(VERBOSE, AG_WILDCARD, max_mismatches, adaptor_sequence, 
	       reads_file, read_start_index, n_reads_to_process,
	       fast_reads, read_index, read_words, read_width);
    
    if (max_mismatches == numeric_limits<size_t>::max())
      max_mismatches = static_cast<size_t>(0.07*read_width);
      
    if (VERBOSE)
      cerr << "MAX MISMATCHES: " << max_mismatches << endl;
	
    //////////////////////////////////////////////////////////////
    // INITIALIZE THE SEED STRUCTURES
    //
    vector<size_t> the_seeds;
    load_seeds(VERBOSE, read_width, the_seeds);
    
    //////////////////////////////////////////////////////////////
    // INITIALIZE THE STRUCTURES THAT HOLD THE RESULTS
    //
    MultiMapResult::init(max_mappings);
    vector<MultiMapResult> best_maps(read_words.size());
    
    //////////////////////////////////////////////////////////////
    // THIS IS WHERE THE ACTUAL MAPPING HAPPENS
    //
    vector<size_t> chrom_sizes;
    vector<string> chrom_names;
    vector<pair<unsigned int, unsigned int> > ambigs;
    
    const wildcard_score specialized_score;
    iterate_over_seeds(VERBOSE, AG_WILDCARD, specialized_score,
		       the_seeds, chrom_files, ambigs, chrom_names, 
		       chrom_sizes, fast_reads, read_words, read_index,
		       best_maps, max_mismatches, read_width);
    
    // First make sure the chrom names don't have spaces (cause
    // problems for later processing)
    for (size_t i = 0; i < chrom_names.size(); ++i) {
      const size_t chr_name_end = chrom_names[i].find_first_of(" \t");
      if (chr_name_end != string::npos)
	chrom_names[i].erase(chr_name_end);
    }
    
    invert_bests_list(read_index, best_maps);
    iterate_over_reads(VERBOSE, adaptor_sequence, reads_file, outfile, 
		       read_width, read_index, best_maps, 
		       chrom_sizes, chrom_names);
    
    //////////////////////////////////////////////////////////////
    // IF IDENTITIES OF AMBIGUOUS READS ARE DESIRED, WRITE THEM
    if (!ambiguous_file.empty()) {
      //////////////////////////////////////////////////////////////
      // LOAD THE NAMES OF READS AGAIN (THEY WILL BE NEEDED)
      vector<string> read_names;
      load_read_names(reads_file, read_names);
      if (VERBOSE)
	cerr << "[WRITING AMBIGS] ";
      write_non_uniques(ambiguous_file, ambigs, read_names);
      if (VERBOSE)
	cerr << "[DONE]" << endl
	     << "TOTAL AMBIGS: " << ambigs.size() << endl;
      ambigs.clear();
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
