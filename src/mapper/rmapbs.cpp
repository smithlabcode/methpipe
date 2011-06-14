/*    rmapbs: a program for mapping bisulfite treated Solexa reads
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
#include <cmath>
#include <tr1/unordered_map>

#include "FastRead.hpp"
#include "FastReadWC.hpp"
#include "FastReadQuality.hpp"
#include "rmap_os.hpp"
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

enum { FASTA_FILE, FASTQ_FILE, FASTA_AND_PRB };
enum { RUN_MODE_MISMATCH, RUN_MODE_WILDCARD, RUN_MODE_WEIGHT_MATRIX };

static void
bisulfite_treatment(bool AG_WC, size_t &x) 
{
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
template <typename T> struct SeedHash 
{
    typedef
    unordered_map<size_t,
                  pair<typename vector<T>::const_iterator,
                       typename vector<T>::const_iterator> > type;
};
typedef vector<pair<size_t, unsigned int> > SeedHashSorter;


static void
load_seeds(const bool VERBOSE, const bool FASTER_MODE,
           const size_t read_width, const size_t n_seeds, 
           const size_t seed_weight, vector<size_t> &the_seeds) 
{
    if (FASTER_MODE) 
    {
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
    }
    else SeedMaker::first_last_seeds(min(read_width, SeedMaker::max_seed_part),
                                     n_seeds, seed_weight, the_seeds);
    if (VERBOSE) 
    {
        cerr << endl << "SEED STRUCTURES:" << endl;
        for (size_t i = 0; i < the_seeds.size(); ++i)
            cerr << bits2string_masked(rmap_bits::all_ones, the_seeds[i]) << endl;
        cerr << endl;
    }
}

inline size_t
bisulfite_treatment_tc(char c) 
{
    switch(c) 
    {
    case 1 : return 3;
    default : return c;
    }
}

struct regular_score 
{
    regular_score() 
    {}
    template <class T, class U> size_t score_tc(T a, U b) const 
    {
        return a->score_tc(b);
    }
    template <class T, class U> size_t score_ag(T a, U b) const 
    {
        return a->score_ag(b);
    }
};

struct wildcard_score 
{
    wildcard_score() {}
    template <class T, class U> size_t score_tc(T a, U b) const 
    {
        return a->score_tc_wild_n(b);
    }
    template <class T, class U> size_t score_ag(T a, U b) const 
    {
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
             const bool strand, vector<MultiMapResult> &best_maps) 
{
  
    MASK_t bad_bases = rmap_bits::all_ones;
    MASK_t read_word = rmap_bits::all_zeros;
    T fast_read;
  
    const size_t key_diff = read_width - min(read_width, SeedMaker::max_seed_part);
    const size_t chrom_size = chrom.size();
  
    size_t chrom_offset = 0;
    while(chrom_offset < min(chrom_size, key_diff))
        fast_read.shift(chrom[chrom_offset++]);
  
    while (chrom_offset < min(chrom_size, read_width - 1)) 
    {
        const size_t key_base = chrom[chrom_offset - key_diff];
        fast_read.shift(chrom[chrom_offset++]);
        SeedMaker::update_bad_bases(key_base, bad_bases);
        SeedMaker::update_read_word(key_base, read_word);
    }
  
    string::const_iterator key_pos(chrom.begin() + (chrom_offset - key_diff));
    const string::const_iterator chrom_lim(chrom.end());
    for (string::const_iterator chrom_pos(chrom.begin() + chrom_offset); 
         chrom_pos < chrom_lim; ++chrom_pos, ++key_pos) 
    {
        const size_t key_base = bisulfite_treatment_tc(*key_pos);
        fast_read.shift(*chrom_pos);
        SeedMaker::update_bad_bases(key_base, bad_bases);
        SeedMaker::update_read_word(key_base, read_word);
        if ((bad_bases & profile) == 0) 
        {
            typename SeedHash<T>::type::const_iterator 
                bucket(seed_hash.find(read_word & profile));
            if (bucket != seed_hash.end()) 
            {
                pair<typename vector<T>::const_iterator, 
                    typename vector<T>::const_iterator> tmp(bucket->second);
                const typename vector<T>::const_iterator limit(tmp.second);
                for (typename vector<T>::const_iterator to_test(tmp.first); 
                     to_test != limit; ++to_test) 
                {
                    // const size_t score = to_test->score_tc(fast_read);
                    const size_t score = specialized_score.score_tc(to_test, fast_read);
                    if (score <= max_diffs) 
                    {
                        const vector<MultiMapResult>::iterator 
                            current(best_maps.begin() + (to_test - fast_reads.begin()));
                        if (score <= current->score)
                            current->add(score, chrom_id, distance(chrom.begin(), chrom_pos)
                                         - read_width + 1, strand);
                    }
                }
        }
    }
}
    }


inline size_t
bisulfite_treatment_ag(char c) 
{
    switch(c) 
    {
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
             const bool strand, vector<MultiMapResult> &best_maps) 
{
  
    MASK_t bad_bases = rmap_bits::all_ones;
    MASK_t read_word = rmap_bits::all_zeros;
    T fast_read;
  
    const size_t key_diff = read_width - min(read_width, SeedMaker::max_seed_part);
    const size_t chrom_size = chrom.size();
  
    size_t chrom_offset = 0;
    while(chrom_offset < min(chrom_size, key_diff))
        fast_read.shift(chrom[chrom_offset++]);

    while (chrom_offset < min(chrom_size, read_width - 1)) 
    {
        const size_t key_base = chrom[chrom_offset - key_diff];
        fast_read.shift(chrom[chrom_offset++]);
        SeedMaker::update_bad_bases(key_base, bad_bases);
        SeedMaker::update_read_word(key_base, read_word);
    }
  
    string::const_iterator key_pos(chrom.begin() + (chrom_offset - key_diff));
    const string::const_iterator chrom_lim(chrom.end());
    for (string::const_iterator chrom_pos(chrom.begin() + chrom_offset); 
         chrom_pos < chrom_lim; ++chrom_pos, ++key_pos) 
    {
        const size_t key_base = bisulfite_treatment_ag(*key_pos);
        fast_read.shift(*chrom_pos);
        SeedMaker::update_bad_bases(key_base, bad_bases);
        SeedMaker::update_read_word(key_base, read_word);
        if ((bad_bases & profile) == 0) 
        {
            typename SeedHash<T>::type::const_iterator 
                bucket(seed_hash.find(read_word & profile));
            if (bucket != seed_hash.end()) 
            {
                pair<typename vector<T>::const_iterator, 
                    typename vector<T>::const_iterator> tmp(bucket->second);
                const typename vector<T>::const_iterator limit(tmp.second);
                for (typename vector<T>::const_iterator to_test(tmp.first); 
                     to_test != limit; ++to_test) 
                {
                    // const size_t score = to_test->score_ag(fast_read);
                    const size_t score = specialized_score.score_ag(to_test, fast_read);
                    if (score <= max_diffs) 
                    {
                        const vector<MultiMapResult>::iterator 
                            current(best_maps.begin() + distance(fast_reads.begin(), to_test));
                        if (score <= current->score)
                            current->add(score, chrom_id, distance(chrom.begin(), chrom_pos)
                                         - read_width + 1, strand);
                    }
                }
        }
    }
}
    }


static void
treat_cpgs(const bool AG_WILDCARD, string &chrom) 
{
    const size_t lim = chrom.length() - (!AG_WILDCARD);
    if (AG_WILDCARD) 
    {
        for (size_t i = 1; i < lim; ++i)
            if (chrom[i] == 0 and chrom[i - 1] == 1) chrom[i] = 2;
    }
    else for (size_t i = 0; i < lim; ++i)
             if (chrom[i] == 3 and chrom[i + 1] == 2) chrom[i] = 1;
}


static void
get_read_matches(const size_t the_seed, const vector<size_t> &read_words,
                 SeedHashSorter &sh_sorter) 
{
    const size_t lim = read_words.size();
    sh_sorter.resize(read_words.size());
    for (size_t i = 0; i < lim; ++i)
        sh_sorter[i] = make_pair(the_seed & read_words[i], i);
    sort(sh_sorter.begin(), sh_sorter.end());
}


template <class T> void
sort_by_key(const SeedHashSorter &sh, vector<T> &in) 
{
    vector<T> tmp(in);
    size_t j = 0;
    for (SeedHashSorter::const_iterator i(sh.begin()); i != sh.end(); ++i, ++j)
        in[j] = tmp[i->second];
}


template <class T> void
sort_by_key(SeedHashSorter &sh_sorter, vector<MultiMapResult> &best_maps,
            vector<size_t> &reads, vector<unsigned int> &read_index, 
            vector<T> &fast_reads) 
{
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
                typename SeedHash<T>::type &seed_hash) 
{
    seed_hash.clear();
    typename vector<T>::const_iterator frb(fast_reads.begin());
    size_t prev_key = 0, prev_idx = 0, curr_idx = 0;
    for (SeedHashSorter::const_iterator shs(sh_sorter.begin()); 
         shs != sh_sorter.end(); ++shs) 
    {
        curr_idx = shs->second;
        if (shs->first != prev_key) 
        {
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
             typename SeedHash<T>::type &seed_hash) 
{
    seed_hash.clear();
    SeedHashSorter sh_sorter;
    get_read_matches(the_seed, read_words, sh_sorter);
    sort_by_key(sh_sorter, best_maps, read_words, read_index, fast_reads);
    build_seed_hash(sh_sorter, fast_reads, seed_hash);
}


unsigned char
b2i(char c) 
{
    switch(c) 
    {
    case 'A' : return 0;
    case 'C' : return 1;
    case 'G' : return 2;
    case 'T' : return 3;
    default : return 4;
    }
}


unsigned char
b2i_rc(char c) 
{
    switch(c) 
    {
    case 'A' : return 3;
    case 'C' : return 2;
    case 'G' : return 1;
    case 'T' : return 0;
    default : return 4;
    }
}


template <class T, class U> void
iterate_over_seeds(const bool VERBOSE, const bool AG_WILDCARD, 
                   const bool ALLOW_METH_BIAS,
                   const U &specialized_score,
                   const vector<size_t> &the_seeds, 
                   const vector<string> &chrom_files,
                   vector<pair<unsigned int, unsigned int> > &ambigs, 
                   vector<string> &chrom_names, 
                   vector<size_t> &chrom_sizes,
                   vector<T> &fast_reads,  vector<size_t> &read_words, 
                   vector<unsigned int> &read_index,
                   vector<MultiMapResult> &best_maps,
                   const size_t max_mismatches, const size_t read_width) 
{
  
    if (VERBOSE)
        cerr << "[SCANNING CHROMOSOMES]" << endl;
  
    for (size_t j = 0; j < the_seeds.size() && !fast_reads.empty(); ++j) 
    {
        if (VERBOSE)
            cerr << "[SEED:" << j + 1 << "/" << the_seeds.size() << "] "
                 << "[FORMATTING READS]" << endl;
    
        typename SeedHash<T>::type seed_hash;
        resort_reads(the_seeds[j], fast_reads, read_words, 
                     read_index, best_maps, seed_hash);
    
        size_t prev_chrom_count = 0;
        for (size_t i = 0; i < chrom_files.size() && !fast_reads.empty(); ++i) 
        {
      
            vector<string> tmp_chrom_names, chroms;
            if (VERBOSE)
                cerr << "[SEED:" << j + 1 << "/" << the_seeds.size() << "] "
                     << "[LOADING CHROM] ";
            read_fasta_file(chrom_files[i].c_str(), tmp_chrom_names, chroms);
      
            if (VERBOSE)
                cerr << "[SCANNING=" << tmp_chrom_names.front() << "] ";
      
            const clock_t start(clock());
            for (size_t k = 0; k < chroms.size(); ++k) 
            {
                if (j == 0) 
                {
                    chrom_sizes.push_back(chroms[k].length());
                    chrom_names.push_back(tmp_chrom_names[k]);
                }
                transform(chroms[k].begin(), chroms[k].end(), chroms[k].begin(), 
                          ptr_fun(&toupper));
                string tmp_chrom(chroms[k]);
                transform(chroms[k].begin(), chroms[k].end(), tmp_chrom.begin(),
                          ptr_fun(&b2i));
                if (!ALLOW_METH_BIAS)
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
                if (!ALLOW_METH_BIAS)
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
                cerr << "[" << static_cast<float>(end - start)/
                    CLOCKS_PER_SEC << " SEC]" << endl;
            if (j == 0 && (i + 1) < chrom_files.size()) 
            {
                if (VERBOSE)
                    cerr << "[CLEANING] ";
                eliminate_ambigs(0, the_seeds[j], best_maps, read_index, 
                                 read_words, ambigs, fast_reads, seed_hash);
                if (VERBOSE)
                    cerr << "[AMBIG=" << ambigs.size() << "] " << endl;
            }
            prev_chrom_count += chroms.size();
        }
        if (j == 0) 
        {
            if (VERBOSE)
                cerr << "[CLEANING] ";
            eliminate_ambigs(1, best_maps, read_index, read_words, ambigs, fast_reads);
            if (VERBOSE)
                cerr << "[AMBIG=" << ambigs.size() << "] " << endl;
        }
    }
    if (VERBOSE)
        cerr << "[FINAL CLEANING] ";
    eliminate_ambigs(max_mismatches, best_maps, read_index, 
                     read_words, ambigs, fast_reads);
    if (VERBOSE)
        cerr << "[AMBIG=" << ambigs.size() << "] " << endl;
    vector<T>().swap(fast_reads);
    vector<size_t>().swap(read_words);
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////
////  CODE TO PREPARE AND WRITE THE OUTPUT
////
////
static void
sites_to_regions(const bool VERBOSE, const size_t RUN_MODE, const size_t read_len,
                 const vector<string> &chrom, const vector<size_t> &chrom_sizes, 
                 const vector<unsigned int> &read_index, 
                 vector<string> &read_names, vector<MultiMapResult> &bests, 
                 const string &outfile) 
{
    if (VERBOSE)
        cerr << "[WRITING MAPS] ";
    std::ostream* out = (!outfile.empty()) ? 
        new std::ofstream(outfile.c_str()) : &cout;
    size_t n_reads_mapped = 0;
    for (size_t i = 0; i < bests.size(); ++i)
        if (!bests[i].empty()) 
        {
            bests[i].sort();
            for (size_t j = 0; j < bests[i].mr.size(); ++j)
                if (j == 0 || bests[i].mr[j - 1] < bests[i].mr[j]) 
                {
                    const size_t chrom_id = bests[i].mr[j].chrom;
                    size_t start = bests[i].mr[j].strand ? 
                        bests[i].mr[j].site : 
                        chrom_sizes[chrom_id] - bests[i].mr[j].site - read_len;
                    size_t end = start + read_len;
                    const double score = (RUN_MODE == RUN_MODE_WEIGHT_MATRIX) ?
                        FastReadQuality::value_to_quality(bests[i].score) :
                        bests[i].score;
                    const char strand = ((bests[i].mr[j].strand) ? '+' : '-');
                    size_t real_read_len = (read_names[read_index[i]]).length();
                    if(strand == '+')
			end = start + real_read_len;
		    else
			start = end - real_read_len;
		    *out << GenomicRegion(chrom[chrom_id], start, end,
                                          read_names[read_index[i]], 
                                          score, strand) << '\n';
                    n_reads_mapped++;
                }
            string().swap(read_names[read_index[i]]);
            bests[i].clear();
        }
    if (out != &cout) delete out;
    if (VERBOSE)
        cerr << "[DONE] " << endl
             << "TOTAL READS MAPPED: " << n_reads_mapped << endl;
}


static void
write_non_uniques(string filename, 
                  const vector<pair<unsigned int, unsigned int> > &ambigs,
                  const vector<string> &read_names) 
{
    std::ofstream out(filename.c_str());
    for (size_t i = 0; i < ambigs.size(); ++i)
        out << read_names[ambigs[i].first] << "\t" << ambigs[i].second << endl;
    out.close();
}


template <class T> void
eliminate_ambigs(const size_t max_mismatches, const size_t the_seed,
                 vector<MultiMapResult> &best_maps, 
                 vector<unsigned int> &read_index, vector<size_t> &read_words, 
                 vector<pair<unsigned int, unsigned int> > &ambigs, vector<T> &fast_reads,
                 typename SeedHash<T>::type &seed_hash) 
{
    size_t prev_idx = 0, j = 0;
    size_t prev_key = 0;
    seed_hash.clear();
    typename vector<T>::const_iterator frb(fast_reads.begin());
    for (size_t i = 0; i < best_maps.size(); ++i) 
    {
        best_maps[i].collapse();
        if (best_maps[i].ambiguous() && best_maps[i].score <= max_mismatches)
            ambigs.push_back(make_pair(read_index[i], best_maps[i].score));
        else 
        {
            best_maps[j] = best_maps[i];
            read_index[j] = read_index[i];
            fast_reads[j] = fast_reads[i];
            read_words[j] = read_words[i];
            const size_t key = (read_words[j] & the_seed);
            if (j > 0 && key != prev_key) 
            {
                seed_hash[prev_key] = make_pair(frb + prev_idx, frb + j);
                prev_idx = j;
            }
            ++j;
            prev_key = key;
        }
    }
    seed_hash[prev_key] = make_pair(frb + prev_idx, frb + j);
    best_maps.erase(best_maps.begin() + j, best_maps.end());
    // This below should work but doesn't... Is there a bug elsewhere?
    // vector<MultiMapResult>(best_maps).swap(best_maps);
    read_index.erase(read_index.begin() + j, read_index.end());
    vector<unsigned int>(read_index).swap(read_index);
    read_words.erase(read_words.begin() + j, read_words.end());
    vector<size_t>(read_words).swap(read_words);
    fast_reads.erase(fast_reads.begin() + j, fast_reads.end());
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


static void
identify_chromosomes(const bool VERBOSE,
                     const string filenames_file,
                     const string fasta_suffix,
                     const string chrom_file, 
                     vector<string> &chrom_files) 
{
    if (VERBOSE)
        cerr << "[IDENTIFYING CHROMS] ";
    if (!filenames_file.empty())
        read_filename_file(filenames_file.c_str(), chrom_files);
    else if (isdir(chrom_file.c_str())) 
        read_dir(chrom_file, fasta_suffix, chrom_files);
    else chrom_files.push_back(chrom_file);
    if (VERBOSE) 
    {
        cerr << "[DONE]" << endl 
             << "chromosome files found (approx size):" << endl;
        for (vector<string>::const_iterator i = chrom_files.begin();
             i != chrom_files.end(); ++i)
            cerr << *i << " (" << roundf(get_filesize(*i)/1e06) << "Mbp)" << endl;
        cerr << endl;
    }
}



static void
load_read_names(const size_t INPUT_MODE, 
                string filename, vector<string> &names) 
{
    static const size_t INPUT_BUFFER_SIZE = 10000;
    std::ifstream in(filename.c_str(), std::ios::binary);
    if (!in)
        throw RMAPException("cannot open input file " + string(filename));
    size_t line_count = 0;
    while (!in.eof()) 
    {
        char buffer[INPUT_BUFFER_SIZE + 1];
        in.getline(buffer, INPUT_BUFFER_SIZE);
        if (in.gcount() == static_cast<int>(INPUT_BUFFER_SIZE))
            throw RMAPException("Line in " + filename + "\nexceeds max length: " +
                                toa(INPUT_BUFFER_SIZE));
        if ((INPUT_MODE == FASTQ_FILE && line_count % 4 == 0) ||
            (INPUT_MODE != FASTQ_FILE && line_count % 2 == 0)) 
        {
            names.push_back(buffer + 1);
            const size_t name_end = names.back().find_first_of(" \t");
            if (name_end != string::npos)
                names.back().erase(names.back().begin() + name_end,
                                   names.back().end());
        }
        ++line_count;
    }
    in.close();
}



static size_t
get_input_mode(const bool VERBOSE, 
               const string &reads_file, const string &prb_file) 
{
    size_t INPUT_MODE = FASTA_FILE;
    if (is_fastq(reads_file)) INPUT_MODE = FASTQ_FILE;
    if (!prb_file.empty()) INPUT_MODE = FASTA_AND_PRB;
    if (VERBOSE)
        cerr << "INPUT MODE: "
             << ((INPUT_MODE == FASTA_FILE) ? 
                 "FASTA" : ((INPUT_MODE == FASTQ_FILE) ? 
                            "FASTQ" : "FASTA+PRB")) << endl;
    return INPUT_MODE;
}  


static size_t
get_run_mode(const bool VERBOSE, const size_t INPUT_MODE, 
             const bool WILDCARD, const bool QUALITY) 
{
    size_t RUN_MODE = RUN_MODE_MISMATCH;
    if (WILDCARD and QUALITY)
        throw RMAPException("wildcard and quality matching: mutually exclusive");
    if (WILDCARD) 
    {
        if (INPUT_MODE == FASTA_FILE)
            throw RMAPException("quality score information "
                                "required to use wildcards");
        RUN_MODE = RUN_MODE_WILDCARD;
    }
    else if (INPUT_MODE == FASTA_AND_PRB || 
             (INPUT_MODE == FASTQ_FILE && QUALITY))
        RUN_MODE = RUN_MODE_WEIGHT_MATRIX;
    if (VERBOSE)
        cerr << "MATCH MODE: "
             << ((RUN_MODE == RUN_MODE_MISMATCH) ? 
                 "MISMATCH" : ((RUN_MODE == RUN_MODE_WILDCARD) ? 
                               "WILDCARD" : "WEIGHT-MATRIX")) << endl;
    return RUN_MODE;
}  


static void
load_reads(const bool VERBOSE, const size_t INPUT_MODE, 
           const size_t RUN_MODE, const bool AG_WILDCARD,
           const size_t max_mismatches, const string &adaptor,
           const string &reads_file, const string &prb_file,
           vector<FastRead> &fast_reads, vector<FastReadWC> &fast_reads_wc,
           vector<FastReadQuality> &fast_reads_q,
           vector<unsigned int> &read_index, vector<size_t> &read_words,
           size_t &read_width) 
{

    //////////////////////////////////////////////////////////////
    // LOAD THE READS (AS SEQUENCES OR PROBABILITIES) FROM DISK
    if (VERBOSE) cerr << "[LOADING READ SEQUENCES] ";
    vector<string> reads;
    if (INPUT_MODE == FASTQ_FILE) 
    {
        if (RUN_MODE == RUN_MODE_WILDCARD)
            load_reads_from_fastq_file(reads_file, adaptor, max_mismatches, read_width,
                                       fast_reads_wc, read_words, read_index);
        else if (RUN_MODE == RUN_MODE_WEIGHT_MATRIX)
            load_reads_from_fastq_file(reads_file, adaptor, max_mismatches, read_width,
                                       fast_reads_q, read_words, read_index);
        else
            load_reads_from_fastq_file(reads_file, adaptor, max_mismatches, read_width,
                                       fast_reads, read_words, read_index);
    }
    else if (INPUT_MODE == FASTA_AND_PRB) 
    {
        if (RUN_MODE == RUN_MODE_WILDCARD)
            load_reads_from_prb_file(prb_file, adaptor, max_mismatches, read_width,
                                     fast_reads_wc, read_words, read_index);
        else load_reads_from_prb_file(prb_file, adaptor, max_mismatches, read_width,
                                      fast_reads_q, read_words, read_index);
    }
    else load_reads_from_fasta_file(reads_file, adaptor, max_mismatches, read_width,
                                    fast_reads, read_words, read_index);
    for (size_t i = 0; i < fast_reads_wc.size(); ++i)
        fast_reads_wc[i].bisulfite_treatment(AG_WILDCARD);
    for (size_t i = 0; i < fast_reads_q.size(); ++i)
        fast_reads_q[i].bisulfite_treatment(AG_WILDCARD);
    for (size_t i = 0; i < read_words.size(); ++i)
        bisulfite_treatment(AG_WILDCARD, read_words[i]);
    if (VERBOSE)
        cerr << "[DONE]" << endl
             << "TOTAL HQ READS: " << read_index.size() << endl
             << "READ WIDTH: " << read_width << endl;
}


struct indexed_best_less 
{
    bool operator()(const pair<unsigned int, MultiMapResult> &a,
                    const pair<unsigned int, MultiMapResult> &b) const 
    {
        return a.first < b.first;
    }
};

static void
invert_bests_list(vector<unsigned int> &read_index, 
                  vector<MultiMapResult> &bests) 
{
    vector<pair<unsigned int, MultiMapResult> > sorter;
    for (size_t i = 0; i < bests.size(); ++i)
        sorter.push_back(make_pair(read_index[i], bests[i]));
    sort(sorter.begin(), sorter.end(), indexed_best_less());
    for (size_t i = 0; i < sorter.size(); ++i) 
    {
        read_index[i] = sorter[i].first;
        bests[i] = sorter[i].second;
    }
}


static void
iterate_over_reads(const bool VERBOSE,
                   const size_t RUN_MODE, 
                   const string &adaptor,
                   const string &filename,
                   const string &outfile,
                   const size_t read_len,
                   const vector<unsigned int> &reads_index, 
                   vector<MultiMapResult> &bests, 
                   const vector<size_t> &chrom_sizes,
                   const vector<string> &chrom) 
{
  
    static const size_t INPUT_BUFFER_SIZE = 10000;
    std::ifstream in(filename.c_str(), std::ios::binary);
    if (!in)
        throw RMAPException("cannot open input file " + string(filename));
    std::ostream* out = (!outfile.empty()) ? 
        new std::ofstream(outfile.c_str()) : &cout;

    if (!outfile.empty() && !(*out))
        throw RMAPException("cannot open output file " + string(outfile));
  
    size_t read_idx = 0, line_count = 0, 
        curr_idx = 0, n_reads_mapped = 0;
    string name, sequence;
  
    while (!in.eof()) 
    {

        char buffer[INPUT_BUFFER_SIZE + 1];
        in.getline(buffer, INPUT_BUFFER_SIZE);
        if (in.gcount() == static_cast<int>(INPUT_BUFFER_SIZE))
            throw RMAPException("Line in " + filename + "\nexceeds max length: " +
                                toa(INPUT_BUFFER_SIZE));
    
        if (line_count % 4 == 0 && read_idx == reads_index[curr_idx]) 
        {
            name = string(buffer + 1);
            const size_t name_end = name.find_first_of(" \t");
            if (name_end != string::npos)
                name.erase(name.begin() + name_end, name.end());
        }
        else if (line_count % 4 == 1 && read_idx == reads_index[curr_idx]) 
        {
            sequence = string(buffer);
        }
        else if (line_count % 4 == 3 && read_idx == reads_index[curr_idx]) 
        {
            if (!bests[curr_idx].empty()) 
            {
                if (!adaptor.empty())
                    clip_adaptor_from_read(adaptor, MIN_ADAPTOR_MATCH_SCORE, sequence);
                bests[curr_idx].sort();
                for (size_t j = 0; j < bests[curr_idx].mr.size(); ++j)
                    if (j == 0 || bests[curr_idx].mr[j - 1] < bests[curr_idx].mr[j]) 
                    {
                        const size_t chrom_id = bests[curr_idx].mr[j].chrom;
                        size_t start = bests[curr_idx].mr[j].strand ? 
                            bests[curr_idx].mr[j].site : 
                            chrom_sizes[chrom_id] - bests[curr_idx].mr[j].site - read_len;
                        size_t end = start + read_len;
                        const double score = 
                            (RUN_MODE == RUN_MODE_WEIGHT_MATRIX) ?
                            FastReadQuality::value_to_quality(bests[curr_idx].score) :
                            bests[curr_idx].score;
                        const char strand = ((bests[curr_idx].mr[j].strand) ? '+' : '-');
			size_t real_read_len = sequence.length();
			if(strand == '+')
			    end = start + real_read_len;
			else
			    start = end - real_read_len;
                        *out << GenomicRegion(chrom[chrom_id], start, end,
                                              name, score, strand) << '\t'
                             << sequence << "\t" << string(buffer) << endl;
                        n_reads_mapped++;
                    }
                bests[curr_idx].clear();
            }
            ++curr_idx;
        }
        if (line_count % 4 == 3) 
            ++read_idx;
        ++line_count;
    }
    in.close();
    if (out != &cout) delete out;
    if (VERBOSE)
        cerr << "[DONE] " << endl
             << "TOTAL READS MAPPED: " << n_reads_mapped << endl;
}


int 
main(int argc, const char **argv) 
{
    try 
    {
    
        string chrom_file;
        string filenames_file;
        string outfile;
        string prb_file;
        string ambiguous_file;
        string fasta_suffix = "fa";
        string adaptor_sequence;
    
        size_t n_seeds = 3;
        size_t seed_weight = 11;
        size_t read_width = 0;
        size_t max_mismatches = 10;
        size_t max_mappings = 1;
        double wildcard_cutoff = numeric_limits<double>::max();
    
        bool VERBOSE = false;
        bool FASTER_MODE = false;
        bool QUALITY = false;
        bool AG_WILDCARD = false;
        bool ALLOW_METH_BIAS = false;
        bool WILDCARD = false;

        bool WILD_N_MODE = false;
        bool ORIGINAL_OUTPUT = false;
    
        /****************** COMMAND LINE OPTIONS ********************/
        OptionParser opt_parse("rmapbs", "The rmapbs mapping tool for Solexa reads"
                               " following bisulfite treatment",
                               "<fast[a/q]-reads-file>");
        opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
                          false , outfile);
        opt_parse.add_opt("chrom", 'c', "FASTA file or dir containing chromosome(s)", 
                          false , chrom_file);
        opt_parse.add_opt("suffix", 's', "suffix of FASTA files "
                          "(assumes -c indicates dir)", false , fasta_suffix);
        opt_parse.add_opt("filenames", 'F', "file listing names of "
                          "chromosome files", false , filenames_file);
        opt_parse.add_opt("prb", 'p', "file with quality scores (prb format)", 
                          false, prb_file);
        opt_parse.add_opt("seeds", 'S', "number of seeds", false , n_seeds);
        opt_parse.add_opt("hit", 'h', "weight of hit", false , seed_weight);
        opt_parse.add_opt("width", 'w', "width of reads", false, read_width);
        opt_parse.add_opt("mismatch", 'm', "maximum allowed mismatches", 
                          false , max_mismatches);
        opt_parse.add_opt("ambiguous", 'a', "file to write names of ambiguously "
                          "mapped reads", false , ambiguous_file);
        opt_parse.add_opt("max-map", 'M', "maximum allowed mappings for a read", 
                          false, max_mappings);
        opt_parse.add_opt("wc", 'W', "run in wildcard matching mode", 
                          false, WILDCARD);
        opt_parse.add_opt("prob", 'P', "wildcard cutoff probability", 
                          false, wildcard_cutoff);
        opt_parse.add_opt("qual", 'Q', "use quality scores (input must be FASTQ)", 
                          false, QUALITY);
        opt_parse.add_opt("ag-wild", 'A', "map using A/G bisulfite wildcards", 
                          false, AG_WILDCARD);
        opt_parse.add_opt("bias", 'B', "allow CpG non-conversion to assist", 
                          false, ALLOW_METH_BIAS);
        opt_parse.add_opt("faster", 'f', "faster seeds (sensitive to 2 mismatches)", 
                          false, FASTER_MODE);
        opt_parse.add_opt("clip", 'C', "clip the specified adaptor", 
                          false, adaptor_sequence);
        opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
        opt_parse.add_opt("original-output", 'G', "use original output format", 
                          false, ORIGINAL_OUTPUT);
        opt_parse.add_opt("wild-n-mode", 'y', "allow N in read to match anything", 
                          false, WILD_N_MODE);
        vector<string> leftover_args;
        opt_parse.parse(argc, argv, leftover_args);
        if (argc == 1 || opt_parse.help_requested()) 
        {
            cerr << opt_parse.help_message() << endl;
            return EXIT_SUCCESS;
        }
        if (opt_parse.about_requested()) 
        {
            cerr << opt_parse.about_message() << endl;
            return EXIT_SUCCESS;
        }
        if (opt_parse.option_missing()) 
        {
            cerr << opt_parse.option_missing_message() << endl;
            return EXIT_SUCCESS;
        }
        if (leftover_args.empty()) 
        {
            cerr << opt_parse.help_message() << endl;
            return EXIT_SUCCESS;
        }
        if (chrom_file.empty() && filenames_file.empty()) 
        {
            cerr << "must specify chroms file/dir or filenames file" << endl;
            return EXIT_FAILURE;
        }
        const string reads_file = leftover_args.front();
        /****************** END COMMAND LINE OPTIONS *****************/
    
        FastReadWC::set_cutoff(wildcard_cutoff);
        if (WILDCARD && wildcard_cutoff != numeric_limits<double>::max() &&
            (wildcard_cutoff > 1.0 || wildcard_cutoff < 0)) 
            throw RMAPException("wildcard cutoff must be in [0, 1]");
    
        //////////////////////////////////////////////////////////////
        //  CHECK HOW QUALITY SCORES ARE USED
        //
        const size_t INPUT_MODE = get_input_mode(VERBOSE, reads_file, prb_file);
        if (INPUT_MODE == FASTA_FILE && !ORIGINAL_OUTPUT)
            throw RMAPException("when reads are given as FASTA, "
                                "original output format is required");
    
        const size_t RUN_MODE = get_run_mode(VERBOSE, INPUT_MODE, WILDCARD, QUALITY);
    
        //////////////////////////////////////////////////////////////
        //  DETERMINE WHICH CHROMOSOMES WILL USED IN MAPPING
        //
        vector<string> chrom_files;
        identify_chromosomes(VERBOSE, filenames_file, fasta_suffix, chrom_file, chrom_files);
    
        //////////////////////////////////////////////////////////////
        // OBTAIN THE READS
        // 
        vector<FastRead> fast_reads;
        vector<FastReadWC> fast_reads_wc;
        vector<FastReadQuality> fast_reads_q;
        vector<unsigned int> read_index;
        vector<size_t> read_words;
        load_reads(VERBOSE, INPUT_MODE, RUN_MODE, AG_WILDCARD,
                   max_mismatches, adaptor_sequence, reads_file, prb_file, 
                   fast_reads, fast_reads_wc, fast_reads_q, 
                   read_index, read_words, read_width);
    
        double max_match_score = max_mismatches*FastReadQuality::get_scaler();
    
        //////////////////////////////////////////////////////////////
        // INITIALIZE THE SEED STRUCTURES
        //
        vector<size_t> the_seeds;
        load_seeds(VERBOSE, FASTER_MODE,
                   read_width, n_seeds, seed_weight, the_seeds);
    
        //////////////////////////////////////////////////////////////
        // INITIALIZE THE STRUCTURES THAT HOLD THE RESULTS
        //
        MultiMapResult::init(max_mappings);
        vector<MultiMapResult> 
            best_maps(read_words.size(), MultiMapResult((RUN_MODE == RUN_MODE_WEIGHT_MATRIX) ?
                                                        max_match_score : max_mismatches));
    
        //////////////////////////////////////////////////////////////
        // THIS IS WHERE THE ACTUAL MAPPING HAPPENS
        //
        vector<size_t> chrom_sizes;
        vector<string> chrom_names;
        vector<pair<unsigned int, unsigned int> > ambigs;
    
        if (RUN_MODE == RUN_MODE_MISMATCH) 
        {
            if (WILD_N_MODE) 
            {
                const wildcard_score specialized_score;
                iterate_over_seeds(VERBOSE, AG_WILDCARD, ALLOW_METH_BIAS,
                                   specialized_score,
                                   the_seeds, chrom_files, ambigs, 
                                   chrom_names, chrom_sizes,
                                   fast_reads, // USE REGULAR FAST READS
                                   read_words, read_index,
                                   best_maps, max_mismatches, read_width);
            }
            else 
            {
                const regular_score specialized_score;
                iterate_over_seeds(VERBOSE, AG_WILDCARD, ALLOW_METH_BIAS,
                                   specialized_score,
                                   the_seeds, chrom_files, ambigs, 
                                   chrom_names, chrom_sizes,
                                   fast_reads, // USE REGULAR FAST READS
                                   read_words, read_index,
                                   best_maps, max_mismatches, read_width);
            }
        }
        if (RUN_MODE == RUN_MODE_WILDCARD) 
        {
            const regular_score specialized_score;
            iterate_over_seeds(VERBOSE, AG_WILDCARD, ALLOW_METH_BIAS,
                               specialized_score,
                               the_seeds, chrom_files, ambigs, 
                               chrom_names, chrom_sizes,
                               fast_reads_wc, // USE FAST READS FOR WILDCARD MATCHING
                               read_words, read_index,
                               best_maps, max_mismatches, read_width);
        }
        if (RUN_MODE == RUN_MODE_WEIGHT_MATRIX) 
        {
            const regular_score specialized_score;
            iterate_over_seeds(VERBOSE, AG_WILDCARD, ALLOW_METH_BIAS,
                               specialized_score,
                               the_seeds, chrom_files, ambigs, 
                               chrom_names, chrom_sizes,
                               fast_reads_q, // USE FAST READS FOR QUALITY MATCHING
                               read_words, read_index,
                               best_maps, max_match_score, read_width);
        }

        // First make sure the chrom names don't have spaces (cause
        // problems for later processing)
        for (size_t i = 0; i < chrom_names.size(); ++i) 
        {
            const size_t chr_name_end = chrom_names[i].find_first_of(" \t");
            if (chr_name_end != string::npos)
                chrom_names[i].erase(chr_name_end);
        }
    
        if (!ORIGINAL_OUTPUT) 
        {
            invert_bests_list(read_index, best_maps);
            iterate_over_reads(VERBOSE, RUN_MODE, adaptor_sequence, reads_file, 
                               outfile, read_width, read_index, best_maps, 
                               chrom_sizes, chrom_names);
            //////////////////////////////////////////////////////////////
            // IF IDENTITIES OF AMBIGUOUS READS ARE DESIRED, WRITE THEM
            if (!ambiguous_file.empty()) 
            {
                //////////////////////////////////////////////////////////////
                // LOAD THE NAMES OF READS AGAIN (THEY WILL BE NEEDED)
                vector<string> read_names;
                load_read_names(INPUT_MODE, reads_file, read_names);
                if (VERBOSE)
                    cerr << "[WRITING AMBIGS] ";
                write_non_uniques(ambiguous_file, ambigs, read_names);
                if (VERBOSE)
                    cerr << "[DONE]" << endl
                         << "TOTAL AMBIGS: " << ambigs.size() << endl;
                ambigs.clear();
            }    
        }
        else 
        {
            //////////////////////////////////////////////////////////////
            // LOAD THE NAMES OF READS AGAIN (THEY WILL BE NEEDED)
            //
            vector<string> read_names;
            load_read_names(INPUT_MODE, reads_file, read_names);
      
            //////////////////////////////////////////////////////////////
            // IF IDENTITIES OF AMBIGUOUS READS ARE DESIRED, WRITE THEM
            if (!ambiguous_file.empty()) 
            {
                if (VERBOSE)
                    cerr << "[WRITING AMBIGS] ";
                write_non_uniques(ambiguous_file, ambigs, read_names);
                if (VERBOSE)
                    cerr << "[DONE]" << endl
                         << "TOTAL AMBIGS: " << ambigs.size() << endl;
                ambigs.clear();
            }    
            //////////////////////////////////////////////////////////////
            // TRANSFORM THE RESULT STRUCTURES INTO BED FORMAT FOR OUTPUT
            sites_to_regions(VERBOSE, RUN_MODE, read_width, chrom_names, chrom_sizes, 
                             read_index, read_names, best_maps, outfile);
        }
    }
    catch (const RMAPException &e) 
    {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
    catch (std::bad_alloc &ba) 
    {
        cerr << "ERROR: could not allocate memory" << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
