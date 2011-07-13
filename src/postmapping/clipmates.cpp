/*    checkoverlap: a program for identifying overlapping ends of
 *                  mapped paired end reads.
 *
 *    Copyright (C) 2010 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith, Song Qiang
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

#include <cassert>

#include <string>
#include <vector>
#include <unistd.h>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>

#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "MappedRead.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::ofstream;
using std::ifstream;

static void
check_sorted_by_ID(MappedRead &prev_mr,
                   const MappedRead &mr)
{
    if(prev_mr.r.get_name() > mr.r.get_name())
    {
        cerr << "CLIPMATES ERROR: "
             << "reads are not sorted by reads' ID" << endl
             << "---------------------------------------------" << endl
             << prev_mr << endl
             << mr << endl 
             << "---------------------------------------------" << endl;
        exit(EXIT_FAILURE);
    }
    else
        prev_mr = mr;
}

static void
revcomp(MappedRead &mr)
{
    if (mr.r.get_strand() == '+')
        mr.r.set_strand('-');
    else
        mr.r.set_strand('+');
    
    revcomp_inplace(mr.seq);
    std::reverse(mr.scr.begin(), mr.scr.end());
}


// inline static size_t 
// get_distance(const MappedRead &a, const MappedRead &b) 
// {
//     return (a.r.pos_strand()) ? b.r.get_end() - a.r.get_start() :
//         a.r.get_end() - b.r.get_start();
// }


// static void
// fix_overlap(const MappedRead &a, MappedRead &b) 
// {
//     if (a.r.get_strand() == b.r.get_strand()) 
//     {
//         if (a.r.get_start() < b.r.get_start()) 
//         {
//             const size_t overlap = a.r.get_end() - b.r.get_start();
//             assert(overlap > 0);
//             b.seq = (a.r.pos_strand()) ?
//                 string(overlap, 'N') + b.seq.substr(overlap) :
//                 b.seq.substr(0, b.seq.length() - overlap) + string(overlap, 'N');
//         }
//         else 
//         {
//             const size_t overlap = b.r.get_end() - a.r.get_start();
//             assert(overlap > 0);
//             b.seq = (a.r.pos_strand()) ? 
//                 b.seq.substr(0, b.seq.length() - overlap) + string(overlap, 'N') :
//                 string(overlap, 'N') + b.seq.substr(overlap);
//         }
//     }
// }


inline static bool
same_read(const MappedRead &a, const MappedRead &b) 
{
    const string sa(a.r.get_name());
    const string sb(b.r.get_name());
    return std::equal(sa.begin(), sa.end() - 1, sb.begin());
}

inline static bool
name_smaller(const MappedRead &a, const MappedRead &b) 
{
    const string sa(a.r.get_name());
    const string sb(b.r.get_name());
    return std::lexicographical_compare(sa.begin(), sa.end() - 1, 
                                        sb.begin(), sb.end() - 1);
}


// static void
// mask_less_informative(MappedRead &one, MappedRead &two) 
// {
//     const size_t informative_one =
//         one.seq.length() - 
//         count(one.seq.begin(), one.seq.end(), 'N') -
//         static_cast<size_t>(one.r.get_score());
//     const size_t informative_two =
//         two.seq.length() - 
//         count(two.seq.begin(), two.seq.end(), 'N') -
//         static_cast<size_t>(two.r.get_score());
//     if (informative_one > informative_two) fix_overlap(one, two);
//     else fix_overlap(two, one);
// }

static void
merge_mates(const MappedRead &one, const MappedRead &two,
            MappedRead &merged, int &len,
            const int MAX_SEGMENT_LENGTH)
{
    const size_t info_one =
        one.seq.length() - 
        count(one.seq.begin(), one.seq.end(), 'N') -
        static_cast<size_t>(one.r.get_score());
    const size_t info_two =
        two.seq.length() - 
        count(two.seq.begin(), two.seq.end(), 'N') -
        static_cast<size_t>(two.r.get_score());
    
    merged = one;
    size_t start_one = 0, end_one = 0, start_two = 0, 
        end_two = 0, start_overlap = 0, end_overlap = 0;
    if (merged.r.pos_strand())
    {
        start_overlap = std::max(one.r.get_start(), two.r.get_start());
        end_overlap = std::min(one.r.get_end(), two.r.get_end());
        start_one = one.r.get_start();
        end_one = std::min(start_overlap, one.r.get_end());
        start_two = std::max(end_overlap, two.r.get_start());
        end_two = two.r.get_end();
        len = end_two - start_one;
        merged.r.set_start(start_one);
        merged.r.set_end(end_two);
    }
    else if (merged.r.neg_strand())
    {
        start_overlap = std::max(one.r.get_start(), two.r.get_start());
        end_overlap = std::min(one.r.get_end(), two.r.get_end());
        start_one = std::max(end_overlap, one.r.get_start());
        end_one = one.r.get_end();
        start_two = two.r.get_start();
        end_two = std::min(start_overlap, two.r.get_end());
        len = end_one - start_two;
        merged.r.set_start(start_two);
        merged.r.set_end(end_one);
    }
    assert(end_one >= start_one && end_two >= start_two);
    if (start_overlap < end_overlap )
        assert(static_cast<size_t>(len) == end_one - start_one
               + end_two - start_two + end_overlap - start_overlap);

    merged.r.set_score(one.r.get_score() + two.r.get_score());
    if (len > 0 && len < MAX_SEGMENT_LENGTH)
    {
        merged.seq = string(len, 'N');
        merged.scr = string(len, 'B');
        const string name = one.r.get_name();
        merged.r.set_name("FRAG:" + name.substr(0, name.size() - 1));
            
        std::copy(one.seq.begin(), one.seq.begin() + end_one - start_one,
                  merged.seq.begin());
        std::copy(one.scr.begin(), one.scr.begin() + end_one - start_one,
                  merged.scr.begin());
        std::copy(two.seq.end() + start_two - end_two, two.seq.end(),
                  merged.seq.end() + start_two - end_two);
        std::copy(two.scr.end() + start_two - end_two, two.scr.end(),
                  merged.scr.end() + start_two - end_two);

        // deal with overlapping part
        if (start_overlap < end_overlap && info_one >= info_two)
        {
            std::copy(one.seq.begin() + start_overlap - one.r.get_start(),
                      one.seq.begin() + end_overlap - one.r.get_start(),
                      merged.seq.begin() + end_one  - start_one);
            std::copy(one.scr.begin() + start_overlap - one.r.get_start(),
                      one.scr.begin() + end_overlap - one.r.get_start(),
                      merged.scr.begin() + end_one  - start_one);
        }
        else if (start_overlap < end_overlap && info_two > info_one)
        {
            std::copy(two.seq.begin() + start_overlap - two.r.get_start(),
                      two.seq.begin() + end_overlap -  two.r.get_start(),
                      merged.seq.begin() + end_one  - start_one);
            std::copy(two.scr.begin() + start_overlap -  two.r.get_start(),
                      two.scr.begin() + end_overlap - two.r.get_start(),
                      merged.scr.begin() + end_one  - start_one);
        }
    }
}


static inline string
collapse_mapped_reads(const MappedRead &mr,
                      const string delimiter = "\26")
{
    return
        mr.r.get_chrom() + delimiter
        + rmap::toa(mr.r.get_start()) + delimiter
        + rmap::toa(mr.r.get_end()) + delimiter
        + mr.r.get_name() + delimiter
        + rmap::toa(mr.r.get_score()) + delimiter
        + (mr.r.get_strand() == '+' ? "+" : "-") + delimiter
        + mr.seq + delimiter
        + mr.scr;
}

int 
main(int argc, const char **argv) 
{
    try 
    {

        int MAX_SEGMENT_LENGTH = 1000;
        
        bool VERBOSE = false;
        bool REVCOMP = true;
        string histogram_file;
        string end_one_file; 
        string end_two_file;
        string end_one_out;
        string end_two_out; 
        string out_stat;
        string outfile;
        string delimiter = "\26";
  
        /****************** COMMAND LINE OPTIONS ********************/
        OptionParser opt_parse(argv[0], "a program for identifying overlapping ends of "
                               "mapped paired end reads.",
                               "");
        opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
        opt_parse.add_opt("max-frag-len", 'L', "Maximum fragment length", 
                          false, MAX_SEGMENT_LENGTH); 
        opt_parse.add_opt("out_stat", 'S', "Name of file with statistics output", 
                          false, out_stat); 
        opt_parse.add_opt("T_rich_mates", 'T', "Name of input file with T-rich mates, mates1",
                          true, end_one_file); 
        opt_parse.add_opt("A_rich_mates", 'A', "Name of input file with A-rich mates, mates2",
                          true, end_two_file); 

        opt_parse.add_opt("outfile", 'o', "Output file name", 
                          false, outfile);
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

        /****************** END COMMAND LINE OPTIONS *****************/

        ifstream in_one(end_one_file.c_str());
        ifstream in_two(end_two_file.c_str());
    
        std::ostream *out = outfile.empty() ?
            &std::cout : new std::ofstream(outfile.c_str());

        /*************Used for collecting statistics***************/
        vector<size_t> fragm_len_distr(MAX_SEGMENT_LENGTH + 1, 0);
        size_t correctly_pairs = 0;
        size_t incorrectly_chr = 0;
        size_t incorrectly_strand = 0;
        size_t incorrectly_orient = 0;
        size_t incorrectly_fragm_size = 0;
        size_t broken_pairs = 0;
        /*************End for collecting statistics**************/

        MappedRead one, two, prev_one, prev_two;
        prev_one.r.set_name("");
        prev_two.r.set_name("");
        
        bool one_is_good = true, two_is_good = true;
        try { in_one >> one; check_sorted_by_ID(prev_one, one);}
        catch (const RMAPException &e) { one_is_good = false;}
        

        try { in_two >> two; check_sorted_by_ID(prev_two, two); }
        catch (const RMAPException &e) { two_is_good = false;}
        if (REVCOMP) revcomp(two);
        while (one_is_good && two_is_good) 
        {
            if (same_read(one, two)) // one and tow are mates
            {
                // if (one.r.overlaps(two.r))
                //     mask_less_informative(one, two);

                if (one.r.get_chrom() != two.r.get_chrom())
                {   
                    incorrectly_chr++;
                    *out << one << endl << two << endl;
                }
                else if (one.r.get_strand() != two.r.get_strand())
                {
                    incorrectly_strand++;
                    *out << one << endl << two << endl;
                }
                else
                {
                    MappedRead merged;
                    int len;
                    merge_mates(one, two, merged, len, MAX_SEGMENT_LENGTH);
                    if (len > 0 && len < MAX_SEGMENT_LENGTH){
                        correctly_pairs++;
                        fragm_len_distr[len]++;
                        *out << merged << endl;
                    }
                    else{
                        *out << one << endl << two << endl;
                        if(len < 0 )
                            incorrectly_orient++;
                        else
                            incorrectly_fragm_size++;
                    }

                }
                
                try { in_one >> one; check_sorted_by_ID(prev_one, one); }
                catch (const RMAPException &e) { one_is_good = false;}
                try { in_two >> two; check_sorted_by_ID(prev_two, two); }
                catch (const RMAPException &e) { two_is_good = false;}
                if (REVCOMP) revcomp(two);
            } 
            else if (name_smaller(one, two))
            {
                *out << one << endl;
                broken_pairs++;

                try { in_one >> one; check_sorted_by_ID(prev_one, one); }
                catch (const RMAPException &e) { one_is_good = false;}
            }
            else // one comes after two
            {
                *out << two << endl;
                broken_pairs++;

                try { in_two >> two; check_sorted_by_ID(prev_two, two); }
                catch (const RMAPException &e) { two_is_good = false;}
                if (REVCOMP) revcomp(two);
            }
        }
        while (one_is_good) 
        {
            *out << one << endl;
            broken_pairs++;            
            
            try { in_one >> one; check_sorted_by_ID(prev_one, one); }
            catch (const RMAPException &e) { one_is_good = false;}
        }
        while (two_is_good) 
        {
            *out << two << endl;
            broken_pairs++;
            
            try { in_two >> two; check_sorted_by_ID(prev_two, two); }
            catch (const RMAPException &e) { two_is_good = false;}
            if (REVCOMP) revcomp(two);
        }

        if(!out_stat.empty()){
            ofstream outst(out_stat.c_str());
            outst << "TOTAL CORRECTLY MAPPED PAIRS (count in pairs):\t" << correctly_pairs << endl;
            outst << "INCORRECTLY MAPPED TO DIFFERENT CHROM:\t" << incorrectly_chr << endl;
            outst << "INCORRECTLY MAPPED DUE TO STRAND INCOMPATIBILITY:\t" << incorrectly_strand << endl;
            outst << "INCORRECTLY MAPPED DUE TO ORIENTATION:\t" << incorrectly_orient << endl;
            outst << "INCORRECTLY MAPPED DUE TO FRAGMENT SIZE:\t" << incorrectly_fragm_size << endl;
            outst << "TOTAL MAPPED BROKEN PAIRS (with missing mates, single-end count):\t" << broken_pairs << endl;
            outst << "FRAGM_LEN:\t" << "PAIRS_WITH_THIS_FRAGM_SIZE:" << endl;
            long int i = 0;
            for(i = 0; i < MAX_SEGMENT_LENGTH + 1; i++){
                outst << i << "\t" << fragm_len_distr[i] << endl;
            }//for
        }//if stat file
        if (out != &cout) delete out;     
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

