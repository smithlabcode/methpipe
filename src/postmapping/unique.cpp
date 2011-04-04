/*    unique: a program to select unique reads mapped to the same location  
 *    and the same strand  
 *
 *    Copyright (C) 2009 University of Southern California and
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

#include <string>
#include <vector>
#include <unistd.h>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <queue>

#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"
#include "RNG.hpp"
#include "FileIterator.hpp"
#include "MappedRead.hpp"
#include "QualityScore.hpp"

#include <sys/time.h>
#include <sys/resource.h>

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::numeric_limits;
using std::ofstream;


static inline bool 
is_duplicate_read(const MappedRead &lhs,
                  const MappedRead &rhs)
{
    return
        lhs.r.get_chrom() == rhs.r.get_chrom() &&
        lhs.r.get_start() == rhs.r.get_start() &&
        lhs.r.get_strand() == rhs.r.get_strand();
}

static inline bool
is_paired_fragment(const MappedRead &mr)
{
    return mr.r.get_name() == "FRAG:PAIRED";
}

// This is really messy. This is approximately right
// The unpaired reads on the negative strand is really messy
// and this progrma tends to retain them
static inline bool 
is_duplicate_fragment(const MappedRead &lhs,
                      const MappedRead &rhs)
{
    if (is_paired_fragment(lhs) && is_paired_fragment(rhs))
    {
        return 
            lhs.r.get_chrom() == rhs.r.get_chrom() &&
            lhs.r.get_start() == rhs.r.get_start() &&
            lhs.r.get_end() == rhs.r.get_end() &&
            lhs.r.get_strand() == rhs.r.get_strand();
    }
    else if (!is_paired_fragment(lhs) && is_paired_fragment(rhs))
    {
        return
            (lhs.r.get_strand() == '+' &&   // being strigent on positive reads
             rhs.r.get_strand() == '+' &&
             lhs.r.get_chrom() == rhs.r.get_chrom()  &&
             lhs.r.get_start() == rhs.r.get_start())
            ||
            (lhs.r.get_strand() == '-' &&  // being a little conservative on
             rhs.r.get_strand() == '-' &&  // negative reads
             lhs.r.get_chrom() == rhs.r.get_chrom()  &&
             lhs.r.get_end() == rhs.r.get_end());
    } 
    else if (is_paired_fragment(lhs) && !is_paired_fragment(rhs))
    {
        return
            (lhs.r.get_strand() == '+' &&   // being strigent on positive reads
             rhs.r.get_strand() == '+' &&
             lhs.r.get_chrom() == rhs.r.get_chrom()  &&
             lhs.r.get_start() == rhs.r.get_start())
            ||
            (lhs.r.get_strand() == '-' &&  // being a little conservative on
             rhs.r.get_strand() == '-' &&  // negative reads
             lhs.r.get_chrom() == rhs.r.get_chrom()  &&
             lhs.r.get_end() == rhs.r.get_end());
    } 
    else // (!is_paired_fragment(lhs) && !is_paired_fragment(rhs))
    {
        return // the same behavior as without fragments
            lhs.r.get_chrom() == rhs.r.get_chrom() &&
            lhs.r.get_start() == rhs.r.get_start() &&
            lhs.r.get_strand() == rhs.r.get_strand();
    }
}


void
consensus_mappedread(vector<MappedRead> &mapped_ties,
                     MappedRead &consensus_mr)
{
    static const size_t random_number_seed = numeric_limits<size_t>::max();
    static const Runif rng(random_number_seed);

    vector<MappedRead>::iterator iter =
        std::partition(mapped_ties.begin(), mapped_ties.end(),
                       is_paired_fragment);
    if (iter != mapped_ties.begin())
    {
        const size_t rand_idx = rng.runif(0, iter - mapped_ties.begin());
        consensus_mr = mapped_ties[rand_idx];
    }
    else
    {
        const size_t rand_idx =
            rng.runif(0, static_cast<int>(mapped_ties.size()));
        consensus_mr = mapped_ties[rand_idx];
    }
}

// void
// consensus_mappedread(const vector<MappedRead> &mapped_ties,
//                MappedRead &consensus_mr)
// {
//     consensus_mr.r = mapped_ties.front().r;
//     const size_t num = mapped_ties.size();
//     const size_t read_len = mapped_ties.front().seq.size();

    
    
//     for (size_t i = 0; i < read_len; ++i)
//     {
//         vector<double> quality_scores(5, 0);
//         for (size_t j = 0; j < num; ++j)
//         {
//             const int base_id = base2int(mapped_ties[j].seq[i]);
//             const double prob = char2prob_solexa(mapped_ties[j].scr[i]);
//             quality_scores[base_id] += prob;
//         }
        
//         const double all_prob = std::accumulate(quality_scores.begin(),
//                                                 quality_scores.end(),
//                                                 0.0);
//         const size_t idx =
//             std::max_element(quality_scores.begin(),
//                              quality_scores.end())
//             - quality_scores.begin();
        
//         consensus_mr.seq[i] = int2base(idx);
//         consensus_mr.scr[i] = prob2char_solexa(quality_scores[idx] / all_prob);
//     }
// }


int 
main(int argc, const char **argv) 
{

    try 
    {

        bool VERBOSE = false;
        string infile;
        string outfile;
        
        /****************** COMMAND LINE OPTIONS ********************/
        OptionParser opt_parse(argv[0], 
                               "A program to select a unique read from multiple reads "
                               "mapped to the same location and the same strand"
                               "<file with names of bed-files>");
        opt_parse.add_opt("output", 'o', "Name of maps output file", 
                          false, outfile);
        opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
        vector<string> leftover_args;
        opt_parse.parse(argc, argv, leftover_args);
        if (opt_parse.help_requested()) 
        {
            cerr << opt_parse.help_message() << endl
                 << opt_parse.about_message() << endl;
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
        
        if (leftover_args.size() > 0)
            infile = leftover_args.front();
        
        std::istream *in = infile.empty() ?
            &std::cin : new std::ifstream(infile.c_str());

        std::ostream *out = (outfile.empty()) ?
            &std::cout : new ofstream(outfile.c_str());
    
        vector<MappedRead> mapped_ties;
        MappedRead mr;

        while (in->good() && (*in >> mr)) 
        {
            if (!mapped_ties.empty() && 
                !is_duplicate_fragment(mapped_ties.front(), mr)) 
            {
                static MappedRead consensus_mr(mapped_ties.front());
                consensus_mappedread(mapped_ties, consensus_mr);
                *out << consensus_mr << endl;
                mapped_ties.clear();
            }
            mapped_ties.push_back(mr);
        }
        if (!mapped_ties.empty()) 
        {
            static  MappedRead consensus_mr(mapped_ties.front());
            consensus_mappedread(mapped_ties, consensus_mr);
            *out << consensus_mr << endl;
        }

        if (in != &std::cin) delete in;
        if (out != &std::cout) delete out;
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
