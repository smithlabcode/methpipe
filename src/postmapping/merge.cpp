/*    mergelanes: a program for sorting read sequences relative to an order
 *    given in a BED file (i.e. mapping locations)
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
//      lhs.r.get_end() == rhs.r.get_end() &&
        lhs.r.get_strand() == rhs.r.get_strand();
}

static inline bool 
is_duplicate_read(const GenomicRegion &lhs,
                  const GenomicRegion &rhs)
{
    return
        lhs.get_chrom() == rhs.get_chrom() &&
        lhs.get_start() == rhs.get_start() &&
//      lhs.get_end() == rhs.r.get_end() &&
        lhs.get_strand() == rhs.get_strand();
}

void
consensus_mappedread(const vector<MappedRead> &mapped_ties,
               MappedRead &consensus_mr)
{
    consensus_mr.r = mapped_ties.front().r;
    const size_t num = mapped_ties.size();
    const size_t read_len = mapped_ties.front().seq.size();

    
    
    for (size_t i = 0; i < read_len; ++i)
    {
        vector<double> quality_scores(5, 0);
        for (size_t j = 0; j < num; ++j)
        {
            const int base_id = base2int(mapped_ties[j].seq[i]);
            const double prob = char2prob_solexa(mapped_ties[j].scr[i]);
            quality_scores[base_id] += prob;
        }
        
        const double all_prob = std::accumulate(quality_scores.begin(),
                                                quality_scores.end(),
                                                0.0);
        const size_t idx =
            std::max_element(quality_scores.begin(),
                             quality_scores.end())
            - quality_scores.begin();
        
        consensus_mr.seq[i] = int2base(idx);
        consensus_mr.scr[i] = prob2char_solexa(quality_scores[idx] / all_prob);
    }
}

struct ComparePairs 
{
    bool operator()(const pair<GenomicRegion, size_t> &a,
                    const pair<GenomicRegion, size_t> &b) const 
    {
        const GenomicRegion &lhs  = a.first;
        const GenomicRegion &rhs  = b.first;
        return (lhs.get_chrom() < rhs.get_chrom()) ||
            (lhs.get_chrom() == rhs.get_chrom() &&
             lhs.get_start() < rhs.get_start()) || 
            (lhs.get_chrom() == rhs.get_chrom() &&
             lhs.get_start() == rhs.get_start() &&
             lhs.get_strand() < rhs.get_strand());
    }
};

int 
main(int argc, const char **argv) 
{

    try 
    {

        bool VERBOSE = false;
        string map_outfile;
        string read_outfile;
        
        size_t BUFFER_SIZE = 10000ul;
        
        bool ALLOW_DUPLICATES = false;
        
        /****************** COMMAND LINE OPTIONS ********************/
        OptionParser opt_parse(argv[0], 
                               "A program to merge different lanes in a shotgun bisulfite "
                               "sequencing experiment producing one large sorted file of all read "
                               "mapping locations (without duplicate 5' ends) and one large file "
                               "sorted similarly containing the corresponding sequences.",
                               "<file with names of bed-files>");
        opt_parse.add_opt("output", 'o', "Name of maps output file", 
                           false, map_outfile);
        opt_parse.add_opt("dups", 'D', "Allow duplicate fragments",
                          false, ALLOW_DUPLICATES);
        opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
        vector<string> leftover_args;
        opt_parse.parse(argc, argv, leftover_args);
        if (argc == 1 || opt_parse.help_requested()) 
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
        
        vector<string> mapped_files;
        const string mapped_file_names_file(leftover_args.front());
        if (mapped_file_names_file.find(".mr") != string::npos
            || mapped_file_names_file.find(".bed") != string::npos)
            std::copy(leftover_args.begin(), leftover_args.end(),
                      std::back_inserter(mapped_files));
        else 
            read_filename_file(mapped_file_names_file.c_str(), mapped_files);

        vector<FileIterator<MappedRead> *> itrs;
        if (VERBOSE)
            cerr << "[MAPPED READS FILES:]" << endl;
        for (size_t i = 0; i < mapped_files.size(); ++i) 
        {
            if (VERBOSE)
                cerr << mapped_files[i] << endl;
            itrs.push_back(new FileIterator<MappedRead>(mapped_files[i], BUFFER_SIZE));
        }

        std::priority_queue<pair<GenomicRegion, size_t>, 
            vector<pair<GenomicRegion, size_t> >, ComparePairs> a;
        for (size_t i = 0; i < itrs.size(); ++i)
            a.push(make_pair((*itrs[i]->get_first()).r, i));
    
        std::ostream *out = (map_outfile.empty()) ?
            &cout : new ofstream(map_outfile.c_str());
    
        vector<MappedRead> mapped_ties;
        
        while (!a.empty()) 
        {
            const size_t file_id = a.top().second;
            if (ALLOW_DUPLICATES) 
            {
                const MappedRead mr = (*itrs[file_id]->get_first());
                *out << mr << endl;
            }
            else 
            {
                if (!mapped_ties.empty() && 
                    !is_duplicate_read(mapped_ties.front().r, a.top().first)) 
                {
                    static MappedRead mr(mapped_ties.front());
                    consensus_mappedread(mapped_ties, mr);
                    *out << mr << endl;
                    mapped_ties.clear();
                }
                mapped_ties.push_back((*itrs[file_id]->get_first()));
            }
            a.pop();
            itrs[file_id]->increment();
            if (itrs[file_id]->first_is_good())
                a.push(make_pair((*itrs[file_id]->get_first()).r, file_id));
        }
        if (!ALLOW_DUPLICATES) 
        {
            static  MappedRead mr(mapped_ties.front());
            consensus_mappedread(mapped_ties, mr);
            *out << mr << endl;
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
