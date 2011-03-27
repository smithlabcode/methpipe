/*    merge-counts.cpp: 
 *    This program is used to merge multiple output file of methcounts
 *
 *    Copyright (C) 2011 University of Southern California and
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
is_same_region(const GenomicRegion &lhs,
               const GenomicRegion &rhs)
{
    return
        lhs.get_chrom() == rhs.get_chrom() &&
        lhs.get_start() == rhs.get_start() &&
//      lhs.get_end() == rhs.r.get_end() &&
        lhs.get_strand() == rhs.get_strand();
}

static inline bool
is_same_cytosine(const GenomicRegion &lhs,
                 const GenomicRegion &rhs)
{
    return is_same_region(lhs, rhs);
}

void
merge_same_cytosine(const vector<GenomicRegion> &same_cytosines,
                    GenomicRegion &cytosine)
{
    double meth = 0, unmeth = 0;
    for (size_t i = 0; i < same_cytosines.size(); ++i)
    {
        const string name = same_cytosines[i].get_name();
        const double n_reads =
            atof(name.substr(name.find_first_of(":") + 1).c_str());
        const double meth_prob = same_cytosines[i].get_score();
        const double n_meth = n_reads * meth_prob;
        const double n_unmeth = n_reads - n_meth;
        meth += n_meth;
        unmeth += n_unmeth;
    }

    cytosine = same_cytosines.front();
    const string name = cytosine.get_name();
    cytosine.set_name(name.substr(0, name.find_first_of(":"))
                      +  ":"
                      + toa(static_cast<size_t>(meth + unmeth)));
    cytosine.set_score(meth + unmeth >= 1.0 ? meth / (meth + unmeth) : 0);
}


struct ComparePairs 
{
    bool operator()(const pair<GenomicRegion, size_t> &a,
                    const pair<GenomicRegion, size_t> &b) const 
    {
        const GenomicRegion &lhs  = a.first;
        const GenomicRegion &rhs  = b.first;
        return (lhs.get_chrom() > rhs.get_chrom()) ||
            (lhs.get_chrom() == rhs.get_chrom() &&
             lhs.get_start() > rhs.get_start()) || 
            (lhs.get_chrom() == rhs.get_chrom() &&
             lhs.get_start() == rhs.get_start() &&
             lhs.get_strand() > rhs.get_strand());
    }
};

int 
main(int argc, const char **argv) 
{

    try 
    {

        bool VERBOSE = false;
        string outfile;
        
        size_t BUFFER_SIZE = 10000ul;
        
        
        /****************** COMMAND LINE OPTIONS ********************/
        OptionParser opt_parse(argv[0], 
                               "A program to merge multiple methcount output files"
                               "<methcount output bedfiles>");
        opt_parse.add_opt("output", 'o', "Name of output file", 
                           false, outfile);
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
        
        vector<string> methcount_files(leftover_args);

        vector<FileIterator<GenomicRegion> *> itrs;
        if (VERBOSE)
            cerr << "[Input methcount files:]" << endl;
        for (size_t i = 0; i < methcount_files.size(); ++i) 
        {
            if (VERBOSE)
                cerr << methcount_files[i] << endl;
            itrs.push_back(
                new FileIterator<GenomicRegion>(methcount_files[i], BUFFER_SIZE));
        }

        std::priority_queue<pair<GenomicRegion, size_t>, 
            vector<pair<GenomicRegion, size_t> >, ComparePairs> a;
        for (size_t i = 0; i < itrs.size(); ++i)
            a.push(make_pair(*itrs[i]->get_first(), i));
    
        std::ostream *out = (outfile.empty()) ?
            &cout : new ofstream(outfile.c_str());
    
        vector<GenomicRegion> same_cytosines;
        
        while (!a.empty()) 
        {
            const size_t file_id = a.top().second;
            if (!same_cytosines.empty() && 
                !is_same_cytosine(same_cytosines.front(), a.top().first)) 
            {
                GenomicRegion cytosine;
                merge_same_cytosine(same_cytosines, cytosine);
                *out << cytosine << endl;
                same_cytosines.clear();
            }
            same_cytosines.push_back((*itrs[file_id]->get_first()));
            a.pop();
            itrs[file_id]->increment();
            if (itrs[file_id]->first_is_good())
                a.push(make_pair(*itrs[file_id]->get_first(), file_id));
        }
        if (!same_cytosines.empty()) 
        {
            GenomicRegion cytosine;
            merge_same_cytosine(same_cytosines, cytosine);
            *out << cytosine << endl;
            same_cytosines.clear();
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
