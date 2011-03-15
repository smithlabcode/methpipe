/*    checkoverlap: a program for identifying overlapping ends of
 *                  mapped paired end reads.
 *
 *    Copyright (C) 2010 University of Southern California and
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

#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "MappedRead.hpp"

using namespace std;

inline static size_t 
get_distance(const MappedRead &a, const MappedRead &b) 
{
    return (a.r.pos_strand()) ? b.r.get_end() - a.r.get_start() :
        a.r.get_end() - b.r.get_start();
}


static void
fix_overlap(const MappedRead &a, MappedRead &b) 
{
    if (a.r.get_strand() == b.r.get_strand()) 
    {
        if (a.r.get_start() < b.r.get_start()) 
        {
            const size_t overlap = a.r.get_end() - b.r.get_start();
            assert(overlap > 0);
            b.seq = (a.r.pos_strand()) ?
                string(overlap, 'N') + b.seq.substr(overlap) :
                b.seq.substr(0, b.seq.length() - overlap) + string(overlap, 'N');
        }
        else 
        {
            const size_t overlap = b.r.get_end() - a.r.get_start();
            assert(overlap > 0);
            b.seq = (a.r.pos_strand()) ? 
                b.seq.substr(0, b.seq.length() - overlap) + string(overlap, 'N') :
                string(overlap, 'N') + b.seq.substr(overlap);
        }
    }
}


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


static void
mask_less_informative(MappedRead &one, MappedRead &two) 
{
    const size_t informative_one =
        one.seq.length() - 
        count(one.seq.begin(), one.seq.end(), 'N') -
        static_cast<size_t>(one.r.get_score());
    const size_t informative_two =
        two.seq.length() - 
        count(two.seq.begin(), two.seq.end(), 'N') -
        static_cast<size_t>(two.r.get_score());
    if (informative_one > informative_two) fix_overlap(one, two);
    else fix_overlap(two, one);
}

inline size_t
frag_len(const MappedRead &lhs, const MappedRead &rhs)
{
    return max(lhs.r.get_end(), rhs.r.get_end()) - 
        min(lhs.r.get_start(), rhs.r.get_start());
}

int 
main(int argc, const char **argv) 
{
    try 
    {
        bool VERBOSE = false;
        bool REVCOMP = false;
        size_t max_distance = 500;
        string histogram_file;
        string end_one_file; 
        string end_two_file;
        string end_one_out;
        string end_two_out; 
        string fraglen_file;
  
        /****************** COMMAND LINE OPTIONS ********************/
        OptionParser opt_parse(argv[0], "a program for identifying overlapping ends of "
                               "mapped paired end reads.",
                               "<end-1-in> <end-2-in> <end-1-out> <end-2-out>");
        opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
        opt_parse.add_opt("ahist", 'h', "print frag size histogram in this file", 
                          false, histogram_file);
        opt_parse.add_opt("revcomp", 'r',
                          "Reverse complement A-rich strand before masking", 
                          false, histogram_file);
        opt_parse.add_opt("max-dist", 'm', "max distance to print", 
                          false, max_distance);
        opt_parse.add_opt("frag-len", 'f', "File name of fragment length", 
                          false, fraglen_file);
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
        if (leftover_args.size() < 2) 
        {
            cerr << opt_parse.help_message() << endl;
            return EXIT_SUCCESS;
        }

        end_one_file = leftover_args[0]; 
        end_two_file = leftover_args[1];
        end_one_out = leftover_args[2];
        end_two_out = leftover_args[3]; 
        /****************** END COMMAND LINE OPTIONS *****************/

        ifstream in_one(end_one_file.c_str());
        ifstream in_two(end_two_file.c_str());
    
        ofstream out_one(end_one_out.c_str());
        ofstream out_two(end_two_out.c_str());
        
        ofstream out_fraglen;
        if (!fraglen_file.empty())
            out_fraglen.open(fraglen_file.c_str());

        bool is_paired = false;
        MappedRead one, two;
    
        in_one >> one;
        in_two >> two;
        if (REVCOMP) revcomp(two);
        while (in_one.good() && in_two.good()) 
        {
            if (same_read(one, two)) // one and tow are mates
            {
                if (out_fraglen.is_open())
                {
                    if (one.r.get_chrom() != two.r.get_chrom())
                        out_fraglen << one.r.get_name() << "\t"
                                    << "WARNING: different chromosomes" << endl;
                    else if (one.r.get_strand() != two.r.get_strand())
                        out_fraglen << one.r.get_name() << "\t"
                                    << "WARNING: different strands" << endl;
                    else
                        out_fraglen << one.r.get_name() << "\t"
                                    << frag_len(one, two) << endl;
                }
                
                if (one.r.overlaps(two.r))
                    mask_less_informative(one, two);
                
                out_one << one << endl; 
                out_two << two << endl;
                
                in_one >> one;
                in_two >> two;
                if (REVCOMP) revcomp(two);
            } 
            else
            {
                if (name_smaller(one, two))
                {
                    out_one << one << endl;
                    if (out_fraglen.is_open())
                        out_fraglen << one.r.get_name() << "\t"
                                    << "WARNING: missed A-rich read" << endl;
                    in_one >> one;
                }
                else
                {
                    out_two << two << endl;
                    if (out_fraglen.is_open())
                        out_fraglen << two.r.get_name() << "\t"
                                    << "WARNING: missed T-rich read" << endl;
                    in_two >> two;
                    if (REVCOMP) revcomp(two);
                }
            }
        }
        while (in_one.good()) 
        {
            out_one << one << endl;
            if (out_fraglen.is_open())
                out_fraglen << one.r.get_name() << "\t"
                            << "WARNING: missed A-rich read" << endl;
            in_one >> one;
        }
        while (in_two.good()) 
        {
            out_two << two << endl;
            if (out_fraglen.is_open())
                out_fraglen << two.r.get_name() << "\t"
                            << "WARNING: missed T-rich read" << endl;
            in_two >> two;
            if (REVCOMP) revcomp(two);
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
