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

inline string
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
        static const int WARNING_DIFF_CHROM = -1;
        static const int WARNING_DIFF_STRAND = -2;
        static const int WARNING_MISS_T_MATE = -3;
        static const int WARNING_MISS_A_MATE = -4;
        
        bool VERBOSE = false;
        bool REVCOMP = false;
        size_t max_distance = 500;
        string histogram_file;
        string end_one_file; 
        string end_two_file;
        string end_one_out;
        string end_two_out; 
        string fraglen_file;
        string outfile;
        string delimiter = "\26";
  
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
                          false, outfile); // to be discontinued. only for 
                                           // backward compatibilities
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
        if (leftover_args.size() < 2) 
        {
            cerr << opt_parse.help_message() << endl;
            return EXIT_SUCCESS;
        }

        end_one_file = leftover_args[0]; 
        end_two_file = leftover_args[1];
        end_one_out = leftover_args.size() >= 3 ? leftover_args[2] : "/dev/null";
        end_one_out = leftover_args.size() >= 4 ? leftover_args[3] : "/dev/null";
        /****************** END COMMAND LINE OPTIONS *****************/

        ifstream in_one(end_one_file.c_str());
        ifstream in_two(end_two_file.c_str());
    
        ofstream out_one(end_one_out.c_str());
        ofstream out_two(end_two_out.c_str());
        
        std::ostream *out = outfile.empty() ?
            &std::cout : new std::ofstream(outfile.c_str());
        
        MappedRead one, two;
        bool one_is_good = true, two_is_good = true;
        
        try { in_one >> one; }
        catch (const RMAPException &e) { one_is_good = false;}
        try { in_two >> two; }
        catch (const RMAPException &e) { two_is_good = false;}
        if (REVCOMP) revcomp(two);
        while (one_is_good && two_is_good) 
        {
            if (same_read(one, two)) // one and tow are mates
            {
                if (one.r.overlaps(two.r))
                    mask_less_informative(one, two);

                out_one << one << endl; // for backward compatibilities
                out_two << two << endl;
                
                if (one.r.get_chrom() != two.r.get_chrom())
                {
                    *out << one.r.get_chrom() << "\t"
                         << one.r.get_start() << "\t"
                         << one.r.get_end() << "\t"
                         <<  "FRAG:UNPAIRED" << "\t"
                         << WARNING_DIFF_CHROM << "\t"
                         << one.r.get_strand() << "\t"
                         << collapse_mapped_reads(one) << "\t"
                         << "different-chromosomes" << endl;
                    *out << two.r.get_chrom() << "\t"
                         << two.r.get_start() << "\t"
                         << two.r.get_end() << "\t"
                         <<  "FRAG:UNPAIRED" << "\t"
                         << WARNING_DIFF_CHROM << "\t"
                         << two.r.get_strand() << "\t"
                         << collapse_mapped_reads(two) << "\t"
                         << "different-chromosomes" << endl;
                }
                else if (one.r.get_strand() != two.r.get_strand())
                {
                    *out << one.r.get_chrom() << "\t"
                         << one.r.get_start() << "\t"
                         << one.r.get_end() << "\t"
                         <<  "FRAG:UNPAIRED" << "\t"
                         << WARNING_DIFF_STRAND << "\t"
                         << one.r.get_strand() << "\t"
                         << collapse_mapped_reads(one) << "\t"
                         << "different-strands" << endl;
                    *out << two.r.get_chrom() << "\t"
                         << two.r.get_start() << "\t"
                         << two.r.get_end() << "\t"
                         <<  "FRAG:UNPAIRED" << "\t"
                         << WARNING_DIFF_STRAND << "\t"
                         << two.r.get_strand() << "\t"
                         << collapse_mapped_reads(two) << "\t"
                         << "different-strands" << endl;
                }
                else
                {
                    *out << one.r.get_chrom() << "\t"
                         << std::min(one.r.get_start(),
                                     two.r.get_start()) << "\t"
                         << std::max(one.r.get_end(),
                                     two.r.get_end()) << "\t"
                         << "FRAG:PAIRED" << "\t"
                         << frag_len(one, two) << "\t"
                         << one.r.get_strand() << "\t"
                         << collapse_mapped_reads(one) << "\t"
                         << collapse_mapped_reads(two) << endl;
                }
                
                try { in_one >> one; }
                catch (const RMAPException &e) { one_is_good = false;}
                try { in_two >> two; }
                catch (const RMAPException &e) { two_is_good = false;}
                if (REVCOMP) revcomp(two);
            } 
            else if (name_smaller(one, two))
            {
                *out << one.r.get_chrom() << "\t"
                     << one.r.get_start() << "\t"
                     << one.r.get_end() << "\t"
                     <<  "FRAG:UNPAIRED" << "\t"
                     << WARNING_MISS_A_MATE << "\t"
                     << one.r.get_strand() << "\t"
                     << collapse_mapped_reads(one) << "\t"
                     << "miss-A-mate" << endl;
                
                out_one << one << endl; // for backward compatibilities

                try { in_one >> one; }
                catch (const RMAPException &e) { one_is_good = false;}
            }
            else // one comes after two
            {
                *out << two.r.get_chrom() << "\t"
                     << two.r.get_start() << "\t"
                     << two.r.get_end() << "\t"
                     <<  "FRAG:UNPAIRED" << "\t"
                     << WARNING_MISS_T_MATE << "\t"
                     << two.r.get_strand() << "\t"
                     << collapse_mapped_reads(two) << "\t"
                     << "miss-T-mate" << endl;
                out_two << two << endl; // for backward compatibilities
                
                try { in_two >> two; }
                catch (const RMAPException &e) { two_is_good = false;}
                if (REVCOMP) revcomp(two);
            }
        }
        while (one_is_good) 
        {
            *out << one.r.get_chrom() << "\t"
                 << one.r.get_start() << "\t"
                 << one.r.get_end() << "\t"
                 <<  "FRAG:UNPAIRED" << "\t"
                 << WARNING_MISS_A_MATE << "\t"
                 << one.r.get_strand() << "\t"
                 << collapse_mapped_reads(one) << "\t"
                 << "miss-A-mate" << endl;
            out_one << one << endl; // for backward compatibilities
            
            try { in_one >> one; }
            catch (const RMAPException &e) { one_is_good = false;}
        }
        while (two_is_good) 
        {
            *out << two.r.get_chrom() << "\t"
                 << two.r.get_start() << "\t"
                 << two.r.get_end() << "\t"
                 <<  "FRAG:UNPAIRED" << "\t"
                 << WARNING_MISS_T_MATE << "\t"
                 << two.r.get_strand() << "\t"
                 << collapse_mapped_reads(two) << "\t"
                 << "miss-T-mate" << endl;
            out_two << two << endl; // for backward compatibilities
            
            try { in_two >> two; }
            catch (const RMAPException &e) { two_is_good = false;}
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

