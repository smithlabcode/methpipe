/*    extract cytosines in typical sequence context
 *    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith, Song Qiang
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
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cctype>

#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"

using namespace std;

void
find_cytosines(const string & sequences,
               const string & outfile,
               const string & chrom,
               const bool CYTOSINE_INDEX)
{
    std::ofstream outf(outfile.c_str());
    size_t nc = 0;
    for (size_t i = 0; i < sequences.size(); ++i)
        if (sequences[i] == 'C')
        {
            ++nc;
            outf << chrom << "\t" << i << "\t" << i + 1 << "\t"
                 << "C:" << (CYTOSINE_INDEX ? nc : 0) << "\t"
                 << 0 << "\t" << "+" << endl;
        }
        else if (sequences[i] == 'G')
        {
            ++nc;
            outf << chrom << "\t" << i << "\t" << i + 1 << "\t"
                 << "C:" << (CYTOSINE_INDEX ? nc : 0) << "\t"
                 << 0 << "\t" << "-" << endl;
        }
    outf.close();
}

void
find_cytosines_CG(const string & sequences,
                  const string & outfile,
                  const string & chrom,
                  const bool CYTOSINE_INDEX)
{
    std::ofstream outf(outfile.c_str());
    size_t ncg = 0;
    for (size_t i = 0; i < sequences.size() - 1; ++i)
        if (sequences[i] == 'C' && sequences[i + 1] == 'G')
        {
            ++ncg;
            outf << chrom << "\t" << i << "\t" << i + 1 << "\t"
                 << "CpG:" << (CYTOSINE_INDEX ? ncg : 0) << "\t"
                 << 0 << "\t" << "+" << endl;
            ++ncg;
            outf << chrom << "\t" << i + 1 << "\t" << i + 2 << "\t"
                 << "CpG:" << (CYTOSINE_INDEX ? ncg : 0) << "\t"
                 << 0 << "\t" << "-" << endl;
        }
    
    outf.close();
}

void
find_cytosines_CHG(const string & sequences,
                   const string & outfile,
                   const string & chrom,
                   const bool CYTOSINE_INDEX)
{
    std::ofstream outf(outfile.c_str());
    size_t nchg = 0;
    for (size_t i = 0; i < sequences.size(); ++i)
        if (i + 2 < sequences.size() &&  sequences[i] == 'C'
            && sequences[i + 1] != 'G' && sequences[i + 2] == 'G')
        {
            ++nchg;
            outf << chrom << "\t" << i << "\t" << i + 1 << "\t"
                 << "C" << sequences[i+1] << "G:"
                 << (CYTOSINE_INDEX ? nchg : 0) << "\t"
                 << 0 << "\t" << "+" << endl;
        } 
        else if (i >= 2 && sequences[i] == 'G'
                 && sequences[i - 1] != 'C' && sequences[i - 2] == 'C')
        {
            ++nchg;
            outf << chrom << "\t" << i << "\t" << i + 1 << "\t"
                 << "C" << complement(sequences[i-1]) << "G:"
                 << (CYTOSINE_INDEX ? nchg : 0) << "\t"
                 << 0 << "\t" << "-" << endl;
        } 
    
    outf.close();
}

void
find_cytosines_CHH(const string & sequences,
                   const string & outfile,
                   const string & chrom,
                   const bool CYTOSINE_INDEX)
{
    std::ofstream outf(outfile.c_str());
    size_t nchh = 0;
    for (size_t i = 0; i < sequences.size(); ++i)
        if (i + 2 < sequences.size() &&  sequences[i] == 'C'
            && sequences[i + 1] != 'G' && sequences[i + 2] != 'G')
        {
            ++nchh;
            outf << chrom << "\t" << i << "\t" << i + 1 << "\t"
                 << "C" << sequences[i+1] << sequences[i+2]
                 << ":"  << (CYTOSINE_INDEX ? nchh : 0) << "\t" 
                 << 0 << "\t" << "+" << endl;
        } 
        else if (i >= 2 && sequences[i] == 'G'
                 && sequences[i - 1] != 'C' && sequences[i - 2] != 'C')
        {
            ++nchh;
            outf << chrom << "\t" << i << "\t" << i + 1 << "\t"
                 << "C" << complement(sequences[i-1])
                 << complement(sequences[i-2]) << ":"
                 << (CYTOSINE_INDEX ? nchh : 0) << "\t"
                 << 0 << "\t" << "-" << endl;
        } 
    
    outf.close();
}

int 
main(int argc, const char **argv) {

    static const size_t max_line_size = 1000;
    
    try {
        bool VERBOSE = false;
        string infile, outfile;
        // string pattern = "CG";
        string cfile, cgfile, chgfile, chhfile;

        bool CYTOSINE_INDEX = true;

        // size_t ws = 200;

        /****************** GET COMMAND LINE ARGUMENTS ***************************/
        OptionParser opt_parse("CpG", "Compute CpG of genomes");
        opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
                          false , outfile);
        opt_parse.add_opt("input", 'i', "Name of input file (default: stdin)", 
                          false , infile);
        opt_parse.add_opt("cytosine", '\0', "Name of output file for cytosines", 
                          false , cfile);
        opt_parse.add_opt("cg-cytosine", '\0', "Name of output file for cytosines in CG", 
                          false , cgfile);
        opt_parse.add_opt("chg-cytosine", '\0', "Name of output file for cytosines in CHG", 
                          false , chgfile);
        opt_parse.add_opt("chh-cytosine", '\0', "Name of output file for cytosines in CHH", 
                          false , chhfile);
        // opt_parse.add_opt("window-size", 'w', "Window size", 
        //                   false , ws);
        // opt_parse.add_opt("pattern", 'p', "Pattern to search (default: CpG)", 
        //                   false , pattern);
        opt_parse.add_opt("index", 'I', "Output index of cytosines (default yes)", 
                          false, CYTOSINE_INDEX);
        opt_parse.add_opt("verbose", 'v', "print run info", 
                          false, VERBOSE);
        vector<string> leftover_args;
        opt_parse.parse(argc, argv, leftover_args);
        if (opt_parse.help_requested())
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

        if (!leftover_args.empty())
            infile = leftover_args.front();

        /**********************************************************************/

        if (VERBOSE)
            cerr << "Reading sequence file " << infile << endl;
        istream *in  = (infile.empty()) ?
            & std::cin : new std::ifstream(infile.c_str());
        string chrom = (infile.empty()) ?
            "chrU" : strip_path_and_suffix(infile);

        string sequences;
        char line[max_line_size];
        size_t n = 0;
        
        while (in->getline(line, max_line_size))
        {
            size_t ll = in->gcount() - 1;
            if (ll > 0 && line[0] != '>')
            {
                std::copy(line, line + ll, std::back_inserter(sequences));
                n += ll;
            }
        }
        
        for (size_t i = 0; i < sequences.size(); ++i)
            sequences[i] = toupper(sequences[i]);
        
        if (VERBOSE && !cfile.empty())
            cerr << "Extracting all cytosines" << endl;
        if (!cfile.empty())
            find_cytosines(sequences, cfile, chrom, CYTOSINE_INDEX);

        if (VERBOSE && !cgfile.empty())
            cerr << "Extracting cytosines in CpG context" << endl;
        if (!cgfile.empty())
            find_cytosines_CG(sequences, cgfile, chrom, CYTOSINE_INDEX);

        if (VERBOSE && !chgfile.empty())
            cerr << "Extracting cytosines in CHG context" << endl;
        if (!chgfile.empty())
            find_cytosines_CHG(sequences, chgfile, chrom, CYTOSINE_INDEX);

        if (VERBOSE && !chhfile.empty())
            cerr << "Extracting cytosines in CHH context" << endl;
        if (!chhfile.empty())
            find_cytosines_CHH(sequences, chhfile, chrom, CYTOSINE_INDEX);
        
        if (in != &std::cin) delete in;
    }
    catch (const RMAPException &e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
    catch (std::bad_alloc &ba) {
        cerr << "ERROR: could not allocate memory" << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
