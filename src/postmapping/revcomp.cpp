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

void
revcomp(MappedRead &mr)
{
    if (mr.r.get_strand() == '+')
        mr.r.set_strand('-');
    else
        mr.r.set_strand('+');
    
    revcomp_inplace(mr.seq);
    std::reverse(mr.scr.begin(), mr.scr.end());
}

int 
main(int argc, const char **argv) 
{
    try 
    {
        
        string infile;
        string outfile;
        bool VERBOSE = false;

        /****************** COMMAND LINE OPTIONS ********************/
        OptionParser opt_parse(argv[0], "A program to revcomp A-rich reads");
        opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
        opt_parse.add_opt("output", 'o', "Name of output file", 
                          false, outfile);
        vector<string> leftover_args;
        opt_parse.parse(argc, argv, leftover_args);
        if (opt_parse.help_requested()) {
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

        /****************** END COMMAND LINE OPTIONS *****************/

        std::istream *in = infile.empty() ?
            &std::cin : new std::ifstream(infile.c_str());

        std::ostream *out = outfile.empty() ?
            &std::cout : new std::ofstream(outfile.c_str());
        
        MappedRead mr;
        bool read_is_good = true;

        try { *in >> mr; }
        catch (const RMAPException &e) { read_is_good = false;}
        
        while (read_is_good)
        {
            revcomp(mr);
            *out << mr << endl;
            
            try { *in >> mr; }
            catch (const RMAPException &e) { read_is_good = false;}
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
