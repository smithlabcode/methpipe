/*    fastLiftOver-CpG: lift over CpGs 
 *
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith
 *                       Song Qiang
 *
 *    Authors: Song Qiang
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
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <tr1/unordered_map>
#include <stdexcept>


#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "GenomicRegion.hpp"
#include "MethpipeFiles.hpp"


using std::string;
using std::ios_base;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::tr1::unordered_map;

struct GenomicSite
{
    string chrom;
    size_t pos;
    GenomicSite(const string &c = "", const size_t p = 0): chrom(c), pos(p) {}
};

static void
read_index_file(const string &indexFile,
                unordered_map<string, unordered_map<size_t, GenomicSite> > &index)
{
    std::ifstream in(indexFile.c_str());

    string toChrom;
    size_t toPos;
    size_t toEnd;
    string toName;
    while (in >> toChrom >> toPos >> toEnd >> toName)
    {
        const size_t dim = toName.find_first_of(":");
        const string chrom = toName.substr(0, dim);
        const size_t pos = atoi(toName.substr(dim + 1).c_str());
        index[chrom][pos] = GenomicSite(toChrom, toPos);
    }
}

int 
main(int argc, const char **argv) 
{
    try 
    {
        string indexfile;
        string tofile;
        string fromfile;

        bool VERBOSE = false;
    
        /****************** COMMAND LINE OPTIONS ********************/
        OptionParser opt_parse(strip_path(argv[0]), "Fast liftOver" );
        opt_parse.add_opt("indexfile", 'i', "index file", true, indexfile);
        opt_parse.add_opt("from", 'f', "Original file", true, fromfile);
        opt_parse.add_opt("to", 't', "Output file liftovered", true, tofile);
        opt_parse.add_opt("verbose", 'v', "print more information",
                          false, VERBOSE);

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
        //////////////////////////////////////////////////////////////
    
        unordered_map<string, unordered_map<size_t, GenomicSite> > index;
        if (VERBOSE)
            cerr << "Loading index file " << indexfile << endl;
        read_index_file(indexfile, index);
        
        std::ifstream from(fromfile.c_str());
        std::ofstream to(tofile.c_str());

        string chrom;
        size_t pos;
        string strand;
        string seq;
        double meth;
        size_t coverage;
        
        if (VERBOSE)
            cerr << "Lifting " << fromfile << " to " << tofile << endl;
        while (methpipe::read_site(from, chrom, pos, strand, seq, meth, coverage))
        {
            typedef unordered_map<string, unordered_map<size_t, GenomicSite> >::iterator HTItor;
            typedef unordered_map<size_t, GenomicSite>::iterator HTItor2;
            
            if ((ito = index.find(chrom)) != index.end()
                && (iti = ito->second.find(pos)) != ito->second.end())
            {
                    chrom = iti->second.chrom;
                    pos = iti->second.pos;
                    methpipe::write_site(to, chrom, pos, strand,
                                         seq, meth, coverage);
            }
        }
    }
    catch (const SMITHLABException &e) 
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
