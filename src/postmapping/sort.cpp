/*    sortmapped: a program for sorting mapped read files
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

#include <fstream>
#include <iterator>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "GenomicRegion.hpp"
#include "MappedRead.hpp"

using std::string;
using std::vector;
using std::ostream;
using std::istream;
using std::endl;
using std::cerr;
using std::pair;
using std::make_pair;
using std::sort;

struct genome_location_cmp 
{
    bool operator()(const MappedRead &a, 
                    const MappedRead &b) const 
    {
        const GenomicRegion &lhs  = a.r;
        const GenomicRegion &rhs  = b.r;
        return (lhs.get_chrom() < rhs.get_chrom()) ||
            (lhs.get_chrom() == rhs.get_chrom() &&
             lhs.get_start() < rhs.get_start()) || 
            (lhs.get_chrom() == rhs.get_chrom() &&
             lhs.get_start() == rhs.get_start() &&
             lhs.get_strand() < rhs.get_strand());
    }
};

struct read_name_cmp 
{
    bool operator()(const MappedRead &a, 
                    const MappedRead &b) const 
    {
        return a.r.get_name() < b.r.get_name();
    }
};


int main(int argc, const char **argv) 
{

    try 
    {
        /* FILES */
        string outfile;
        bool SORT_ON_NAME = false;

        /****************** GET COMMAND LINE ARGUMENTS ***************************/
        OptionParser opt_parse(argv[0], "A program for sorting mapped read files",
                               "<mapped-reads>");
        opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
                          false , outfile);
        opt_parse.add_opt("name", 'N', "Sort by the names of reads", 
                          false , SORT_ON_NAME);
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
        const string input_file_name =
            leftover_args.size() ? leftover_args.front() : "";
        /**********************************************************************/
        vector<MappedRead> mrs;
        std::istream *in = input_file_name.empty() ?
            &std::cin : new std::ifstream(input_file_name.c_str());
        /*std::copy(std::istream_iterator<MappedRead>(*in),
                  std::istream_iterator<MappedRead>(),
                  std::back_inserter(mrs));*/
	MappedRead tmp;
	while(in->good() && (*in >> tmp)){
		mrs.push_back(tmp);
	}
        if (in != &std::cin) delete in;
        
        if (SORT_ON_NAME)
            std::sort(mrs.begin(), mrs.end(), read_name_cmp());
        else
            std::sort(mrs.begin(), mrs.end(), genome_location_cmp());
        
        ostream* out = (outfile.empty()) ? 
            &std::cout : new std::ofstream(outfile.c_str());
        std::copy(mrs.begin(), mrs.end(),
                  std::ostream_iterator<MappedRead>(*out, "\n"));
        if (out != &std::cout) delete out;
    }
    catch (RMAPException &e) 
    {
        cerr << "ERROR:\t" << e.what() << endl;
        return EXIT_FAILURE;
    }
    catch (std::bad_alloc &ba) 
    {
        cerr << "ERROR: could not allocate memory" << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
