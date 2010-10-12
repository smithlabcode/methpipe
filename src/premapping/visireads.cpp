/*    visigenes convert MR format reads to bed format highlight C's
 *    Song Qiang <qiang.song@usc.edu> 2010  
 */

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"

#include "bsutils.hpp"
#include "MappedRead.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::cin;

int 
main(int argc, const char ** argv)
{
    string pos_color = "0,0,255";
    string neg_color = "255,0,255";
    string infile  = "";
    string outfile = "";
    
    bool VERBOSE = false;
    bool TRIM_TRAINING_N = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], "a program for converting MR file to BED file"
                           "highlighting C's ",
                           ".mr file");
    opt_parse.add_opt("input", 'i', "Name of input file (default: stdin)", 
                      false, outfile);
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
                      false, outfile);
    opt_parse.add_opt("pos_color", '\0', "Color of reads on positive strand", 
                      false , pos_color);
    opt_parse.add_opt("neg_color", '\0', "Color of reads on negative strand", 
                      false , neg_color);
    opt_parse.add_opt("TRIM_TRAINING_N", 'n', "Trim trailing N's",
                      false, TRIM_TRAINING_N);
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

    if (leftover_args.size() >= 1)
        infile = leftover_args[0];
    if (leftover_args.size() >= 2)
        outfile = leftover_args[1];
    /****************** END COMMAND LINE OPTIONS *****************/

    std::istream *in =
        (infile.empty()) ? &cin : new std::ifstream(infile.c_str());

    std::ostream *out =
        (outfile.empty()) ? &cout : new std::ofstream(outfile.c_str());

    MappedRead mr;

    while (*in >> mr && in->good())
    {
        size_t start = 0;
        size_t end = mr.seq.size();
        if (TRIM_TRAINING_N)
            while (end > start && mr.seq[end -1] == 'N')
                --end;

        size_t nblock = 1;
        string block_starts = "0,";
        string block_sizes = "1,";
        for (size_t i = start + 1; i < end - 1; ++i)
            if (mr.seq[i] == 'C')
            {
                ++nblock;
                if (mr.r.get_strand() == '+')
                    block_starts += rmap::toa(i) + ",";
                else
                    block_starts += rmap::toa(end - 1 - i) + ",";
                block_sizes += rmap::toa(1) + ",";
            }

        ++nblock;
        block_starts += rmap::toa(end - 1) + ",";
        block_sizes += rmap::toa(1) + ",";
        string color(mr.r.get_strand() == '+' ? pos_color : neg_color);
        *out << mr.r << "\t"
             << mr.r.get_start() << "\t"
             << mr.r.get_end() << "\t"
             << color;
        if (nblock > 0)
            *out << "\t" << nblock << "\t"
                 << block_sizes << "\t"
                 << block_starts;
        *out << endl;
    }

//     if (!infile.empty()) in->close();
//     if (!outfile.empty()) out->close();
}
