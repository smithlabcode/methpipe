/* read_quality_prof.cpp 
 * Song Qiang <qiang.song@usc.edu> 2010
 * Home implementation of the same program in fastx_toolkit
 */

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "QualityScore.hpp"

using namespace std;

inline size_t 
effective_read_length(const string &seq,
                      const size_t NS,
                      const bool COUNT_N)
{
    if (COUNT_N)
        return seq.size() - NS;
    else
    {
        const size_t found = seq.find_last_not_of("N");
        return (found == string::npos) ? 0 : found + 1;
    }
}

int 
main(int argc, const char **argv) 
{
    try 
    {
        string infile;
        string outfile;

        /****************** COMMAND LINE OPTIONS ********************/
        OptionParser opt_parse(argv[0], "This program profile base composition"
                               " and quality score of reads"
                               "<fastq-reads-file>");
        opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
                          false , outfile);
        
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
        
        std::istream *in = infile.empty() ?
            &std::cin : new std::ifstream(infile.c_str());
        
        std::ostream *out = outfile.empty() ?
            &std::cout : new std::ofstream(outfile.c_str());
        
        string name_1, seq, name_2, scr;
        getline(*in, name_1);
        getline(*in, seq);
        getline(*in, name_2);
        getline(*in, scr);
        const size_t readlen = seq.size();
        
        vector<size_t>
            a_counts(readlen, 0), c_counts(readlen, 0),
            g_counts(readlen, 0), t_counts(readlen, 0), n_counts(readlen, 0);
        
        vector<double>
            min_scrs(readlen, 999), max_scrs(readlen, -999), sum_scrs(readlen, 0);
        
        size_t readnum = 0;
        while (in->good())
        {
            ++readnum;
            for (size_t i = 0; i < readlen; ++i)
            {
                switch (seq[i])
                {
                case 'A': ++a_counts[i]; break;
                case 'C': ++c_counts[i]; break;
                case 'G': ++g_counts[i]; break;
                case 'T': ++t_counts[i]; break;
                case 'N': ++n_counts[i]; break;
                }
                const char s = scr[i] - 64;
                if (s < min_scrs[i])
                    min_scrs[i] = s;
                if (s > max_scrs[i])
                    max_scrs[i] = s;
                sum_scrs[i] += s;
            }
            getline(*in, name_1);
            getline(*in, seq);
            getline(*in, name_2);
            getline(*in, scr);
        }
        
        *out << "column" << "\t" << "min" << "\t"
             << "max" << "\t" << "sum" << "\t"
             << "mean" << "\t" << "A-count" << "\t"
             << "C-count" << "\t" << "G-count" << "\t"
             << "T-count" << "\t" << "N-count" << endl;
        for (size_t i = 0; i < readlen; ++i)
            *out << i+1 << "\t" << min_scrs[i] << "\t"
                 << max_scrs[i] << "\t" << sum_scrs[i] << "\t"
                 << sum_scrs[i] / static_cast<double>(readnum) << "\t"
                 << a_counts[i] << "\t" << c_counts[i] << "\t"
                 << g_counts[i] << "\t" << t_counts[i] << "\t"
                 << n_counts[i] << endl;
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

