/* trim-adapter.cpp 
 * Song Qiang <qiang.song@usc.edu> 2010
 */

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include "clip_adaptor_from_reads.hpp"
#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"

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
        string adaptor;
        string readlen_file;

        const char LOWEST_SOLEXA_QUAL_SCORE = 59;
        bool COUNT_N = false;

        /****************** COMMAND LINE OPTIONS ********************/
        OptionParser opt_parse(argv[0], "This program trims adapter sequences from reads"
                               "<fast[a/q]-reads-file>");
        opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
                          false , outfile);
        opt_parse.add_opt("readlen", 'l',
                          "Name of output file of effective read lengths"  
                          "after triming adapters", 
                          false , readlen_file);
        opt_parse.add_opt("count-N", 'N', "Count trailing N's when report read length", 
                          false , COUNT_N);
        opt_parse.add_opt("adapter", 'a', "Adapter seqeunce", 
                          false , adaptor);
        
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
        
        std::ostream *out_readlen = readlen_file.empty() ?
            0 : new std::ofstream(readlen_file.c_str());

	const size_t MAX_READ_LEN = 1000; 
	vector<size_t> len_distr(MAX_READ_LEN + 1, 0)

        string name_1, seq, name_2, scr;
        while (in->good())
        {
            getline(*in, name_1);
            getline(*in, seq);
            getline(*in, name_2);
            getline(*in, scr);
            if (in->good())
            {
                const size_t NS =
                    clip_adaptor_from_read(adaptor,
                                           MIN_ADAPTOR_MATCH_SCORE, seq);
                fill_n(scr.rbegin(), NS, LOWEST_SOLEXA_QUAL_SCORE);
                *out << name_1 << endl
                     << seq << endl
                     << name_2 << endl
                     << scr << endl;
                if (out_readlen){
                   size_t cur_len = effective_read_length(seq, NS, COUNT_N);
		   if(cur_len < MAX_READ_LEN)
			len_distr[cur_len]++;
		   else
			len_distr[MAX_READ_LEN]++;
		}//if
            }
        }
	
	if(out_readlen){
		size_t i = MAX_READ_LEN;
		for( i = MAX_READ_LEN; i >=0 ; i--)
		{
			if(len_distr[i] > 0)
				break;
		}//for
		size_t stop_iter = i;
		*out_readlen << "EFFECTIVE_READ_LENGTH:\t" << "COUNT_OF_READS_WITH_THIS_LENGTH:" << endl;

		for(i = 0; i <= stop_iter; i++)
		   *out_readlen << i << len_distr[i] << endl; 
	}
        if (in != &std::cin) delete in;
        if (out != &std::cout) delete out;
        if (out_readlen != 0) delete out_readlen;
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

