/* diffseg 
 * Song Qiang <qiang.song@usc.edu> 2009
 *
 * This program reads output from methdiff, first convert it to a 0/1 sequence
 * and then get a sequence of intervals between 1's in basepairs. It then
 * segments this interval sequence using HMM method. We use two-state HMM with
 * geometric distribution.
 *
 */

#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "StringTool"

#include "FileIterator.hpp"

#include <gsl/gsl_sf_gamma.h>


#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>


using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;

using std::numeric_limits;
using std::ostream_iterator;
using std::ofstream;

// get reads number from CpG name which is in the form CpG:integer
int
get_reads_num(const string &name)
{
		return atoi(name.substr(name.find(':') + 1));
}

// seperate CpGs by chromosome and remove deserts

void
seperate_regions(vector<GenomicRegion> &cpgs,
				 vector<size_t> &reset_points,
				 const size_t desert_size)
{
		size_t j = 0;
		for (size_t i = 0; i < cpgs.size(); ++i)
				if (get_reads_num(cpgs[i].get_name()) > 0)
				{
						cpgs[j] = cpgs[i];
						++j;
				}
		cpgs.erase(cpgs.begin() + j, cpgs.end());
		
		reset_points.push_back(0);
		for (size_t i = 1; i < cpgs.size(); ++i)
		{
				const size_t dist = cpgs[i].distance(cpgs[i - 1]);
				if (dist > desert_size)
						reset_points.push_back(i);
		}
}

// get significant CpG's and get intervals
class pred_significant_cpg
{
		double crit;
		bool LOWER_TAIL;
		bool UPPER_TAIL;
public:
		pred_significant_cpg(double c, short si) :
				crit(c), LOWER_TAIL(si & 1), UPPER_TAIL(si & 2) {}
		inline bool operator() (const GenomicRegion &cpg)
		{
				return (LOWER_TAIL && cpg.get_score() < crit) ||
						(UPPER_TAIL && cpg.get_score() > 1 - crit);
		}
}

void
get_intervals(const vector<GenomicRegion> &cpgs,
			  vector<size_t> &reset_points,
			  vector<GenomicRegion> &sig_cpgs,
			  vector<SimpleGenomicRegion> &intervals,
			  const double crit,
			  const short TAIL_MARKER)
{
		for (size_t i = 0; i < sig_cpgs.size(); i++)
				sig_cpgs[i].set_score(0);
		const pred_significant_cpg is_significant_cpg(crit, TAIL_MARKER);		

		// deal with the end of CpG sequence
		reset_points.push_back(cpgs.begin() - cpgs.end());
		
		size_t j = 0;			// counter for reset_points of intervals
		for (size_t i = 0; i < reset_points.size() - 1; ++i)
		{
				bool is_first_interval = true; 
				const vector<GenomicRegion>::const_iterator first
						= cpgs.begin() + reset_points[i];
				const vector<GenomicRegion>::const_iterator last
						= cpgs.begin() + reset_points[i + 1];
				vector<GenomicRegion>::const_iterator iter
						= std::find_if(first, last, is_significant_cpg);
				while (iter != last)
				{
						sig_cpgs[iter - cpgs.begin()].set_score(1);
						const GenomicRegion last_sig_cpg = *iter;
						iter = (iter + 1 == last) ?
								last : std::find_if(iter + 1, last, is_significant_cpg);
						if (iter != last && iter->same_chrom(last_sig_cpg))
						{
								intervals.push_back(
										SimpleGenomicRegion(last_sig_cpg.get_chrom(),
															last_sig_cpg.get_end(),
															iter->get_start()));
								if (is_first_interval)
								{
										reset_points[j] = intervals.size() - 1;
										is_first_interval = false;
										++j;
								}
						}
				}
		}
		reset_points.erase(reset_points.begin() + j, reset_points.end());
}

void 
build_domains(const vector<SimpleGenomicRegion> &intervals,
			  const vector<size_t> &reset_points,
			  const vector<bool> &classes,
			  vector<SimpleGenomicRegion> &domains)
{
		assert(intervals.size() == classes.size());
		
		const bool FG_CLASS = true;
		const bool BG_CLASS = false;
		for (size_t i = 0; i < reset_points.size() - 1; ++i)
		{
				const vector<bool>::const_iterator first
						= classes.begin() + reset_points[i];
				const vector<bool>::const_iterator last
						= classes.begin() + reset_points[i + 1];
				vector<bool>::const_iterator iter
						=  std::find(first, last, FG_CLASS);
				while (iter != last)
				{
						const size_t start_index = iter - classes.begin();
						iter = (iter + 1 == last) ?
								last : std::find(iter + 1, last, BG_CLASS);
						const size_t end_index = iter - classes.begin() - 1;
						domains.push_back(
								SimpleGenomicRegion(intervals[start_index].get_chrom(),
													intervals[start_index].get_start(),
													intervals[end_index].get_end()));
						iter = (iter == last) ? last : std::find(iter, last, FG_CLASS);
				}
		}
}					

int
main(int argc, const char **argv) {

		try {

				string outfile;
				double crit = 0.05;
				short TAIL_MARKER = 1;
				size_t deserts = 2000;

				// mode for HMM
				bool USE_VITERBI = false;
				bool VERBOSE = false;
				bool BROWSER = false;
				
				// corrections for small values (not parameters):
				double tolerance = 1e-10;
				double min_prob  = 1e-10;

				bool VERBOSE = false;

				size_t BUFFER_SIZE = 100000;
    
				/****************** COMMAND LINE OPTIONS ********************/
				OptionParser opt_parse("diffseg", 
									   "A program for identifying DMRs (differentially "
									   "methylated regions) based on a file showing "
									   "probability of differential methylation at "
									   "each CpG or base. ",
									   "<cpg_meth_diffs_file>");
				opt_parse.add_opt("crit", 'c', "critical value (default: 0.05)", 
								  false, crit);
				opt_parse.add_opt("tail", 't',
								  "use which tail to determine signif.\n" + 
								  "1: LOWER TAIL; 2: UPPER_TAIL; 3: BOTH",
								  false, TAIL_MARKER);
				opt_parse.add_opt("out", 'o', "output file (BED format)", 
								  false, outfile);
				opt_parse.add_opt("buffer", 'B', "buffer size (in records, not bytes)", 
								  false , BUFFER_SIZE);
				opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
				vector<string> leftover_args;
				opt_parse.parse(argc, argv, leftover_args);
				if (argc == 1 || opt_parse.help_requested()) {
						cerr << opt_parse.help_message() << endl
							 << opt_parse.about_message() << endl;
						return EXIT_SUCCESS;
				}
				if (opt_parse.about_requested()) {
						cerr << opt_parse.about_message() << endl;
						return EXIT_SUCCESS;
				}
				if (opt_parse.option_missing()) {
						cerr << opt_parse.option_missing_message() << endl;
						return EXIT_SUCCESS;
				}
				if (leftover_args.size() != 1) {
						cerr << opt_parse.help_message() << endl;
						return EXIT_SUCCESS;
				}
				const string cpgs_file = leftover_args.front();
				/****************** END COMMAND LINE OPTIONS *****************/

				if (VERBOSE)
						cerr << "[READING CPGS]";
				vector<GenomicRegion> cpgs;
				ReadBEDFile(cpgs_file, cpgs);
				if (!check_sorted(cpgs))
						throw RMAPException("file not sorted: \"" + cpgs_file + "\"");
				if (VERBOSE)
						cerr << "[DONE]" << endl;

				// seperate regions by chromosome and remove deserts
				vector<size_t> reset_points;
				seperate_regions(cpgs, reset_points, desert_size);

				// discretize p-values and get intervals
				vector<GenomicRegion> sig_cpgs(cpgs);
				vector<SimpleGenomicRegion> intervals;
				vector<size_t> widths(intervals.size());
				if (VERBOSE)
						cerr << "[Discretizing p-values ...]"
				get_intervals(cpgs, reset_points, sig_cpgs, intervals, crit, TAIL_MARKER);
				for (size_t i = 0; i < intervals.size(); ++i)
						widths[i] = intervals[i].get_width();
				if (VERBOSE)
						cerr << "[done]" << endl;

				
				/******************  HMM *************************************/

				// HMM setup
				vector<double> start_trans(2, 0.5), end_trans(2, 1e-10);
				vector<vector<double> > trans(2, vector<double>(2, 0.25));
				trans[0][0] = trans[1][1] = 0.75;

				const double mean_width
						= std::accumulate(widths.begin(), widths.end(), 0.0) / widths.size();
				double fg_lambda = 1 / (0.25 * mean_width);
				double bg_lambda = 1 / (4 * mean_width);

				const TwoStateHMM hmm(min_prob, tolerance, max_iterations, VERBOSE);

				// EM training
				hmm.BaumWelchTraining(widths, reset_points,
									  start_trans, trans, end_trans,
									  fg_lambda, bg_lambda);

				// Decoding: using either Viterbi or posterior
				vector<bool> classes;
				vector<double> scores;
				if (USE_VITERBI)
						hmm.ViterbiDecoding(widths, reset_points,
											start_trans, trans, end_trans,
											fg_lambda, bg_lambda,
											classes);
				else 
						hmm.PosteriorDecoding(widths, reset_points,
											  start_trans, trans, end_trans,
											  fg_lambda, bg_lambda,
											  classes, scores);
				
				// Build domains
				vector<SimpleGenomicRegion> domains;
				build_domains(intervals, reset_points, classes, domains); 

				// output result
				std::ostream *out = (outfile.empty()) ? &cout : 
						new std::ofstream(outfile.c_str());
				copy(domains.begin(), domains.end(), 
					 ostream_iterator<SimpleGenomicRegion>(*out, "\n"));
				if (out != &cout) delete out;
		}
		catch (RMAPException &e) {
				cerr << "ERROR:\t" << e.what() << endl;
				return EXIT_FAILURE;
		}
		catch (std::bad_alloc &ba) {
				cerr << "ERROR: could not allocate memory" << endl;
				return EXIT_FAILURE;
		}
		return EXIT_SUCCESS;
}
