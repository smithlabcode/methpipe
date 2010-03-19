/* dmrhmm 
 * Song Qiang <qiang.song@usc.edu> 2009
 *
 * This program reads output from methdiff, and then use three-state HMM to 
 * seek differentially methylated regions
 */

#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"

#include "HMM.hpp"

#include <algorithm>
#include <numeric>
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

void
get_meth_reads_number(const string &cpg_name,
					  size_t &meth_a, size_t &total_a,
					  size_t &meth_b, size_t &total_b)
{
		std::stringstream ss(cpg_name);
		char ch;
		
		ss.ignore(4);
		ss >> meth_a >> ch
		   >> total_a >> ch
		   >> meth_b >> ch
		   >> total_b;
		return;
}

size_t 
get_reads_num(const string &cpg_name)
{
		size_t meth_a, total_a, meth_b, total_b;
		get_meth_reads_number(cpg_name,
							  meth_a, total_a,
							  meth_b, total_b);
		return std::min(total_a - 2, total_b -2);
}


template <class T> string 
tostring(T t)
{
		std::ostringstream oss;
		oss << t;
		return oss.str();
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
		reset_points.push_back(cpgs.size());
}

// get significant CpG's and get intervals
class pred_significant_cpg
{
		double crit;
		bool LOWER_TAIL;
		bool UPPER_TAIL;
public:
		pred_significant_cpg(double c, int TAIL_MARKER) :
				crit(c), LOWER_TAIL(TAIL_MARKER & 1), UPPER_TAIL(TAIL_MARKER & 2) {}
		inline bool operator() (const GenomicRegion &cpg)
		{
				return (LOWER_TAIL && cpg.get_score() < crit) ||
						(UPPER_TAIL && cpg.get_score() > 1 - crit);
		}
};

void 
build_domains(const vector<GenomicRegion> &cpgs,
			  const vector<size_t> &reset_points,
			  const vector<state_type> &classes,
			  vector<GenomicRegion> &domains)
{
		assert(cpgs.size() == classes.size());
		
		const size_t STATE_DMR_LOWER = 0;
		const size_t STATE_NON_DMR = 1;
		const size_t STATE_DMR_UPPER = 2;

		for (size_t i = 0; i < reset_points.size() - 1; ++i)
		{
				size_t j = reset_points[i];
				GenomicRegion domain = cpgs[reset_points[i]];
				double total = domain.get_score();
				size_t cpg_num = 1;
				++j;
				while (j < reset_points[i+1])
				{
						if (classes[j] == classes[j - 1])
						{
								total += cpgs[j].get_score();
								++cpg_num;
						} else {
								// finish previous domain
								domain.set_end(cpgs[j - 1].get_end());
								domain.set_name(tostring(cpg_num)
												+ ":"
												+ tostring(total / cpg_num));
								domain.set_score(classes[j - 1]);
								if (classes[j - 1] == STATE_DMR_UPPER ||
									classes[j -1] == STATE_DMR_LOWER)
										domains.push_back(domain);
								
								// start a new domain
								domain = cpgs[j];
								total = cpgs[j].get_score();
								cpg_num = 1;
						}
						++j;
				}

				// finish final domain
				domain.set_end(cpgs[j - 1].get_end());
				domain.set_name(tostring(cpg_num)
								+ ":"
								+ tostring(total / cpg_num));
				domain.set_score(classes[j - 1]);
				if (classes[j - 1] != STATE_NON_DMR)
						domains.push_back(domain);
				
		}
}					


void 
build_domains(const vector<GenomicRegion> &cpgs,
			  const vector<size_t> &reset_points,
			  const vector<state_type> &classes,
			  const vector<double> &post_scores,
			  vector<GenomicRegion> &domains)
{
		assert(cpgs.size() == classes.size());
		assert(cpgs.size() == post_scores.size());
		
//		const size_t STATE_DMR_LOWER = 0;
		const size_t STATE_NON_DMR = 1;
//		const size_t STATE_DMR_UPPER = 2;

		vector<string> direction_labels;
		direction_labels.push_back("L");
		direction_labels.push_back("N");
		direction_labels.push_back("U");

		for (size_t i = 0; i < reset_points.size() - 1; ++i)
		{
				size_t j = reset_points[i];
				GenomicRegion domain = cpgs[j];
				double total = post_scores[j];
				size_t cpg_num = 1;
				++j;
				while (j < reset_points[i+1])
				{
						if (classes[j] == classes[j - 1])
						{
								total += post_scores[j];
								++cpg_num;
						} else {
								// finish previous domain
								if (classes[j - 1] != STATE_NON_DMR)
								{
										domain.set_end(cpgs[j - 1].get_end());
										domain.set_name("DMR:" +
														direction_labels[classes[j - 1]]);
										domain.set_score(total);
										domains.push_back(domain);
								}
								
								// start a new domain
								domain = cpgs[j];
								total = post_scores[j];
								cpg_num = 1;
						}
						++j;
				}

				// finish final domain
				if (classes[j - 1] != STATE_NON_DMR)
				{
						domain.set_end(cpgs[j - 1].get_end());
						domain.set_name("DMR:" + tostring(classes[j - 1]));
						domain.set_score(total);
						domains.push_back(domain);
				}
		}
}					


int
main(int argc, const char **argv) {

		try {

				string outfile;
				double crit = 0.05;
				int TAIL_MARKER = 1;
				size_t desert_size = 2000;

				// mode for HMM
				bool USE_VITERBI = false;
				
				// corrections for small values (not parameters):
				double tolerance = 1e-10;
				double min_prob  = 1e-10;
				size_t max_iterations = 10;

				bool VERBOSE = false;

				double fdr = 0.05;

				size_t BUFFER_SIZE = 100000;
    
				/****************** COMMAND LINE OPTIONS ********************/
				OptionParser opt_parse("dmrhmm", 
									   "A program for identifying DMRs (differentially "
									   "methylated regions) using three-state HMM ",
									   "<cpg_meth_diffs_file>");
				opt_parse.add_opt("out", 'o', "output file (BED format)", 
								  false, outfile);
				opt_parse.add_opt("desert", 'd', "size of desert",
								  false, desert_size);
				opt_parse.add_opt("buffer", 'B', "buffer size (in records, not bytes)", 
								  false , BUFFER_SIZE);
				opt_parse.add_opt("iteration", 'i',  "Max number of iteration for EM training (defualt 10)",
								  false, max_iterations);
				opt_parse.add_opt("viterbi", 'V', "useing viterbi decoding",
								  false, USE_VITERBI);
				opt_parse.add_opt("fdr", 'C', "Cutoff p-value (default 0.05)",
								  false, fdr);
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
						cerr << "Reading CpG file ... ";
				vector<GenomicRegion> cpgs;
				ReadBEDFile(cpgs_file, cpgs);
				if (!check_sorted(cpgs))
						throw RMAPException("file not sorted: \"" + cpgs_file + "\"");
				if (VERBOSE)
						cerr << "done" << endl;

				// seperate regions by chromosome and remove deserts
				if (VERBOSE)
						cerr << "Seperating regions ... ";
				vector<size_t> reset_points;
				seperate_regions(cpgs, reset_points, desert_size);
				if (VERBOSE)
						cerr << "done" << endl;

				// get the sequence of probabilities
				if (VERBOSE)
						cerr << "Extracting probabilities ... ";
				vector<value_type> prob_vals(cpgs.size());
				const double MIN_METH_PROB = 1e-2;
				for (size_t i = 0; i < cpgs.size(); ++i)
				{
						prob_vals[i] = cpgs[i].get_score();
						if (prob_vals[i] < MIN_METH_PROB)
						{
// 								if (VERBOSE)
// 										cerr << endl  << cpgs[i] << endl
// 											 << "Methylation probablity of "
// 											 << " is smaller than MIN_PROB "
// 											 << ". Adjusted to MIN_PROB" << endl;
								prob_vals[i] = MIN_METH_PROB;
						}
						if (prob_vals[i] > 1 - MIN_METH_PROB)
						{
// 								if (VERBOSE)
// 										cerr << endl << cpgs[i] << endl
// 											 << "Methylation probablity of "
// 											 << " is larger than MAX_PROB "
// 											 << ". Adjusted to MAX_PROB" << endl;
								prob_vals[i] = 1 - MIN_METH_PROB;
						}
				}
				
				if (VERBOSE)
						cerr << "done" << endl;
				
				/******************  HMM *************************************/

				// HMM setup
				if (VERBOSE)
						cerr << "Setting up HMM ... ";

				const size_t STATE_NUM = 3;
				const size_t STATE_DMR_LOWER = 0;
				const size_t STATE_NON_DMR = 1;
				const size_t STATE_DMR_UPPER = 2;

				vector<double> start_trans(STATE_NUM, 1.0 / STATE_NUM), end_trans(STATE_NUM, 1e-10);
				vector<vector<double> > trans(STATE_NUM, vector<double>(STATE_NUM, 0.2));
				trans[0][0] = trans[1][1] = trans[2][2] = 0.6;

				vector<distro_type> distros(STATE_NUM);
				distros[STATE_DMR_LOWER].set_alpha(0.8).set_beta(8);
				distros[STATE_NON_DMR].set_alpha(1.0).set_beta(1.0);
				distros[STATE_DMR_UPPER].set_beta(8).set_beta(0.2);

				const HMM hmm(min_prob, tolerance, max_iterations, VERBOSE);
				if (VERBOSE)
						cerr << "done" << endl;

				// EM training
				if (VERBOSE)
						cerr << "HMM: Baum-Welch Training ... ";
				hmm.BaumWelchTraining(prob_vals, reset_points,
									  start_trans, trans, end_trans,
									  distros);
				if (VERBOSE)
						cerr << "done" << endl;

				// Decoding: using either Viterbi or posterior
				if (VERBOSE)
						cerr << "HMM: Decoding ... ";
				vector<state_type> classes;
				vector<double> scores;
				if (USE_VITERBI)
						hmm.ViterbiDecoding(prob_vals, reset_points,
											start_trans, trans, end_trans,
											distros,
											classes);
				else 
						hmm.PosteriorDecoding(prob_vals, reset_points,
											  start_trans, trans, end_trans,
											  distros,
											  classes, scores);
				
				if (VERBOSE)
						cerr << "done" << endl;
				
				// Build domains
				if (VERBOSE)
						cerr << "HMM: Clustering segnificant domains ... ";
				vector<GenomicRegion> domains;
				build_domains(cpgs, reset_points, classes, scores, domains); 
				if (VERBOSE)
						cerr << "done" << endl;

				
				// output result
				if (VERBOSE)
						cerr << "HMM: output result ... ";
				std::ostream *out = (outfile.empty()) ? &cout : 
						new std::ofstream(outfile.c_str());
				copy(domains.begin(), domains.end(), 
					 ostream_iterator<GenomicRegion>(*out, "\n"));
				if (out != &cout) delete out;
				if (VERBOSE)
						cerr << "done" << endl;
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
