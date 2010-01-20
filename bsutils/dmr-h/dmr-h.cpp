/* dmr-h 
 * Song Qiang <qiang.song@usc.edu> 2009
 *
 * dmr-h extends dmrfinder by allowing nonsignigicant CpGs insides DMR and
 * using HMM to determine the proportion of nonsignigicant CpGs inside and
 * outside DMRs. 
 *
 */

#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"

#include "TwoStateHMM.hpp"

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


// seperate CpGs by chromosome and remove deserts

size_t get_reads_num(const string &name)
{
		return atoi(name.substr(name.find(':') + 1).c_str());
}

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
get_sig_cpgs(const vector<GenomicRegion> &cpgs,
			 vector<size_t> &sig_cpg_labels,
			 const double crit,
			 const int TAIL_MARKER)
{
		pred_significant_cpg is_significant_cpg(crit, TAIL_MARKER);		

		for (size_t i = 0; i < cpgs.size(); ++i)
				sig_cpg_labels[i] = static_cast<size_t>(is_significant_cpg(cpgs[i]));
		return;
}


template <class T> string 
tostring(T t)
{
		std::ostringstream oss;
		oss << t;
		return oss.str();
}

void 
build_domains(const vector<GenomicRegion> &cpgs,
			  const vector<size_t> &reset_points,
			  const vector<bool> &classes,
			  vector<GenomicRegion> &domains,
			const size_t min_cpgs)
{
		assert(cpgs.size() == classes.size());
		
		const bool FG_CLASS = true;
		const bool BG_CLASS = false;

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
								domain.set_name("CpG:"
												+ tostring(cpg_num));
								domain.set_score( total / cpg_num);
								if (classes[j - 1] == FG_CLASS && cpg_num > min_cpgs)
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
				domain.set_name("CpG:"
								+ tostring(cpg_num));
				domain.set_score(static_cast<size_t>(1000 * total / cpg_num));
				if (classes[j - 1] == FG_CLASS && cpg_num > min_cpgs)
						domains.push_back(domain);
				
		}
}					

int
main(int argc, const char **argv) {

		try {

				string outfile;

				size_t min_cpgs = 20;

				double crit = 0.05;
				int TAIL_MARKER = 1;

				size_t desert_size = 2000;
				// mode for HMM
				bool USE_VITERBI = false;
				size_t max_iterations = 10;
				
				// corrections for small values (not parameters):
				double tolerance = 1e-10;
				double min_prob  = 1e-10;


				bool VERBOSE = false;

				/****************** COMMAND LINE OPTIONS ********************/
				OptionParser opt_parse("dmr-h", 
									   "A program for identifying DMRs (differentially "
									   "methylated regions) based on a file showing "
									   "probability of differential methylation at "
									   "each CpG or base. ",
									   "<cpg_meth_diffs_file>");
				opt_parse.add_opt("crit", 'c', "critical value (default: 0.05)", 
								  false, crit);
				opt_parse.add_opt("tail", 't',
								  "use which tail to determine signif.\n 1: LOWER TAIL; 2: UPPER_TAIL; 3: BOTH",
								  false, TAIL_MARKER);
				opt_parse.add_opt("out", 'o', "output file (BED format)", 
								  false, outfile);
				opt_parse.add_opt("desert", 'd', "size of desert",
								  false, desert_size);
				opt_parse.add_opt("width", 'w', "width in terms of CpGs of min DMR size (default: 20)", 
								  false, min_cpgs);
				opt_parse.add_opt("iteration", 'i',  "Max number of iteration for EM training (defualt 10)",
								  false, max_iterations);
				opt_parse.add_opt("viterbi", 'V', "useing viterbi decoding",
								  false, USE_VITERBI);
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

				// discretize p-values and get the 0/1 sequence
				vector<size_t> labels(cpgs.size());
				if (VERBOSE)
						cerr << "Discretizing p-values ... ";
				get_sig_cpgs(cpgs, labels, crit, TAIL_MARKER);
				if (VERBOSE)
						cerr << "done" << endl;

				
				/******************  HMM *************************************/

				// HMM setup
				vector<double> start_trans(2, 0.5), end_trans(2, 1e-10);
				vector<vector<double> > trans(2, vector<double>(2, 0.25));
				trans[0][0] = trans[1][1] = 0.75;

				double fg_prob = 0.8;
				double bg_prob = 0.2;

				const TwoStateHMM hmm(min_prob, tolerance, max_iterations, VERBOSE);

				// EM training
				if (VERBOSE)
						cerr << "HMM: Baum-Welch Training ... ";
				hmm.BaumWelchTraining(labels, reset_points,
									  start_trans, trans, end_trans,
									  fg_prob, bg_prob);
				if (VERBOSE)
						cerr << "done" << endl;

				// Decoding: using either Viterbi or posterior
				if (VERBOSE)
						cerr << "HMM: Decoding ... ";
				vector<bool> classes;
				vector<double> scores;
				if (USE_VITERBI)
						hmm.ViterbiDecoding(labels, reset_points,
											start_trans, trans, end_trans,
											fg_prob, bg_prob,
											classes);
				else 
						hmm.PosteriorDecoding(labels, reset_points,
											  start_trans, trans, end_trans,
											  fg_prob, bg_prob,
											  classes, scores);
				if (VERBOSE)
						cerr << "done" << endl;
				
				// Build domains
				if (VERBOSE)
						cerr << "HMM: clustering segnificant intervals ... ";
				vector<GenomicRegion> domains;
				build_domains(cpgs, reset_points, classes, domains, min_cpgs); 
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
