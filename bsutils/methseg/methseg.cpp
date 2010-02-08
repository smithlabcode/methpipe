/* Copyright (C) 2009 University of Southern California
 *                    Andrew D Smith
 * Author: Andrew D. Smith
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "TwoStateHMM.hpp"

#include <numeric>
#include <cmath>
#include <fstream>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;
using std::max;
using std::min;
using std::pair;
using std::make_pair;

using std::ostream_iterator;
using std::ofstream;

// static string
// format_color(const GenomicRegion domain) {
//   const double score = domain.get_score();
//   string color;
//   if (score < 1.0)
//     color = "255,190,113";
//   else if (score < 10.0)
//     color = "255,163,55";
//   else if (score < 100.0)
//     color = "255,140,5";
//   else color = "184,100,0";
//   return "\t" + rmap::toa(domain.get_start()) + "\t" +
//     rmap::toa(domain.get_end()) + "\t" + color;
// }

// static void
// write_scores_bedgraph(const string &filename,
// 		      const vector<SimpleGenomicRegion> &cpgs,
// 		      const vector<double> &scores) {
//   std::ofstream wigout(filename.c_str());
//   for (size_t i = 0; i < cpgs.size(); ++i)
//     wigout << cpgs[i] << "\t" << scores[i] << "\n";
//   wigout.close();
// }

// static void
// build_domains(const bool VERBOSE, 
// 	      const vector<SimpleGenomicRegion> &cpgs,
// 	      const vector<double> &post_scores,
// 	      const vector<size_t> &reset_points,
// 	      const vector<bool> &classes,
// 	      vector<GenomicRegion> &domains) {
//   static const bool CLASS_ID = true;
//   size_t n_cpgs = 0, n_domains = 0, reset_idx = 1, prev_end = 0;
//   bool in_domain = false;
//   double score = 0;
//   for (size_t i = 0; i < classes.size(); ++i) {
//     if (reset_points[reset_idx] == i) {
//       if (in_domain) {
// 	in_domain = false;
// 	domains.back().set_end(prev_end);
// 	domains.back().set_score(score);
// 	n_cpgs = 0;
// 	score = 0;
//       }
//       ++reset_idx;
//     }
//     if (classes[i] == CLASS_ID) {
//       if (!in_domain) {
// 	in_domain = true;
// 	domains.push_back(GenomicRegion(cpgs[i]));
// 	domains.back().set_name("HYPO" + toa(n_domains++));
//       }
//       ++n_cpgs;
//       score += post_scores[i];
//     }
//     else if (in_domain) {
//       in_domain = false;
//       domains.back().set_end(prev_end);
//       domains.back().set_score(score);//n_cpgs);
//       n_cpgs = 0;
//       score = 0;
//     }
//     prev_end = cpgs[i].get_end();
//   }
//   // Do we miss the final domain?????
// }


// template <class T, class U> static void
// separate_regions(const bool VERBOSE, const size_t desert_size, 
// 		 vector<SimpleGenomicRegion> &cpgs,
// 		 vector<T> &meth, vector<U> &reads,
// 		 vector<size_t> &reset_points) {
//   if (VERBOSE)
//     cerr << "[SEPARATING BY CPG DESERT]" << endl;
//   // eliminate the zero-read cpgs
//   size_t j = 0;
//   for (size_t i = 0; i < cpgs.size(); ++i)
//     if (reads[i] > 0) {
//       cpgs[j] = cpgs[i];
//       meth[j] = meth[i];
//       reads[j] = reads[i];
//       ++j;
//     }
//   cpgs.erase(cpgs.begin() + j, cpgs.end());
//   meth.erase(meth.begin() + j, meth.end());
//   reads.erase(reads.begin() + j, reads.end());
  
//   // segregate cpgs
//   size_t prev_cpg = 0;
//   for (size_t i = 0; i < cpgs.size(); ++i) {
//     const size_t dist = (i > 0 && cpgs[i].same_chrom(cpgs[i - 1])) ? 
//       cpgs[i].get_start() - prev_cpg : numeric_limits<size_t>::max();
//     if (dist > desert_size)
//       reset_points.push_back(i);
//     prev_cpg = cpgs[i].get_start();
//   }
//   if (VERBOSE)
//     cerr << "CPGS RETAINED: " << cpgs.size() << endl
// 	 << "DESERTS REMOVED: " << reset_points.size() << endl << endl;
// }


// static void
// load_cpgs(const bool VERBOSE, string cpgs_file, 
// 	  vector<SimpleGenomicRegion> &cpgs,
// 	  // vector<double> &meth, vector<size_t> &reads) {
// 	  vector<pair<double, double> > &meth, vector<size_t> &reads) {
//   if (VERBOSE)
//     cerr << "[READING CPGS AND METH PROPS]" << endl;
//   vector<GenomicRegion> cpgs_in;
//   ReadBEDFile(cpgs_file, cpgs_in);
//   if (!check_sorted(cpgs_in))
//     throw RMAPException("CpGs not sorted in file \"" + cpgs_file + "\"");
  
//   for (size_t i = 0; i < cpgs_in.size(); ++i) {
//     cpgs.push_back(SimpleGenomicRegion(cpgs_in[i]));
//     meth.push_back(make_pair(cpgs_in[i].get_score(), 0.0));
//     const string r(cpgs_in[i].get_name());
//     reads.push_back(atoi(r.substr(r.find_first_of(":") + 1).c_str()));

//     meth.back().first = int(meth.back().first * reads.back());
//     meth.back().second = int(reads.back() - meth.back().first);
//   }

//   if (VERBOSE)
//     cerr << "TOTAL CPGS: " << cpgs.size() << endl
// 	 << "MEAN COVERAGE: " 
// 	 << accumulate(reads.begin(), reads.end(), 0.0)/reads.size() << endl
// 	 << endl;
// }


int
main(int argc, const char **argv) 
{
		cout << "This program is depreciated. Use hmr in the same directory instead. "
			 << "hmr add -F option, from which you can set False discovery rate"
			 << endl;
		return -1;

// 		try {

//     string outfile;
//     string scores_file;
//     string trans_file;
//     string dataset_name;
    
//     size_t desert_size = 2000;
//     size_t max_iterations = 10;
    
//     // run mode flags
//     bool USE_VITERBI = false;
//     bool VERBOSE = false;
//     bool BROWSER = false;
    
//     // corrections for small values (not parameters):
//     double tolerance = 1e-10;
//     double min_prob  = 1e-10;
    
//     /****************** COMMAND LINE OPTIONS ********************/
//     OptionParser opt_parse("methseg", "A program for segmenting DNA "
// 			   "methylation data"
// 			   "<cpg-BED-file>");
//     opt_parse.add_opt("out", 'o', "output file (BED format)", 
// 		      false, outfile);
//     opt_parse.add_opt("scores", 's', "scores file (WIG format)", 
// 		      false, scores_file);
//     opt_parse.add_opt("trans", 't', "trans file (WIG format)", 
// 		      false, trans_file);
//     opt_parse.add_opt("desert", 'd', "desert size", false, desert_size);
//     opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations); 
//     opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
//     opt_parse.add_opt("browser", 'B', "format for browser", false, BROWSER);
//     opt_parse.add_opt("name", 'N', "data set name", false, dataset_name);
//     opt_parse.add_opt("vit", 'V', "use Viterbi decoding (default: posterior)", 
// 		      false, USE_VITERBI);
//     vector<string> leftover_args;
//     opt_parse.parse(argc, argv, leftover_args);
//     if (argc == 1 || opt_parse.help_requested()) {
//       cerr << opt_parse.help_message() << endl;
//       return EXIT_SUCCESS;
//     }
//     if (opt_parse.about_requested()) {
//       cerr << opt_parse.about_message() << endl;
//       return EXIT_SUCCESS;
//     }
//     if (opt_parse.option_missing()) {
//       cerr << opt_parse.option_missing_message() << endl;
//       return EXIT_SUCCESS;
//     }
//     if (leftover_args.empty()) {
//       cerr << opt_parse.help_message() << endl;
//       return EXIT_SUCCESS;
//     }
//     const string cpgs_file = leftover_args.front();
//     /****************** END COMMAND LINE OPTIONS *****************/
    
//     // separate the regions by chrom and by desert
//     vector<SimpleGenomicRegion> cpgs;
//     // vector<double> meth;
//     vector<pair<double, double> > meth;
//     vector<size_t> reads;
//     load_cpgs(VERBOSE, cpgs_file, cpgs, meth, reads);
    
//     // separate the regions by chrom and by desert, and eliminate
//     // those isolated CpGs
//     vector<size_t> reset_points;
//     separate_regions(VERBOSE, desert_size, cpgs, meth, reads, reset_points);
    
//     vector<double> start_trans(2, 0.5), end_trans(2, 1e-10);
//     vector<vector<double> > trans(2, vector<double>(2, 0.25));
//     trans[0][0] = trans[1][1] = 0.75;
    
//     const double n_reads = accumulate(reads.begin(), reads.end(), 0.0)/reads.size();
//     double fg_alpha = 0.33*n_reads;
//     double fg_beta = 0.67*n_reads;
//     double bg_alpha = 0.67*n_reads;
//     double bg_beta = 0.33*n_reads;
    
//     const TwoStateHMMB hmm(min_prob, tolerance, max_iterations, VERBOSE);
//     hmm.BaumWelchTraining(meth, reset_points, start_trans, trans, 
// 			  end_trans, fg_alpha, fg_beta, bg_alpha, bg_beta);
    
//     /***********************************
//      * STEP 5: DECODE THE DOMAINS
//      */
//     vector<bool> classes;
//     vector<double> scores;
//     if (USE_VITERBI)
//       hmm.ViterbiDecoding(meth, reset_points, start_trans, trans, end_trans,
// 			  fg_alpha, fg_beta, bg_alpha, bg_beta, classes);
//     else hmm.PosteriorDecoding(meth, reset_points, start_trans, trans, end_trans,
// 			       fg_alpha, fg_beta, bg_alpha, bg_beta, classes, scores);
    
//     /***********************************
//      * STEP 6: WRITE THE RESULTS
//      */
//     if (VERBOSE)
//       cerr << "[COLLECTING POSTERIOR SCORES]" << endl;
//     vector<double> post_scores;
//     hmm.PosteriorScores(meth, reset_points, start_trans, trans, end_trans, 
// 			fg_alpha, fg_beta, bg_alpha, bg_beta, true, post_scores);
//     // if scores have been requested calculate and write them
//     if (!scores_file.empty())
//       write_scores_bedgraph(scores_file, cpgs, post_scores);
//     if (!trans_file.empty()) {
//       vector<double> fg_to_bg_scores;
//       hmm.TransitionPosteriors(meth, reset_points, start_trans, trans, end_trans, 
// 			       fg_alpha, fg_beta, bg_alpha, bg_beta, 
// 			       1, fg_to_bg_scores);
//       vector<double> bg_to_fg_scores;
//       hmm.TransitionPosteriors(meth, reset_points, start_trans, trans, end_trans, 
// 			       fg_alpha, fg_beta, bg_alpha, bg_beta, 
// 			       2, bg_to_fg_scores);
//       for (size_t i = 0; i < fg_to_bg_scores.size(); ++i)
// 	fg_to_bg_scores[i] = max(fg_to_bg_scores[i], bg_to_fg_scores[i]);
//       write_scores_bedgraph(trans_file, cpgs, fg_to_bg_scores);
//     }
    
//     vector<GenomicRegion> domains;
//     build_domains(VERBOSE, cpgs, post_scores,
// 		  reset_points, classes, domains);
    
//     std::ostream *out = (outfile.empty()) ? &cout : 
//       new std::ofstream(outfile.c_str());
//     if (BROWSER) {
//       if (dataset_name.empty())
// 	dataset_name = strip_path_and_suffix(cpgs_file);
//       *out << "track name=" << dataset_name << " itemRgb=\"On\"\n";
//     }
//     for (size_t i = 0; i < domains.size(); ++i) {
//       *out << domains[i];
//       if (BROWSER)
// 	*out << format_color(domains[i]);
//       *out << '\n';
//     }
//     if (out != &cout) delete out;
//   }
//   catch (RMAPException &e) {
//     cerr << "ERROR:\t" << e.what() << endl;
//     return EXIT_FAILURE;
//   }
//   catch (std::bad_alloc &ba) {
//     cerr << "ERROR: could not allocate memory" << endl;
//     return EXIT_FAILURE;
//   }
//   return EXIT_SUCCESS;
}
