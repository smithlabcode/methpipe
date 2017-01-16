/* Copyright (C) 2014 University of Southern California
 *                         Andrew D Smith
 * Author: Andrew D. Smith, Song Qiang, Jenny Qu
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

#include <numeric>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <unistd.h>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "TwoStateHMM.hpp"
#include "MethpipeFiles.hpp"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;
using std::max;
using std::min;
using std::pair;


double
get_fdr_cutoff(const vector<double> &scores, const double fdr) {
  if (fdr <= 0)
    return numeric_limits<double>::max();
  else if (fdr > 1)
    return numeric_limits<double>::min();
  vector<double> local(scores);
  std::sort(local.begin(), local.end());
  size_t i = 0;
  for (; i < local.size() - 1 &&
	 local[i+1] < fdr*static_cast<double>(i+1)/local.size(); ++i);
  return local[i];
}


static void
get_domain_scores_rep(const vector<bool> &classes,
		      const vector<vector<pair<double, double> > > &meth,
		      const vector<size_t> &reset_points,
		      vector<double> &scores) {
  size_t NREP = meth.size();
  static const bool CLASS_ID = true;
  size_t reset_idx = 1;
  bool in_domain = false;
  double score = 0;
  for (size_t i = 0; i < classes.size(); ++i) {
    if (reset_points[reset_idx] == i) {
      if (in_domain) {
	in_domain = false;
	scores.push_back(score);
	score = 0;
      }
      ++reset_idx;
    }
    if (classes[i] == CLASS_ID) {
      in_domain = true;
      for (size_t r = 0; r < NREP ; ++r) {
	if (meth[r][i].first + meth[r][i].second >= 1) {
	  score += 1.0 - (meth[r][i].first/(meth[r][i].first +
                                            meth[r][i].second));
	}
      }
    } else if (in_domain) {
      in_domain = false;
      scores.push_back(score);
      score = 0;
    }
  }
}


static void
build_domains(const bool VERBOSE,
	      const vector<SimpleGenomicRegion> &cpgs,
	      const vector<double> &post_scores,
	      const vector<size_t> &reset_points,
	      const vector<bool> &classes,
	      vector<GenomicRegion> &domains) {
  static const bool CLASS_ID = true;
  size_t n_cpgs = 0, n_domains = 0, reset_idx = 1, prev_end = 0;
  bool in_domain = false;
  double score = 0;
  for (size_t i = 0; i < classes.size(); ++i) {
    if (reset_points[reset_idx] == i) {
      if (in_domain) {
	in_domain = false;
	domains.back().set_end(prev_end);
	domains.back().set_score(n_cpgs);
	n_cpgs = 0;
	score = 0;
      }
      ++reset_idx;
    }
    if (classes[i] == CLASS_ID) {
      if (!in_domain) {
	in_domain = true;
	domains.push_back(GenomicRegion(cpgs[i]));
	domains.back().set_name("HYPO" + toa(n_domains++));
      }
      ++n_cpgs;
      score += post_scores[i];
    } else if (in_domain) {
      in_domain = false;
      domains.back().set_end(prev_end);
      domains.back().set_score(n_cpgs);
      n_cpgs = 0;
      score = 0;
    }
    prev_end = cpgs[i].get_end();
  }
}


//Modified to take multiple replicates
template <class T, class U> static void
separate_regions(const bool VERBOSE, const size_t desert_size,
		 vector<vector<SimpleGenomicRegion> > &cpgs,
		 vector<vector<T> > &meth, vector<vector<U> > &reads,
		 vector<size_t> &reset_points) {
  if (VERBOSE)
    cerr << "[SEPARATING BY CPG DESERT]" << endl;
  // eliminate the zero-read cpg sites if no coverage in any replicates
  size_t NREP = cpgs.size();
  size_t j = 0;
  bool all_empty;
  size_t rep_idx;
  for (size_t i = 0; i < cpgs[0].size(); ++i) {
    all_empty = true;
    rep_idx = 0;
    while (all_empty && rep_idx < NREP) {
      if (reads[rep_idx][i]==0) ++rep_idx;
      else all_empty = false;
    }
    if (!all_empty) {
      for (size_t r = 0; r < NREP; ++r) {
	cpgs[r][j] = cpgs[r][i];
	meth[r][j] = meth[r][i];
	reads[r][j] = reads[r][i];
      }
      ++j;
    }
  }
  for (size_t r = 0; r < NREP; ++r) {
    cpgs[r].erase(cpgs[r].begin() + j, cpgs[r].end());
    meth[r].erase(meth[r].begin() + j, meth[r].end());
    reads[r].erase(reads[r].begin() + j, reads[r].end());
  }

  // segregate cpgs
  size_t prev_cpg = 0;
  for (size_t i = 0; i < cpgs[0].size(); ++i) {
    const size_t dist = (i > 0 && cpgs[0][i].same_chrom(cpgs[0][i - 1])) ?
      cpgs[0][i].get_start() - prev_cpg : numeric_limits<size_t>::max();
    if (dist > desert_size)
      reset_points.push_back(i);
    prev_cpg = cpgs[0][i].get_start();
  }
  reset_points.push_back(cpgs[0].size());
  if (VERBOSE)
    cerr << "CPGS RETAINED: " << cpgs[0].size() << endl
	 << "DESERTS REMOVED: " << reset_points.size() - 2 << endl << endl;
}

static void
shuffle_cpgs_rep(const TwoStateHMMB &hmm,
		 vector<vector<pair<double, double> > > meth,
		 vector<size_t> reset_points,
		 const vector<double> &start_trans,
		 const vector<vector<double> > &trans,
		 const vector<double> &end_trans,
		 const vector<double> fg_alpha, const vector<double> fg_beta,
		 const vector<double> bg_alpha, const vector<double> bg_beta,
		 vector<double> &domain_scores) {

  srand(time(0) + getpid());
  size_t NREP = meth.size();

  for (size_t r =0 ; r < NREP; ++r)
    random_shuffle(meth[r].begin(), meth[r].end());

  vector<bool> classes;
  vector<double> scores;
  hmm.PosteriorDecoding_rep(meth, reset_points, start_trans, trans,
			    end_trans, fg_alpha, fg_beta, bg_alpha,
			    bg_beta, classes, scores);
  get_domain_scores_rep(classes, meth, reset_points, domain_scores);
  sort(domain_scores.begin(), domain_scores.end());
}

static void
assign_p_values(const vector<double> &random_scores,
		const vector<double> &observed_scores,
		vector<double> &p_values) {
  const double n_randoms =
    random_scores.size() == 0 ? 1 : random_scores.size();
  for (size_t i = 0; i < observed_scores.size(); ++i)
    p_values.push_back((random_scores.end() -
			upper_bound(random_scores.begin(),
				    random_scores.end(),
				    observed_scores[i]))/n_randoms);
}


static void
read_params_file(const bool VERBOSE,
		 const string &params_file,
		 double &fg_alpha,
		 double &fg_beta,
		 double &bg_alpha,
		 double &bg_beta,
		 vector<double> &start_trans,
		 vector<vector<double> > &trans,
		 vector<double> &end_trans,
		 double &fdr_cutoff) {
  string jnk;
  std::ifstream in(params_file.c_str());

  in >> jnk >> fg_alpha
     >> jnk >> fg_beta
     >> jnk >> bg_alpha
     >> jnk >> bg_beta
     >> jnk >> start_trans[0]
     >> jnk >> start_trans[1]
     >> jnk >> trans[0][0]
     >> jnk >> trans[0][1]
     >> jnk >> trans[1][0]
     >> jnk >> trans[1][1]
     >> jnk >> end_trans[0]
     >> jnk >> end_trans[1]
     >> jnk >> fdr_cutoff;

  if (VERBOSE)
    cerr << "Read in params from " << params_file << endl
	 << "FG_ALPHA\t" << fg_alpha << endl
	 << "FG_BETA\t" << fg_beta << endl
	 << "BG_ALPHA\t" << bg_alpha << endl
	 << "BG_BETA\t" << bg_beta << endl
	 << "S_F\t" << start_trans[0] << endl
	 << "S_B\t" << start_trans[1] << endl
	 << "F_F\t" << trans[0][0] << endl
	 << "F_B\t" << trans[0][1] << endl
	 << "B_F\t" << trans[1][0] << endl
	 << "B_B\t" << trans[1][1] << endl
	 << "F_E\t" << end_trans[0] << endl
	 << "B_E\t" << end_trans[1] << endl
	 << "FDR_CUTOFF\t" << fdr_cutoff << endl;
}


static void
write_params_file(const string &outfile,
		  const vector<double> fg_alpha,
		  const vector<double> fg_beta,
		  const vector<double> bg_alpha,
		  const vector<double> bg_beta,
		  const vector<double> &start_trans,
		  const vector<vector<double> > &trans,
		  const vector<double> &end_trans) {

  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out.precision(30);
  for (size_t r =0; r < fg_alpha.size(); ++r)
    out << "FG_ALPHA_" << r+1 << "\t" << std::setw(14) << fg_alpha[r] << "\t"
	<< "FG_BETA_" << r+1 << "\t" << std::setw(14) << fg_beta[r] << "\t"
	<< "BG_ALPHA_" << r+1 << "\t" << std::setw(14) << bg_alpha[r] << "\t"
	<< "BG_BETA_" << r+1 << "\t" << std::setw(14) << bg_beta[r] << endl;

  out << "S_F\t" << start_trans[0] << endl
      << "S_B\t" << start_trans[1] << endl
      << "F_F\t" << trans[0][0] << endl
      << "F_B\t" << trans[0][1] << endl
      << "B_F\t" << trans[1][0] << endl
      << "B_B\t" << trans[1][1] << endl
      << "F_E\t" << end_trans[0] << endl
      << "B_E\t" << end_trans[1] << endl;
}



int
main(int argc, const char **argv) {

  try {

    const char* sep = ",";
    string outfile;

    bool DEBUG = false;
    size_t desert_size = 1000;
    size_t max_iterations = 10;

    // run mode flags
    bool VERBOSE = false;

    // corrections for small values (not parameters):
    double tolerance = 1e-10;
    double min_prob  = 1e-10;

    string params_in_files;
    string params_out_file;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "Program for identifying "
			   "HMRs from methylomes of replicates ",
			   "<methcount-files separated by comma>");
    opt_parse.add_opt("out", 'o', "output file (default: stdout)",
		      false, outfile);
    opt_parse.add_opt("desert", 'd', "max dist btwn cpgs with reads in HMR",
		      false, desert_size);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("debug", 'D', "print more run info", false, DEBUG);
    opt_parse.add_opt("params-in", 'P', "HMM parameter files for "
                      "individual methylomes (separated with comma)",
		      false, params_in_files);
    opt_parse.add_opt("params-out", 'p', "write HMM parameters to this file",
		      false, params_out_file);

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
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string cpgs_files = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    vector<string> cpgs_file =  smithlab::split(cpgs_files, sep, false);
    vector<string> params_in_file;
    if(!params_in_file.empty()) {
      params_in_file = smithlab::split(params_in_files, sep, false);
      assert(cpgs_file.size() == params_in_file.size());
    }

    size_t NREP = cpgs_file.size();

    // separate the regions by chrom and by desert
    vector<vector<SimpleGenomicRegion> > cpgs;
    vector<vector<pair<double, double> > > meth;
    vector<vector<size_t> > reads;

    for (size_t i = 0; i < NREP; ++i) {
      vector<SimpleGenomicRegion> cpgs_rep;
      vector<pair<double, double> > meth_rep;
      vector<size_t> reads_rep;
      if (VERBOSE)
	cerr << "[READING CPGS AND METH PROPS] from " << cpgs_file[i] << endl;

      methpipe::load_cpgs(cpgs_file[i], cpgs_rep, meth_rep, reads_rep);

      cpgs.push_back(cpgs_rep);
      meth.push_back(meth_rep);
      reads.push_back(reads_rep);
      if (VERBOSE)
        cerr << "TOTAL CPGS: " << cpgs[i].size() << endl
             << "MEAN COVERAGE: "
             << accumulate(reads[i].begin(), reads[i].end(),
                           0.0)/reads[i].size()
             << endl;
    }

    bool aligned = true;
    size_t rep_idx = 0;
    for (; aligned && rep_idx < NREP; ++rep_idx) {
      aligned = (cpgs[rep_idx].size() == cpgs[0].size());
    }
    if (rep_idx < NREP)
      throw SMITHLABException("Inputs contain different number of sites");

    // separate the regions by chrom and by desert, and eliminate
    // those isolated CpGs
    vector<size_t> reset_points;
    separate_regions(VERBOSE, desert_size, cpgs, meth, reads, reset_points);

    /****************** Read in params *****************/
    vector<double> start_trans(2, 0.5), end_trans(2, 1e-10);
    vector<vector<double> > trans(2, vector<double>(2, 0.25));
    trans[0][0] = trans[1][1] = 0.75;
    const TwoStateHMMB hmm(min_prob, tolerance, max_iterations, VERBOSE, DEBUG);
    vector<double> reps_fg_alpha(NREP, 0);
    vector<double> reps_fg_beta(NREP, 0);
    vector<double> reps_bg_alpha(NREP, 0);
    vector<double> reps_bg_beta(NREP, 0);
    double fdr_cutoff = std::numeric_limits<double>::max();

    if (!params_in_file.empty()) {
      // READ THE PARAMETERS FILES
      double fdr_cutoff_rep;
      for (size_t i= 0; i< NREP; ++i) {
	read_params_file(VERBOSE, params_in_file[i], reps_fg_alpha[i],
                         reps_fg_beta[i], reps_bg_alpha[i], reps_bg_beta[i],
			 start_trans, trans, end_trans, fdr_cutoff_rep);
      }
    } else {
      for (size_t r = 0; r < NREP; ++r) {
	// Actually, there are many 0s in reads[r],
	// But the parameter start points don't need to be accurate
	double n_reads =
	  accumulate(reads[r].begin(), reads[r].end(), 0.0)/reads[r].size();
	reps_fg_alpha[r] = 0.33*n_reads;
	reps_fg_beta[r] = 0.67*n_reads;
	reps_bg_alpha[r] = 0.67*n_reads;
	reps_bg_beta[r] = 0.33*n_reads;
      }
    }

    if (max_iterations > 0)
      hmm.BaumWelchTraining_rep(meth, reset_points,
    				start_trans, trans, end_trans,
    				reps_fg_alpha, reps_fg_beta,
    				reps_bg_alpha, reps_bg_beta);

    if (!params_out_file.empty()) {
      // WRITE ALL THE HMM PARAMETERS:
      write_params_file(params_out_file,
			reps_fg_alpha, reps_fg_beta,
			reps_bg_alpha, reps_bg_beta,
			start_trans, trans, end_trans);
    }

    /***********************************/

    vector<bool> classes;
    vector<double> scores;
    hmm.PosteriorDecoding_rep(meth, reset_points, start_trans, trans, end_trans,
			      reps_fg_alpha, reps_fg_beta,
			      reps_bg_alpha, reps_bg_beta, classes, scores);

    vector<double> domain_scores;
    get_domain_scores_rep(classes, meth, reset_points, domain_scores);

    vector<double> random_scores;
    shuffle_cpgs_rep(hmm, meth, reset_points, start_trans, trans, end_trans,
		     reps_fg_alpha, reps_fg_beta, reps_bg_alpha, reps_bg_beta,
                     random_scores);

    vector<double> p_values;
    assign_p_values(random_scores, domain_scores, p_values);

    if (fdr_cutoff == numeric_limits<double>::max())
      fdr_cutoff = get_fdr_cutoff(p_values, 0.01);
    if (!params_out_file.empty()) {
      std::ofstream out(params_out_file.c_str(), std::ios::app);
      out << "FDR_CUTOFF\t"
          << std::setprecision(30) << fdr_cutoff << endl;
      out.close();
    }

    vector<GenomicRegion> domains;
    build_domains(VERBOSE, cpgs[0], scores, reset_points, classes, domains);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    size_t good_hmr_count = 0;
    for (size_t i = 0; i < domains.size(); ++i)
      if (p_values[i] < fdr_cutoff) {
	domains[i].set_name("HYPO" + smithlab::toa(good_hmr_count++));
	out << domains[i] << '\n';
      }

  }
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
