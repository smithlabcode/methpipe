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
#include <stdexcept>

#include <unistd.h>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "zlib_wrapper.hpp"

#include "TwoStateHMM.hpp"
#include "MethpipeSite.hpp"

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
using std::runtime_error;

static GenomicRegion
as_gen_rgn(const MSite &s) {
  return GenomicRegion(s.chrom, s.pos, s.pos + 1);
}

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
  static const bool CLASS_ID = true;

  const size_t n_replicates = meth.size();
  size_t reset_idx = 1;
  bool in_domain = false;
  double score = 0.0;
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
      for (size_t r = 0; r < n_replicates ; ++r) {
        if (meth[r][i].first + meth[r][i].second >= 1) {
          score += 1.0 - (meth[r][i].first/(meth[r][i].first +
                                            meth[r][i].second));
        }
      }
    }
    else if (in_domain) {
      in_domain = false;
      scores.push_back(score);
      score = 0;
    }
  }
}


static void
build_domains(const bool VERBOSE,
              const vector<MSite> &cpgs,
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
        domains.push_back(as_gen_rgn(cpgs[i]));
        domains.back().set_name("HYPO" + toa(n_domains++));
      }
      ++n_cpgs;
      score += post_scores[i];
    }
    else if (in_domain) {
      in_domain = false;
      domains.back().set_end(prev_end);
      domains.back().set_score(n_cpgs);
      n_cpgs = 0;
      score = 0;
    }
    prev_end = cpgs[i].pos + 1;
  }
}


template <class T, class U>
static void
separate_regions(const bool VERBOSE,
                 const size_t desert_size,
                 vector<MSite> &cpgs,
                 vector<vector<T>> &meth,
                 vector<vector<U>> &reads,
                 vector<size_t> &reset_points) {
  if (VERBOSE)
    cerr << "[SEPARATING BY CPG DESERT]" << endl;

  // eliminate the zero-read cpg sites if no coverage in any replicates
  const size_t n_replicates = meth.size();

  size_t j = 0;
  bool all_empty;
  size_t rep_idx;
  for (size_t i = 0; i < cpgs.size(); ++i) {
    all_empty = true;
    rep_idx = 0;
    while (all_empty && rep_idx < n_replicates)
      if (reads[rep_idx][i] == 0)
        ++rep_idx;
      else
        all_empty = false;

    if (!all_empty) {
      // swap(cpgs[r][j], cpgs[r][i]);
      cpgs[j] = cpgs[i];
      for (size_t r = 0; r < n_replicates; ++r) {
        // swap(meth[r][j], meth[r][i]);
        // swap(reads[r][j], reads[r][i]);
        meth[r][j] = meth[r][i];
        reads[r][j] = reads[r][i];
      }
      ++j;
    }
  }

  cpgs.erase(begin(cpgs) + j, end(cpgs));
  for (size_t r = 0; r < n_replicates; ++r) {
    meth[r].erase(begin(meth[r]) + j, end(meth[r]));
    reads[r].erase(begin(reads[r]) + j, end(reads[r]));
  }

  // segregate cpgs
  size_t prev_cpg = 0;
  for (size_t i = 0; i < cpgs.size(); ++i) {
    const size_t dist = (i > 0 && cpgs[i].chrom == cpgs[i - 1].chrom) ?
      cpgs[i].pos - prev_cpg : numeric_limits<size_t>::max();
    if (dist > desert_size)
      reset_points.push_back(i);
    prev_cpg = cpgs[i].pos;
  }
  reset_points.push_back(cpgs.size());

  if (VERBOSE)
    cerr << "CPGS RETAINED: " << cpgs.size() << endl
         << "DESERTS REMOVED: " << reset_points.size() - 2
         << endl << endl;
}

static void
shuffle_cpgs_rep(const size_t seed,
                 const TwoStateHMM &hmm,
                 vector<vector<pair<double, double> > > meth,
                 vector<size_t> reset_points,
                 const vector<double> &start_trans,
                 const vector<vector<double> > &trans,
                 const vector<double> &end_trans,
                 const vector<double> &fg_alpha,
                 const vector<double> &fg_beta,
                 const vector<double> &bg_alpha,
                 const vector<double> &bg_beta,
                 vector<double> &domain_scores) {

  srand(seed);
  const size_t n_replicates = meth.size();

  for (size_t r =0 ; r < n_replicates; ++r)
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
  std::ifstream in(params_file);

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
  if (!outfile.empty()) of.open(outfile);
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

static void
load_cpgs(const string &cpgs_file,
          vector<MSite> &cpgs,
          vector<pair<double, double> > &meth,
          vector<size_t> &reads) {

  igzfstream in(cpgs_file);
  if (!in)
    throw runtime_error("failed opening file: " + cpgs_file);

  MSite the_site;
  while (in >> the_site) {
    cpgs.push_back(the_site);
    reads.push_back(the_site.n_reads);
    const double m = reads.back()*the_site.meth;
    meth.push_back(make_pair(m, reads.back() - m));
  }
}

int
main(int argc, const char **argv) {

  try {

    const char* sep = ",";
    string outfile;
    string hypo_post_outfile;
    string meth_post_outfile;
    size_t seed = 408;

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
    opt_parse.add_opt("post-hypo", '\0', "output file for single-CpG posteiror "
                      "hypomethylation probability (default: NULL)",
                      false, hypo_post_outfile);
    opt_parse.add_opt("post-meth", '\0', "output file for single-CpG posteiror "
                      "methylation probability (default: NULL)",
                      false, meth_post_outfile);
    opt_parse.add_opt("params-in", 'P', "HMM parameter files for "
                      "individual methylomes (separated with comma)",
                      false, params_in_files);
    opt_parse.add_opt("params-out", 'p', "write HMM parameters to this file",
                      false, params_out_file);
    opt_parse.add_opt("seed", 's', "specify random seed", false, seed);
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
    const vector<string> cpgs_files(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    vector<string> params_in_file;
    if (!params_in_files.empty()) {
      params_in_file = smithlab::split(params_in_files, sep, false);
      assert(cpgs_files.size() == params_in_file.size());
    }

    const size_t n_replicates = cpgs_files.size();

    vector<MSite> cpgs;
    vector<vector<pair<double, double>>> meth(n_replicates);
    vector<vector<size_t>> reads(n_replicates);

    for (size_t i = 0; i < n_replicates; ++i) {
      if (VERBOSE)
        cerr << "[READING CPGS AND METH PROPS] from " << cpgs_files[i] << endl;
      vector<MSite> cpgs_rep;
      load_cpgs(cpgs_files[i], cpgs_rep, meth[i], reads[i]);
      if (VERBOSE)
        cerr << "TOTAL CPGS: " << cpgs_rep.size() << endl
             << "MEAN COVERAGE: "
             << accumulate(begin(reads[i]), end(reads[i]), 0.0)/reads[i].size()
             << endl;
      if (i > 0 && cpgs_rep.size() != cpgs.size())
        throw runtime_error("inconsistent number of sites");
      else swap(cpgs, cpgs_rep);
    }

    // separate the regions by chrom and by desert, and eliminate
    // those isolated CpGs
    vector<size_t> reset_points;
    separate_regions(VERBOSE, desert_size, cpgs, meth, reads, reset_points);

    /****************** Read in params *****************/
    vector<double> start_trans(2, 0.5), end_trans(2, 1e-10);
    vector<vector<double> > trans(2, vector<double>(2, 0.25));
    trans[0][0] = trans[1][1] = 0.75;
    const TwoStateHMM hmm(min_prob, tolerance, max_iterations, VERBOSE, DEBUG);
    vector<double> reps_fg_alpha(n_replicates, 0);
    vector<double> reps_fg_beta(n_replicates, 0);
    vector<double> reps_bg_alpha(n_replicates, 0);
    vector<double> reps_bg_beta(n_replicates, 0);
    double fdr_cutoff = std::numeric_limits<double>::max();

    if (!params_in_file.empty()) {
      // READ THE PARAMETERS FILES
      double fdr_cutoff_rep;
      for (size_t i = 0; i< n_replicates; ++i)
        read_params_file(VERBOSE, params_in_file[i], reps_fg_alpha[i],
                         reps_fg_beta[i], reps_bg_alpha[i], reps_bg_beta[i],
                         start_trans, trans, end_trans, fdr_cutoff_rep);
      max_iterations = 0;
    }
    else {
      const size_t n_sites = reads.front().size();
      for (size_t i = 0; i < n_replicates; ++i) {
        // actually, there are many 0s in reads[r], but the parameter
        // start points don't need to be accurate
        const double mean_reads =
          accumulate(begin(reads[i]), end(reads[i]), 0.0)/n_sites;
        reps_fg_alpha[i] = 0.33*mean_reads;
        reps_fg_beta[i] = 0.67*mean_reads;
        reps_bg_alpha[i] = 0.67*mean_reads;
        reps_bg_beta[i] = 0.33*mean_reads;
      }
    }

    if (max_iterations > 0)
      hmm.BaumWelchTraining_rep(meth, reset_points,
                                start_trans, trans, end_trans,
                                reps_fg_alpha, reps_fg_beta,
                                reps_bg_alpha, reps_bg_beta);

    vector<bool> classes;
    vector<double> scores;
    hmm.PosteriorDecoding_rep(meth, reset_points, start_trans, trans, end_trans,
                              reps_fg_alpha, reps_fg_beta,
                              reps_bg_alpha, reps_bg_beta, classes, scores);

    vector<double> domain_scores;
    get_domain_scores_rep(classes, meth, reset_points, domain_scores);

    vector<double> random_scores;
    shuffle_cpgs_rep(seed, hmm, meth, reset_points,
                     start_trans, trans, end_trans,
                     reps_fg_alpha, reps_fg_beta, reps_bg_alpha, reps_bg_beta,
                     random_scores);

    vector<double> p_values;
    assign_p_values(random_scores, domain_scores, p_values);

    if (fdr_cutoff == numeric_limits<double>::max())
      fdr_cutoff = get_fdr_cutoff(p_values, 0.01);

    vector<GenomicRegion> domains;
    build_domains(VERBOSE, cpgs, scores, reset_points, classes, domains);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    size_t good_hmr_count = 0;
    for (size_t i = 0; i < domains.size(); ++i)
      if (p_values[i] < fdr_cutoff) {
        domains[i].set_name("HYPO" + smithlab::toa(good_hmr_count++));
        out << domains[i] << '\n';
      }
    /***********************************
     * STEP 6: (OPTIONAL) WRITE POSTERIOR
     */

    // write all the hmm parameters if requested
    if (!params_out_file.empty()) {
      write_params_file(params_out_file,
                        reps_fg_alpha, reps_fg_beta,
                        reps_bg_alpha, reps_bg_beta,
                        start_trans, trans, end_trans);
      std::ofstream out(params_out_file, std::ios::app);
      out << "FDR_CUTOFF\t"
          << std::setprecision(30) << fdr_cutoff << endl;
    }

    if (!hypo_post_outfile.empty() || !meth_post_outfile.empty()) {
      bool fg_class = true;
      vector<double> fg_posterior;
      hmm.PosteriorScores_rep(meth, reset_points, start_trans, trans,
                              end_trans, reps_fg_alpha, reps_fg_beta,
                              reps_bg_alpha, reps_bg_beta,
                              fg_class, fg_posterior);

      if (!hypo_post_outfile.empty()) {
        if (VERBOSE)
          cerr << "[WRITING " << hypo_post_outfile
               << " (4th field: CpG:<M_reads>:<U_reads>)]" << endl;
        std::ofstream out(hypo_post_outfile);
        for (size_t i = 0; i < cpgs.size(); ++i) {
          GenomicRegion cpg(as_gen_rgn(cpgs[i]));
          size_t M_reads = 0;
          size_t U_reads = 0;
          for (size_t j = 0; j < n_replicates; ++j){
            M_reads += static_cast<size_t>(meth[j][i].first);
            U_reads += static_cast<size_t>(meth[j][i].second);
          }
          cpg.set_name("CpG:" + toa(M_reads) + ":" + toa(U_reads));
          cpg.set_score(scores[i]);
          out << cpg << '\n';
        }
      }

      if (!meth_post_outfile.empty()) {
        std::ofstream out(meth_post_outfile);
        if (VERBOSE)
          cerr << "[WRITING " << meth_post_outfile
               << " (4th field: CpG:<M_reads>:<U_reads>)]" << endl;
        for (size_t i = 0; i < cpgs.size(); ++i) {
          GenomicRegion cpg(as_gen_rgn(cpgs[i]));
          size_t M_reads = 0;
          size_t U_reads = 0;
          for (size_t j = 0; j < n_replicates; ++j){
            M_reads += static_cast<size_t>(meth[j][i].first);
            U_reads += static_cast<size_t>(meth[j][i].second);
          }
          cpg.set_name("CpG:" +toa(M_reads) + ":" + toa(U_reads));
          cpg.set_score(1.0 - scores[i]);
          out << cpg << '\n';
        }
      }
    }
  }
  catch (const runtime_error &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
