/* Copyright (C) 2019 University of Southern California
 *                    Andrew D Smith
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
#include <string>
#include <stdexcept>

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
using std::to_string;
using std::begin;
using std::end;

static GenomicRegion
as_gen_rgn(const MSite &s) {
  return GenomicRegion(s.chrom, s.pos, s.pos + 1);
}

static double
get_stepup_cutoff(vector<double> scores, const double cutoff) {
  if (cutoff <= 0) return numeric_limits<double>::max();
  else if (cutoff > 1) return numeric_limits<double>::min();

  const size_t n = scores.size();
  std::sort(begin(scores), end(scores));
  size_t i = 1;
  while (i < n && scores[i-1] < (cutoff*i)/n) ++i;
  return scores[i - 1];
}

template <class T> T
pair_sum(const std::pair<T, T> &t) {return t.first + t.second;}

static void
get_domain_scores_rep(const vector<bool> &state_ids,
                      const vector<vector<pair<double, double> > > &meth,
                      const vector<size_t> &reset_points,
                      vector<double> &scores) {

  const size_t n_reps = meth.size();
  size_t reset_idx = 1;
  bool in_domain = false;
  double score = 0.0;
  for (size_t i = 0; i < state_ids.size(); ++i) {
    if (reset_points[reset_idx] == i) {
      if (in_domain) {
        in_domain = false;
        scores.push_back(score);
        score = 0;
      }
      ++reset_idx;
    }
    if (state_ids[i]) {
      in_domain = true;
      for (size_t r = 0; r < n_reps; ++r)
        if (pair_sum(meth[r][i]) >= 1)
          score += 1.0 - meth[r][i].first/pair_sum(meth[r][i]);
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
              const vector<bool> &state_ids,
              vector<GenomicRegion> &domains) {

  size_t n_cpgs = 0, n_domains = 0, reset_idx = 1, prev_end = 0;
  bool in_domain = false;
  double score = 0;
  for (size_t i = 0; i < state_ids.size(); ++i) {
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
    if (state_ids[i]) {
      if (!in_domain) {
        in_domain = true;
        domains.push_back(as_gen_rgn(cpgs[i]));
        domains.back().set_name("HYPO" + to_string(n_domains++));
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
                 vector<vector<T> > &meth,
                 vector<vector<U> > &reads,
                 vector<size_t> &reset_points) {
  if (VERBOSE)
    cerr << "[separating by cpg desert]" << endl;

  // eliminate the zero-read cpg sites if no coverage in any replicates
  const size_t n_reps = meth.size();

  size_t j = 0;
  for (size_t i = 0; i < cpgs.size(); ++i) {
    bool has_data = true;
    for (size_t rep_idx = 0; rep_idx < n_reps; ++rep_idx)
      has_data = (has_data && reads[rep_idx][i] > 0);

    if (has_data) {
      cpgs[j] = cpgs[i];
      for (size_t r = 0; r < n_reps; ++r) {
        meth[r][j] = meth[r][i];
        reads[r][j] = reads[r][i];
      }
      ++j;
    }
  }

  cpgs.erase(begin(cpgs) + j, end(cpgs));
  for (size_t r = 0; r < n_reps; ++r) {
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
    cerr << "[cpgs retained: " << cpgs.size() << "]" << endl
         << "[deserts removed: " << reset_points.size() - 2 << "]" << endl;
}

static void
shuffle_cpgs_rep(const size_t seed, const TwoStateHMM &hmm,
                 vector<vector<pair<double, double> > > meth,
                 vector<size_t> reset_points,
                 const double f_to_b_trans, const double b_to_f_trans,
                 const vector<double> &fg_alpha, const vector<double> &fg_beta,
                 const vector<double> &bg_alpha, const vector<double> &bg_beta,
                 vector<double> &domain_scores) {

  srand(seed);
  const size_t n_reps = meth.size();

  for (size_t r = 0 ; r < n_reps; ++r)
    random_shuffle(begin(meth[r]), end(meth[r]));

  vector<bool> state_ids;
  vector<double> scores;
  hmm.PosteriorDecoding(meth, reset_points, f_to_b_trans, b_to_f_trans,
                        fg_alpha, fg_beta, bg_alpha, bg_beta,
                        state_ids, scores);

  get_domain_scores_rep(state_ids, meth, reset_points, domain_scores);
  sort(begin(domain_scores), end(domain_scores));
}


static void
assign_p_values(const vector<double> &random_scores,
                const vector<double> &observed_scores,
                vector<double> &p_values) {
  const double n_randoms = max(random_scores.size(), 1ul);
  for (size_t i = 0; i < observed_scores.size(); ++i)
    p_values.push_back((end(random_scores) -
                        upper_bound(begin(random_scores),
                                    end(random_scores),
                                    observed_scores[i]))/n_randoms);
}


static void
read_params_file(const bool VERBOSE,
                 const string &params_file,
                 double &fg_alpha, double &fg_beta,
                 double &bg_alpha, double &bg_beta,
                 double &f_to_b_trans, double &b_to_f_trans,
                 double &fdr_cutoff) {

  string jnk;
  std::ifstream in(params_file);
  if (!in)
    throw runtime_error("failed to parse params file: " + params_file);

  in >> jnk >> fg_alpha
     >> jnk >> fg_beta
     >> jnk >> bg_alpha
     >> jnk >> bg_beta
     >> jnk >> f_to_b_trans
     >> jnk >> b_to_f_trans
     >> jnk >> fdr_cutoff;

  if (VERBOSE)
    cerr << "read in params from " << params_file << endl
         << "FG_ALPHA\t" << fg_alpha << endl
         << "FG_BETA\t" << fg_beta << endl
         << "BG_ALPHA\t" << bg_alpha << endl
         << "BG_BETA\t" << bg_beta << endl
         << "F_B\t" << f_to_b_trans << endl
         << "B_F\t" << b_to_f_trans << endl
         << "FDR_CUTOFF\t" << fdr_cutoff << endl;
}


static void
write_params_file(const string &outfile,
                  const vector<double> &fg_alpha,
                  const vector<double> &fg_beta,
                  const vector<double> &bg_alpha,
                  const vector<double> &bg_beta,
                  const double f_to_b_trans,
                  const double b_to_f_trans,
                  const double fdr_cutoff) {

  std::ofstream of;
  if (!outfile.empty()) of.open(outfile);
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out.precision(30);
  for (size_t r =0; r < fg_alpha.size(); ++r)
    out << "FG_ALPHA_" << r+1 << '\t' << fg_alpha[r] << '\t'
        << "FG_BETA_" << r+1 << '\t' << fg_beta[r] << '\t'
        << "BG_ALPHA_" << r+1 << '\t' << bg_alpha[r] << '\t'
        << "BG_BETA_" << r+1 << '\t' << bg_beta[r] << endl;

  out << "F_B\t" << f_to_b_trans << endl
      << "B_F\t" << b_to_f_trans << endl
      << "FDR_CUTOFF\t" << fdr_cutoff << endl
    ;
}

static void
load_cpgs(const string &cpgs_file, vector<MSite> &cpgs,
          vector<pair<double, double> > &meth,
          vector<uint32_t> &reads) {

  igzfstream in(cpgs_file);
  if (!in)
    throw runtime_error("failed opening file: " + cpgs_file);

  MSite the_site;
  while (in >> the_site) {
    cpgs.push_back(the_site);
    reads.push_back(the_site.n_reads);
    meth.push_back(make_pair(the_site.n_meth(), the_site.n_unmeth()));
  }
}

static void
check_consistent_sites(const string &expected_filename,
                       const vector<MSite> &expected,
                       const string &observed_filename,
                       const vector<MSite> &observed) {
  if (expected.size() != observed.size()) {
    std::ostringstream err_msg;
    err_msg << "inconsistent number of sites" << endl
            << "file=" << expected_filename << ","
            << "sites=" << expected.size() << endl
            << "file=" << observed_filename << ","
            << "sites=" << observed.size() << endl;
    throw runtime_error(err_msg.str());
  }
}

template <class InputIterator> double
get_mean(InputIterator first, InputIterator last) {
  return accumulate(first, last, 0.0)/std::distance(first, last);
}

static vector<string>
split_comma(const string &orig) {
  string tmp(orig);
  replace(begin(tmp), end(tmp), ',', ' ');
  std::istringstream iss(tmp);
  vector<string> parts;
  while (iss >> tmp)
    parts.push_back(tmp);
  return parts;
}

int
main_hmr_rep(int argc, const char **argv) {

  try {

    string outfile;
    string hypo_post_outfile;
    string meth_post_outfile;

    size_t desert_size = 1000;
    size_t max_iterations = 10;
    size_t seed = 408;

    // run mode flags
    bool VERBOSE = false;

    const double tolerance = 1e-10; // corrections for small values

    string params_in_files;
    string params_out_file;

    const string description =
      "Identify HMRs in a set of replicate methylomes. Methylation must be \
      provided in the methcounts format (chrom, position, strand, context, \
      methylation, reads). See the methcounts documentation for details    \
      for details. This program assumes only data at CpG sites and that    \
      strands are collapsed so only the positive site appears in the file.";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<methcount-file-1> <methcount-file-2> ...");
    opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("desert", 'd', "max dist btwn cpgs with reads in HMR",
                      false, desert_size);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
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
    opt_parse.set_show_defaults();
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
      params_in_file = split_comma(params_in_files);
      assert(cpgs_files.size() == params_in_file.size());
    }

    const size_t n_reps = cpgs_files.size();

    vector<MSite> cpgs;
    vector<vector<pair<double, double> > > meth(n_reps);
    vector<vector<uint32_t> > reads(n_reps);

    if (VERBOSE)
      cerr << "[reading methylation levels]" << endl;
    for (size_t i = 0; i < n_reps; ++i) {
      if (VERBOSE)
        cerr << "[filename=" << cpgs_files[i] << "]" << endl;
      vector<MSite> curr_rep;
      load_cpgs(cpgs_files[i], curr_rep, meth[i], reads[i]);
      if (VERBOSE)
        cerr << "[total_cpgs=" << curr_rep.size() << "]" << endl
             << "[mean_coverage="
             << get_mean(begin(reads[i]), end(reads[i])) << "]" << endl;
      if (i > 0)
        check_consistent_sites(cpgs_files[0], cpgs, cpgs_files[i], curr_rep);
      else swap(cpgs, curr_rep);
    }

    // separate the regions by chrom and by desert, and eliminate
    // those isolated CpGs
    vector<size_t> reset_points;
    separate_regions(VERBOSE, desert_size, cpgs, meth, reads, reset_points);

    /****************** initalize params *****************/
    const TwoStateHMM hmm(tolerance, max_iterations, VERBOSE);
    vector<double> fg_alpha(n_reps), fg_beta(n_reps);
    vector<double> bg_alpha(n_reps), bg_beta(n_reps);
    double fdr_cutoff = std::numeric_limits<double>::max();

    double f_to_b_trans = 0.25;
    double b_to_f_trans = 0.25;

    if (!params_in_file.empty()) { // read parameters files
      double fdr_cutoff_rep; // ignore this cutoff
      for (size_t i = 0; i < n_reps; ++i)
        read_params_file(VERBOSE, params_in_file[i], fg_alpha[i],
                         fg_beta[i], bg_alpha[i], bg_beta[i],
                         f_to_b_trans, b_to_f_trans, fdr_cutoff_rep);
      max_iterations = 0;
    }
    else {
      for (size_t i = 0; i < n_reps; ++i) {
        // JQU: there are many 0s in reads[r], but the parameter start
        // points don't need to be perfect
        const double mean_reads = get_mean(begin(reads[i]), end(reads[i]));
        fg_alpha[i] = 0.33*mean_reads;
        fg_beta[i] = 0.67*mean_reads;
        bg_alpha[i] = 0.67*mean_reads;
        bg_beta[i] = 0.33*mean_reads;
      }
    }

    if (max_iterations > 0)
      hmm.BaumWelchTraining(meth, reset_points, f_to_b_trans, b_to_f_trans,
                            fg_alpha, fg_beta, bg_alpha, bg_beta);

    vector<bool> state_ids;
    vector<double> posteriors;
    hmm.PosteriorDecoding(meth, reset_points, f_to_b_trans, b_to_f_trans,
                          fg_alpha, fg_beta, bg_alpha, bg_beta,
                          state_ids, posteriors);

    vector<double> domain_scores;
    get_domain_scores_rep(state_ids, meth, reset_points, domain_scores);

    vector<double> random_scores;
    shuffle_cpgs_rep(seed, hmm, meth, reset_points, f_to_b_trans, b_to_f_trans,
                     fg_alpha, fg_beta, bg_alpha, bg_beta, random_scores);

    vector<double> p_values;
    assign_p_values(random_scores, domain_scores, p_values);

    if (fdr_cutoff == numeric_limits<double>::max())
      fdr_cutoff = get_stepup_cutoff(p_values, 0.01);

    vector<GenomicRegion> domains;
    build_domains(VERBOSE, cpgs, posteriors, reset_points, state_ids, domains);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    size_t good_hmr_count = 0;
    for (size_t i = 0; i < domains.size(); ++i)
      if (p_values[i] < fdr_cutoff) {
        domains[i].set_name("HYPO" + to_string(good_hmr_count++));
        out << domains[i] << '\n';
      }

    // write all the hmm parameters if requested
    if (!params_out_file.empty())
      write_params_file(params_out_file, fg_alpha, fg_beta,
                        bg_alpha, bg_beta, f_to_b_trans, b_to_f_trans,
                        fdr_cutoff);

    if (!hypo_post_outfile.empty()) {
      if (VERBOSE)
        cerr << "[writing=" << hypo_post_outfile << "]" << endl;
      std::ofstream out(hypo_post_outfile);
      for (size_t i = 0; i < cpgs.size(); ++i) {
        size_t m_reads = 0, u_reads = 0;
        for (size_t j = 0; j < n_reps; ++j){
          m_reads += meth[j][i].first;
          u_reads += meth[j][i].second;
        }
        GenomicRegion cpg(as_gen_rgn(cpgs[i]));
        cpg.set_name("CpG:" + to_string(m_reads) + ":" + to_string(u_reads));
        cpg.set_score(posteriors[i]);
        out << cpg << '\n';
      }
    }

    if (!meth_post_outfile.empty()) {
      std::ofstream out(meth_post_outfile);
      if (VERBOSE)
        cerr << "[writing=" << meth_post_outfile << "]" << endl;
      for (size_t i = 0; i < cpgs.size(); ++i) {
        size_t m_reads = 0, u_reads = 0;
        for (size_t j = 0; j < n_reps; ++j) {
          m_reads += meth[j][i].first;
          u_reads += meth[j][i].second;
        }
        GenomicRegion cpg(as_gen_rgn(cpgs[i]));
        cpg.set_name("CpG:" + to_string(m_reads) + ":" + to_string(u_reads));
        cpg.set_score(1.0 - posteriors[i]);
        out << cpg << '\n';
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
