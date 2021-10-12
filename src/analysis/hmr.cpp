/* Copyright (C) 2009-2019 University of Southern California
 *                         Andrew D Smith
 *
 * Author: Andrew D. Smith, Song Qiang
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

static GenomicRegion
as_gen_rgn(const MSite &s) {
  return GenomicRegion(s.chrom, s.pos, s.pos + 1);
}

static string
format_cpg_meth_tag(const pair<double, double> &m) {
  return "CpG:" + to_string(static_cast<size_t>(m.first)) +
    ":" + to_string(static_cast<size_t>(m.second));
}

double
get_stepup_cutoff(vector<double> scores, const double cutoff) {
  if (cutoff <= 0) return numeric_limits<double>::max();
  else if (cutoff > 1) return numeric_limits<double>::min();

  const size_t n = scores.size();
  std::sort(begin(scores), end(scores));
  size_t i = 1;
  while (i < n && scores[i-1] < (cutoff*i)/n) ++i;
  return scores[i - 1];
}

static void
get_domain_scores(const vector<bool> &state_ids,
                  const vector<pair<double, double> > &meth,
                  const vector<size_t> &reset_points,
                  vector<double> &scores) {

  size_t n_cpgs = 0, reset_idx = 1;
  bool in_domain = false;
  double score = 0;
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
      score += 1.0 - (meth[i].first/(meth[i].first + meth[i].second));
      ++n_cpgs;
    }
    else if (in_domain) {
      in_domain = false;
      scores.push_back(score);
      score = 0;
    }
  }

  if (in_domain) {
    scores.push_back(score);
  }
}


static void
build_domains(const bool VERBOSE,
              const vector<MSite> &cpgs,
              const vector<double> &post_scores,
              const vector<size_t> &reset_points,
              const vector<bool> &state_ids,
              vector<GenomicRegion> &domains) {

  size_t n_cpgs = 0, reset_idx = 1, prev_end = 0;
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
    if (state_ids[i]) { // currently in an hmr
      if (!in_domain) {
        in_domain = true;
        domains.push_back(as_gen_rgn(cpgs[i]));
        domains.back().set_name("HYPO" + to_string(domains.size()));
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
  if (in_domain) {
    domains.back().set_end(prev_end);
    domains.back().set_score(n_cpgs);
  }
}


template <class T, class U>
static void
separate_regions(const bool VERBOSE,
                 const size_t desert_size,
                 vector<MSite> &cpgs,
                 vector<T> &meth, vector<U> &reads,
                 vector<size_t> &reset_points) {
  if (VERBOSE)
    cerr << "[separating by cpg desert]" << endl;
  // eliminate the zero-read cpgs
  size_t j = 0;
  for (size_t i = 0; i < cpgs.size(); ++i)
    if (reads[i] > 0) {
      cpgs[j] = cpgs[i];
      meth[j] = meth[i];
      reads[j] = reads[i];
      ++j;
    }
  cpgs.erase(begin(cpgs) + j, end(cpgs));
  meth.erase(begin(meth) + j, end(meth));
  reads.erase(begin(reads) + j, end(reads));

  double total_bases = 0;
  double bases_in_deserts = 0;
  // segregate cpgs
  size_t prev_pos = 0;
  for (size_t i = 0; i < cpgs.size(); ++i) {
    const size_t dist = (i > 0 && cpgs[i].chrom == cpgs[i - 1].chrom) ?
      cpgs[i].pos - prev_pos : numeric_limits<size_t>::max();
    if (dist > desert_size) {
      reset_points.push_back(i);
      if (dist < numeric_limits<size_t>::max())
        bases_in_deserts += dist;
    }
    if (dist<numeric_limits<size_t>::max())
      total_bases += dist;

    prev_pos = cpgs[i].pos;
  }
  reset_points.push_back(cpgs.size());

  if (VERBOSE)
    cerr << "[cpgs retained: " << cpgs.size() << "]" << endl
         << "[deserts removed: " << reset_points.size() - 2 << "]" << endl
         << "[genome fraction covered: "
         << 1.0 - (bases_in_deserts/total_bases) << "]" << endl;
}


/* function to "fold" the methylation profile so that the middle
 * methylation becomes lower methylation, and both the low and high
 * methylation become high. this method actually seems to work.
 */
static void
make_partial_meth(const vector<uint32_t> &reads,
                  vector<pair<double, double> > &meth) {
  for (size_t i = 0; i < reads.size(); ++i) {
    double m = meth[i].first/reads[i];
    m = (m <= 0.5) ? (1.0 - 2*m) : (1.0 - 2*(1.0 - m));
    meth[i].first = reads[i]*m;
    meth[i].second = (reads[i] - meth[i].first);
  }
}

static void
shuffle_cpgs(const size_t seed,
             const TwoStateHMM &hmm,
             vector<pair<double, double> > meth,
             vector<size_t> reset_points,
             const double p_fb, const double p_bf,
             const double fg_alpha, const double fg_beta,
             const double bg_alpha, const double bg_beta,
             vector<double> &domain_scores) {
  srand(seed);
  random_shuffle(begin(meth), end(meth));
  vector<bool> state_ids;
  vector<double> scores;
  hmm.PosteriorDecoding(meth, reset_points, p_fb, p_bf,
                        fg_alpha, fg_beta, bg_alpha,
                        bg_beta, state_ids, scores);
  get_domain_scores(state_ids, meth, reset_points, domain_scores);
  sort(begin(domain_scores), end(domain_scores));
}

static void
assign_p_values(const vector<double> &random_scores,
                const vector<double> &observed_scores,
                vector<double> &p_values) {
  const double n_randoms = random_scores.empty() ? 1 : random_scores.size();
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
                 double &p_fb, double &p_bf,
                 double &domain_score_cutoff) {
  string jnk;
  std::ifstream in(params_file);
  if (!in)
    throw runtime_error("failed to parse params file: " + params_file);
  in >> jnk >> fg_alpha
     >> jnk >> fg_beta
     >> jnk >> bg_alpha
     >> jnk >> bg_beta
     >> jnk >> p_fb
     >> jnk >> p_bf
     >> jnk >> domain_score_cutoff;
  if (VERBOSE)
    cerr << "FG_ALPHA\t" << fg_alpha << endl
         << "FG_BETA\t" << fg_beta << endl
         << "BG_ALPHA\t" << bg_alpha << endl
         << "BG_BETA\t" << bg_beta << endl
         << "F_B\t" << p_fb << endl
         << "B_F\t" << p_bf << endl
         << "DOMAIN_SCORE_CUTOFF\t" << domain_score_cutoff << endl;
}

static void
write_params_file(const string &outfile,
                  const double fg_alpha,
                  const double fg_beta,
                  const double bg_alpha,
                  const double bg_beta,
                  const double p_fb,
                  const double p_bf,
                  const double domain_score_cutoff) {

  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out.precision(30);
  out << "FG_ALPHA\t" << fg_alpha << endl
      << "FG_BETA\t" << fg_beta << endl
      << "BG_ALPHA\t" << bg_alpha << endl
      << "BG_BETA\t" << bg_beta << endl
      << "F_B\t" << p_fb << endl
      << "B_F\t" << p_bf << endl
      << "DOMAIN_SCORE_CUTOFF\t" << domain_score_cutoff << endl;
}

static void
load_cpgs(const string &cpgs_file, vector<MSite> &cpgs,
          vector<pair<double, double> > &meth,
          vector<uint32_t> &reads) {

  igzfstream in(cpgs_file);
  if (!in)
    throw runtime_error("failed opening file: " + cpgs_file);

  MSite prev_site, the_site;
  while (in >> the_site) {
    if (!the_site.is_cpg() || distance(prev_site, the_site) < 2)
      throw runtime_error("error: input is not symmetric-CpGs: " + cpgs_file);
    cpgs.push_back(the_site);
    reads.push_back(the_site.n_reads);
    meth.push_back(make_pair(the_site.n_meth(), the_site.n_unmeth()));
    prev_site = the_site;
  }
}

template <class InputIterator> double
get_mean(InputIterator first, InputIterator last) {
  return accumulate(first, last, 0.0)/std::distance(first, last);
}

int
main(int argc, const char **argv) {

  try {

    string outfile;
    string hypo_post_outfile;
    string meth_post_outfile;

    size_t desert_size = 1000;
    size_t max_iterations = 10;
    size_t seed = 408;

    // run mode flags
    bool VERBOSE = false;
    bool PARTIAL_METH = false;

    // corrections for small values
    const double tolerance = 1e-10;

    string params_in_file;
    string params_out_file;

    const string description =
      "Identify HMRs in methylomes. Methylation must be provided in the \
      methcounts format (chrom, position, strand, context,              \
      methylation, reads). See the methcounts documentation for         \
      details. This program assumes only data at CpG sites and that     \
      strands are collapsed so only the positive site appears in the    \
      file.";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           description, "<methylation-file>");
    opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("desert", 'd', "max dist btwn covered cpgs in HMR",
                      false, desert_size);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("partial", '\0', "identify PMRs instead of HMRs",
                      false, PARTIAL_METH);
    opt_parse.add_opt("post-hypo", '\0', "output file for single-CpG posterior "
                      "hypomethylation probability (default: none)",
                      false, hypo_post_outfile);
    opt_parse.add_opt("post-meth", '\0', "output file for single-CpG posteiror "
                      "methylation probability (default: none)",
                      false, meth_post_outfile);
    opt_parse.add_opt("params-in", 'P', "HMM parameter file "
                      "(override training)", false, params_in_file);
    opt_parse.add_opt("params-out", 'p', "write HMM parameters to this "
                      "file (default: none)", false, params_out_file);
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
    const string cpgs_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    // separate the regions by chrom and by desert
    vector<MSite> cpgs;
    vector<pair<double, double> > meth;
    vector<uint32_t> reads;
    if (VERBOSE)
      cerr << "[reading methylation levels]" << endl;
    load_cpgs(cpgs_file, cpgs, meth, reads);

    if (!std::is_sorted(begin(cpgs), end(cpgs)))
      throw runtime_error("error: input is not properly sorted: " + cpgs_file);

    if (PARTIAL_METH)
      make_partial_meth(reads, meth);

    if (VERBOSE)
      cerr << "[total_cpgs=" << cpgs.size() << "]" << endl
           << "[mean_coverage="
           << get_mean(begin(reads), end(reads)) << "]" << endl;

    // separate the regions by chrom and by desert, and eliminate
    // those isolated CpGs
    vector<size_t> reset_points;
    separate_regions(VERBOSE, desert_size, cpgs, meth, reads, reset_points);

    const TwoStateHMM hmm(tolerance, max_iterations, VERBOSE);

    double p_fb = 0.25;
    double p_bf = 0.25;

    double fg_alpha = 0, fg_beta = 0;
    double bg_alpha = 0, bg_beta = 0;
    double domain_score_cutoff = std::numeric_limits<double>::max();

    if (!params_in_file.empty()) { // read parameters file
      read_params_file(VERBOSE, params_in_file,
                       fg_alpha, fg_beta, bg_alpha, bg_beta,
                       p_fb, p_bf, domain_score_cutoff);
      max_iterations = 0;
    }
    else {
      const double n_reads = get_mean(begin(reads), end(reads));
      fg_alpha = 0.33*n_reads;
      fg_beta = 0.67*n_reads;
      bg_alpha = 0.67*n_reads;
      bg_beta = 0.33*n_reads;
    }

    if (max_iterations > 0)
      hmm.BaumWelchTraining(meth, reset_points, p_fb, p_bf,
                            fg_alpha, fg_beta, bg_alpha, bg_beta);

    // DECODE THE DOMAINS
    vector<bool> state_ids;
    vector<double> posteriors;
    hmm.PosteriorDecoding(meth, reset_points, p_fb, p_bf, fg_alpha, fg_beta,
                          bg_alpha, bg_beta, state_ids, posteriors);

    vector<double> domain_scores;
    get_domain_scores(state_ids, meth, reset_points, domain_scores);

    vector<double> random_scores;
    shuffle_cpgs(seed, hmm, meth, reset_points, p_fb, p_bf,
                 fg_alpha, fg_beta, bg_alpha, bg_beta, random_scores);

    vector<double> p_values;
    assign_p_values(random_scores, domain_scores, p_values);

    if (domain_score_cutoff == numeric_limits<double>::max())
      domain_score_cutoff = get_stepup_cutoff(p_values, 0.01);

    // write parameters if requested
    if (!params_out_file.empty())
      write_params_file(params_out_file, fg_alpha, fg_beta, bg_alpha, bg_beta,
                        p_fb, p_bf, domain_score_cutoff);

    vector<GenomicRegion> domains;
    build_domains(VERBOSE, cpgs, posteriors, reset_points, state_ids, domains);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    size_t good_hmr_count = 0;
    for (size_t i = 0; i < domains.size(); ++i)
      if (p_values[i] < domain_score_cutoff) {
        domains[i].set_name("HYPO" + to_string(good_hmr_count++));
        out << domains[i] << '\n';
      }

    if (!hypo_post_outfile.empty()) {
      if (VERBOSE)
        cerr << "[writing=" << hypo_post_outfile << "]" << endl;
      std::ofstream out(hypo_post_outfile);
      for (size_t i = 0; i < cpgs.size(); ++i) {
        GenomicRegion cpg(as_gen_rgn(cpgs[i]));
        cpg.set_name(format_cpg_meth_tag(meth[i]));
        cpg.set_score(posteriors[i]);
        out << cpg << '\n';
      }
    }

    if (!meth_post_outfile.empty()) {
      std::ofstream out(meth_post_outfile);
      if (VERBOSE)
        cerr << "[writing=" << meth_post_outfile << "]" << endl;
      for (size_t i = 0; i < cpgs.size(); ++i) {
        GenomicRegion cpg(as_gen_rgn(cpgs[i]));
        cpg.set_name(format_cpg_meth_tag(meth[i]));
        cpg.set_score(1.0 - posteriors[i]);
        out << cpg << '\n';
      }
    }
  }
  catch (runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
