/* Copyright (C) 2013 University of Southern California
 *                    Andrew D Smith and Jenny Qu
 *
 * Author: Andrew D. Smith and Jenny Qu
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
 */

/*****************************************
 *This program adopts the same TwoStateHMM functions as the hmr
 *program does.  First, the read counts at single CpGs are pooled in
 *non-overlapping bins.  The pooled counts for each bin is treated as
 *a single observation from a beta-binomial distribution.  Then the
 *two state HMM is built on these pooled counts.  If the bin size is
 *set to 1, this program would be the same as the hmr program.
 ****************************************/

#include <cmath>
#include <numeric>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <unistd.h>

#include <gsl/gsl_sf.h>

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
using std::make_pair;
using std::round;


static void
get_adjacent_distances(const vector<GenomicRegion> &pmds,
                       vector<size_t> &dists) {
  for (size_t i = 1; i < pmds.size(); ++i)
    if (pmds[i].same_chrom(pmds[i - 1]))
      dists.push_back(pmds[i].get_start() - pmds[i - 1].get_end());
}

static void
merge_nearby_pmd(const size_t max_merge, vector<GenomicRegion> &pmds) {
  size_t j = 0;
  for (size_t i = 1; i < pmds.size(); ++i) {
    if (pmds[j].same_chrom(pmds[i]) &&
        pmds[i].get_start() - pmds[j].get_end() <= max_merge)
      pmds[j].set_end(pmds[i].get_end());
    else pmds[++j] = pmds[i];
  }
  pmds.erase(pmds.begin() + j, pmds.end());
}

static bool
precedes(const string &chrom, const size_t position,
         const GenomicRegion &r) {
  return chrom < r.get_chrom() ||
                 (chrom == r.get_chrom() && position < r.get_start());
}

static bool
succeeds(const string &chrom, const size_t position,
         const GenomicRegion &r) {
  return r.get_chrom() < chrom ||
                         (chrom == r.get_chrom() && r.get_end() <= position);
}

static size_t
find_best_bound(const bool IS_RIGHT_BOUNDARY,
                const vector<pair<size_t, size_t> > &meth_tot,
                const vector<size_t> &positions) {

  size_t cumu_meth = 0, cumu_tot = 0;
  for (size_t i = 0; i < meth_tot.size(); ++i) {
    cumu_meth += meth_tot[i].first;
    cumu_tot += meth_tot[i].second;
  }

  size_t best_idx = 0;
  double best_score = -std::numeric_limits<double>::max();
  for (size_t i = 0; i < meth_tot.size(); ++i) {

    const size_t N_low = IS_RIGHT_BOUNDARY ?
      meth_tot[i].second : cumu_tot - meth_tot[i].second;
    const size_t k_low = IS_RIGHT_BOUNDARY ?
      meth_tot[i].first : cumu_meth - meth_tot[i].first;

    const size_t N_hi = IS_RIGHT_BOUNDARY ?
      cumu_tot - meth_tot[i].second : meth_tot[i].second;
    const size_t k_hi = IS_RIGHT_BOUNDARY ?
      cumu_meth - meth_tot[i].first : meth_tot[i].first;

    if (N_hi > 0 && N_low > 0) {

      const double p_low = static_cast<double>(k_low)/N_low;
      const double p_hi = static_cast<double>(k_hi)/N_hi;

      const double score =
        (gsl_sf_lnchoose(N_hi, k_hi) +
         k_hi*log(p_hi) + (N_hi - k_hi)*log(1.0 - p_hi)) +
        (gsl_sf_lnchoose(N_low, k_low) +
         k_low*log(p_low) + (N_low - k_low)*log(1.0 - p_low));

      if (p_hi > p_low && score > best_score) {
        best_idx = i;
        best_score = score;
      }
    }
  }
  return (best_score > -std::numeric_limits<double>::max()) ?
                 positions[best_idx] : -std::numeric_limits<size_t>::max();
}



static void
optimize_boundaries(const size_t bin_size,
                    const string &cpgs_file,
                    vector<GenomicRegion> &pmds, bool FORMAT) {

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ///// PULL OUT THE BOUNDARY POSITIONS
  /////
  vector<GenomicRegion> bounds;
  for (size_t i = 0; i < pmds.size(); ++i) {
    bounds.push_back(pmds[i]);
    if (pmds[i].get_start() > 0)
      bounds.back().set_start(pmds[i].get_start() - bin_size);
    bounds.back().set_end(pmds[i].get_start() + bin_size);

    bounds.push_back(pmds[i]);
    bounds.back().set_start(pmds[i].get_end() - bin_size);
    bounds.back().set_end(pmds[i].get_end() + bin_size);
  }

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ///// FIND THE EXACT CPG BOUNDARY SITES
  /////
  std::ifstream in(cpgs_file.c_str());

  vector<size_t> bound_site, positions;
  vector<pair<size_t, size_t> > meth_tot;
  string prevChrom;
  size_t bound_idx = 0;
  string chrom, site_name, strand;
  double meth_level;
  size_t position = 0ul, coverage = 0ul, n_meth = 0ul, n_unmeth = 0ul;

  while (methpipe::read_site(in, chrom, position,
                             strand, site_name, meth_level, coverage)
         && bound_idx < bounds.size()) {

    n_meth = round(meth_level*coverage);
    n_unmeth = coverage - n_meth;

    if (succeeds(chrom, position, bounds[bound_idx])) {
      // find the boundary CpG position
      bound_site.push_back(find_best_bound(bound_idx % 2, meth_tot, positions));
      positions.clear();
      meth_tot.clear();
      ++bound_idx;
    }
    // add CpG positions to the potential boundary set
    else if (!precedes(chrom, position, bounds[bound_idx])) {
      positions.push_back(position);
      meth_tot.push_back(make_pair(n_meth, n_meth + n_unmeth));
    }
    prevChrom.swap(chrom);
  }
  // find the boundary CpG position
  bound_site.push_back(find_best_bound(bound_idx % 2, meth_tot, positions));

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ///// NOW RE-ASSEMBLE THE GIANTS
  /////

  for (size_t i = 0; i < pmds.size(); ++i) {
    const size_t start_site = bound_site[2*i];
    if (start_site != -std::numeric_limits<size_t>::max())
      pmds[i].set_start(start_site);
    const size_t end_site = bound_site[2*i + 1];
    if (end_site != -std::numeric_limits<size_t>::max())
      pmds[i].set_end(end_site + 1);
  }

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ///// NOW MERGE THE PMDS THAT ARE TOO CLOSE
  /////

  vector<size_t> dists;
  get_adjacent_distances(pmds, dists);

  sort(dists.begin(), dists.end());

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  // NEED TO USE SOME RANDOMIZATION METHOD HERE TO FIGURE OUT THE
  // MERGING DISTANCE
  vector<pair<size_t, size_t> > dist_hist;
  size_t first = 0;
  for (size_t i = 1; i < dists.size(); ++i)
    if (dists[i] != dists[i - 1]) {
      dist_hist.push_back(make_pair(dists[i - 1], i - first));
      first = i;
    }

  merge_nearby_pmd(2*bin_size, pmds);
}



double
get_score_cutoff_for_fdr(const vector<double> &scores, const double fdr) {
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
get_domain_scores(const vector<bool> &classes,
                  const vector<pair<double, double> > &meth,
                  const vector<size_t> &reset_points,
                  vector<double> &scores) {
  static const bool CLASS_ID = true;
  size_t n_cpgs = 0, reset_idx = 1;
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
      score += 1.0 - (meth[i].first/(meth[i].first + meth[i].second));
      ++n_cpgs;
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
    }
    else if (in_domain) {
      in_domain = false;
      domains.back().set_end(prev_end);
      domains.back().set_score(n_cpgs);
      n_cpgs = 0;
      score = 0;
    }
    prev_end = cpgs[i].get_end();
  }
  if(in_domain){
    domains.back().set_end(prev_end);
    domains.back().set_score(n_cpgs);
  }
}



template <class T, class U> static void
separate_regions(const bool VERBOSE, const size_t desert_size,
                 vector<SimpleGenomicRegion> &cpgs,
                 vector<T> &meth, vector<U> &reads,
                 vector<size_t> &reset_points) {
  if (VERBOSE)
    cerr << "[SEPARATING BY CPG DESERT]" << endl;
  // eliminate the zero-read cpgs
  size_t j = 0;
  for (size_t i = 0; i < cpgs.size(); ++i)
    if (reads[i] > 0) {
      cpgs[j] = cpgs[i];
      meth[j] = meth[i];
      reads[j] = reads[i];
      ++j;
    }
  cpgs.erase(cpgs.begin() + j, cpgs.end());
  meth.erase(meth.begin() + j, meth.end());
  reads.erase(reads.begin() + j, reads.end());

  // segregate cpgs
  size_t prev_cpg = 0;
  for (size_t i = 0; i < cpgs.size(); ++i) {
    const size_t dist = (i > 0 && cpgs[i].same_chrom(cpgs[i - 1])) ?
      cpgs[i].get_start() - prev_cpg : numeric_limits<size_t>::max();
    if (dist > desert_size)
      reset_points.push_back(i);
    prev_cpg = cpgs[i].get_start();
  }
  reset_points.push_back(cpgs.size());
  if (VERBOSE)
    cerr << "CPGS RETAINED: " << cpgs.size() << endl
         << "DESERTS REMOVED: " << reset_points.size() - 2 << endl << endl;
}



static void
shuffle_cpgs(const TwoStateHMMB &hmm,
             vector<pair<double, double> > meth,
             vector<size_t> reset_points,
             const vector<double> &start_trans,
             const vector<vector<double> > &trans,
             const vector<double> &end_trans,
             const double fg_alpha, const double fg_beta,
             const double bg_alpha, const double bg_beta,
             vector<double> &domain_scores) {
  srand(time(0) + getpid());

  static const size_t n_shuffles = 100;
  for (size_t i = 0; i < n_shuffles; ++i) {
    random_shuffle(meth.begin(), meth.end());
    vector<bool> classes;
    vector<double> scores;
    hmm.PosteriorDecoding(meth, reset_points, start_trans, trans,
                          end_trans, fg_alpha, fg_beta, bg_alpha,
                          bg_beta, classes, scores);
    get_domain_scores(classes, meth, reset_points, domain_scores);
  }
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
                 double &score_cutoff_for_fdr) {
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
     >> jnk >> score_cutoff_for_fdr;
  if (VERBOSE)
    cerr << "FG_ALPHA\t" << fg_alpha << endl
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
         << "SCORE_CUTOFF_FOR_FDR\t" << score_cutoff_for_fdr << endl;
}



static void
write_params_file(const string &outfile,
                  const double fg_alpha,
                  const double fg_beta,
                  const double bg_alpha,
                  const double bg_beta,
                  const vector<double> &start_trans,
                  const vector<vector<double> > &trans,
                  const vector<double> &end_trans) {

  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out.precision(30);
  out << "FG_ALPHA\t" << fg_alpha << endl
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
    ;
}



static void
load_intervals(const size_t bin_size,
               const string &cpgs_file,
               vector<SimpleGenomicRegion> &intervals,
               vector<pair<double, double> > &meth,
               vector<size_t> &reads, bool FORMAT) {

  std::ifstream in(cpgs_file.c_str());

  string curr_chrom;
  size_t prev_pos = 0ul, curr_pos = 0ul;
  size_t n_meth_bin = 0ul, n_unmeth_bin = 0ul;
  string chrom, site_name, strand;
  double meth_level;
  size_t position = 0ul, coverage = 0ul, n_meth = 0ul, n_unmeth = 0ul;

  while (methpipe::read_site(in, chrom, position,
                             strand, site_name, meth_level, coverage)) {

    n_meth = round(meth_level * coverage);
    n_unmeth = coverage - n_meth;

    // no range, or no chrom
    if (curr_chrom != chrom) {
      if (!curr_chrom.empty()) {
        if(chrom < curr_chrom)
          throw SMITHLABException("CpGs not sorted in file \"" + cpgs_file + "\"");
        intervals.push_back(SimpleGenomicRegion(chrom, curr_pos, prev_pos + 1));
        reads.push_back(n_meth_bin + n_unmeth_bin);
        meth.push_back(make_pair(n_meth_bin, n_unmeth_bin));
      }
      curr_chrom = chrom;
      curr_pos = position;
      n_meth_bin = 0ul;
      n_unmeth_bin = 0ul;
    }
    else if (curr_pos > position){
      throw SMITHLABException("CpGs not sorted in file \"" + cpgs_file + "\"");
    }
    else if (position > curr_pos + bin_size) {
      intervals.push_back(SimpleGenomicRegion(curr_chrom, curr_pos,
                                              curr_pos + bin_size));
      reads.push_back(n_meth_bin + n_unmeth_bin);
      meth.push_back(make_pair(n_meth_bin, n_unmeth_bin));

      n_meth_bin = 0ul;
      n_unmeth_bin = 0ul;
      curr_pos += bin_size;
      while (curr_pos + bin_size < position) {
        intervals.push_back(SimpleGenomicRegion(curr_chrom, curr_pos,
                                                curr_pos + bin_size));
        reads.push_back(0);
        meth.push_back(make_pair(0.0, 0.0));

        curr_pos += bin_size;
      }
    }
    n_meth_bin += n_meth;
    n_unmeth_bin += n_unmeth;
    prev_pos = position;
  }
  if (!curr_chrom.empty()) {
    intervals.push_back(SimpleGenomicRegion(curr_chrom, curr_pos, prev_pos + 1));
    reads.push_back(n_meth_bin + n_unmeth_bin);
    meth.push_back(make_pair(n_meth_bin, n_unmeth_bin));
  }
}



int
main(int argc, const char **argv) {

  try {

    string outfile;

    size_t desert_size = 20000;
    size_t max_iterations = 10;

    // run mode flags
    bool VERBOSE = false;

    // corrections for small values (not parameters):
    double tolerance = 1e-10;
    double min_prob  = 1e-10;

    string params_in_file;
    string params_out_file;

    double fdr_cutoff = 0.01;

    size_t bin_size = 1000;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "identify PMDs in a methylome",
                           "<cpg-meth-file>");
    opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("desert", 'd', "max allowed unmapped region size",
                      false, desert_size);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations);
    opt_parse.add_opt("bin", 'b', "bin size", false, bin_size);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("fdr", '\0', "fdr cutoff", false, fdr_cutoff);
    opt_parse.add_opt("params-in", 'P', "read HMM parameters from file",
                      false, params_in_file);
    opt_parse.add_opt("params-out", 'p', "write HMM parameters to file",
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string cpgs_file = leftover_args.front();
    bool FORMAT;
    FORMAT = methpipe::is_methpipe_file_single(cpgs_file);
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "[READING CPGS AND METH PROPS]" << endl;
    vector<pair<double, double> > meth;
    vector<size_t> reads;
    vector<SimpleGenomicRegion> cpgs;
    load_intervals(bin_size, cpgs_file, cpgs, meth, reads, FORMAT);

    if (VERBOSE)
      cerr << "CPG SITES LOADED: " << cpgs.size() << endl;

    // separate the regions by chrom and by desert, and eliminate
    // those isolated CpGs
    vector<size_t> reset_points;
    separate_regions(VERBOSE, desert_size, cpgs, meth, reads, reset_points);

    vector<double> start_trans(2, 0.5), end_trans(2, 1e-10);
    vector<vector<double> > trans(2, vector<double>(2, 0.25));
    trans[0][0] = trans[1][1] = 0.75;

    const TwoStateHMMB hmm(min_prob, tolerance, max_iterations, VERBOSE);

    double fg_alpha = 0, fg_beta = 0;
    double bg_alpha = 0, bg_beta = 0;
    double score_cutoff_for_fdr = std::numeric_limits<double>::max();

    if (!params_in_file.empty()) {
      // READ THE PARAMETERS FILE
      read_params_file(VERBOSE, params_in_file,
                       fg_alpha, fg_beta, bg_alpha, bg_beta,
                       start_trans, trans, end_trans, score_cutoff_for_fdr);
    }
    else {
      const double n_reads =
        accumulate(reads.begin(), reads.end(), 0.0)/reads.size();
      fg_alpha = 0.33*n_reads;
      fg_beta = 0.67*n_reads;
      bg_alpha = 0.67*n_reads;
      bg_beta = 0.33*n_reads;
    }

    if (max_iterations > 0)
      hmm.BaumWelchTraining(meth, reset_points, start_trans, trans,
                            end_trans, fg_alpha, fg_beta, bg_alpha, bg_beta);

    if (!params_out_file.empty()) {
      // WRITE ALL THE HMM PARAMETERS:
      write_params_file(params_out_file, fg_alpha, fg_beta, bg_alpha, bg_beta,
                        start_trans, trans, end_trans);
    }

    /***********************************
     * STEP 5: DECODE THE DOMAINS
     */
    vector<bool> classes;
    vector<double> scores;
    hmm.PosteriorDecoding(meth, reset_points, start_trans, trans,
                          end_trans, fg_alpha, fg_beta, bg_alpha,
                          bg_beta, classes, scores);

    vector<double> domain_scores;
    get_domain_scores(classes, meth, reset_points, domain_scores);

    if (VERBOSE)
      cerr << "[RANDOMIZING SCORES FOR FDR]" << endl;

    vector<double> random_scores;
    shuffle_cpgs(hmm, meth, reset_points, start_trans, trans, end_trans,
                 fg_alpha, fg_beta, bg_alpha, bg_beta, random_scores);

    vector<double> p_values;
    assign_p_values(random_scores, domain_scores, p_values);

    if (score_cutoff_for_fdr == numeric_limits<double>::max())
      score_cutoff_for_fdr = get_score_cutoff_for_fdr(p_values, fdr_cutoff);

    if (!params_out_file.empty()) {
      std::ofstream out(params_out_file.c_str(), std::ios::app);
      out << "SCORE_CUTOFF_FOR_FDR\t"
          << std::setprecision(30) << score_cutoff_for_fdr << endl;
      out.close();
    }

    vector<GenomicRegion> domains;
    build_domains(VERBOSE, cpgs, scores, reset_points, classes, domains);

    size_t good_hmr_count = 0;
    vector<GenomicRegion> good_domains;
    for (size_t i = 0; i < domains.size(); ++i) {
      if (p_values[i] < score_cutoff_for_fdr) {
        good_domains.push_back(domains[i]);
        good_domains.back().set_name("PMD" + smithlab::toa(good_hmr_count++));
      }
    }

    optimize_boundaries(bin_size, cpgs_file, good_domains, FORMAT);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    copy(good_domains.begin(), good_domains.end(),
         std::ostream_iterator<GenomicRegion>(out, "\n"));

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
