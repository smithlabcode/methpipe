/*
 * Copyright (C) 2011 University of Southern California
 *                    Andrew D Smith and Qiang Song
 * Author: Qiang Song and Andrew D. Smith
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

#include <cmath>
#include <fstream>
#include <algorithm>
#include <utility>
#include <numeric>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "SplitDistro.hpp"
#include "ReadCounts.hpp"
#include "SelectBinSize.hpp"
#include "EvaluateBoundaries.hpp"
#include "LoadReadsByRegion.hpp"
#include "OptionParser.hpp"

#include "ModelParams.hpp"

#include "TwoStateScaleSplitResolveMixture.hpp"
#include "TwoStateScaleSplitHMM.hpp"

#include "ThreeStateScaleSplitResolveMixture.hpp"
#include "ThreeStateScaleSplitHMM.hpp"

#include "rseg_utils.hpp"

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

// functions for two-state modes

void
output_boundaries(const vector<double> &tmp_read_bins,
                  const vector<double> &scales,
                  const vector<size_t> &reset_points,
                  const vector<vector<SimpleGenomicRegion> > &bin_bounds,
                  const TwoStateScaleSplitHMM &hmm,
                  const vector<SplitDistro> &distros,
                  const vector<double> &start_trans,
                  const vector<vector<double> > &trans,
                  const vector<double> &end_trans,
                  const vector<bool> &tmp_classes,
                  const string &boundary_file,
                  const string &boundary_score_file,
                  const bool VERBOSE) {

  vector<vector<vector<double> > > post_trans_scores;
  hmm.TransitionPosteriors(tmp_read_bins, scales, reset_points,
                           start_trans, trans, end_trans,
                           distros.front(), distros.back(),
                           post_trans_scores);

  vector<double>& f_to_f_scores(post_trans_scores[0][0]);
  vector<double>& f_to_b_scores(post_trans_scores[0][1]);
  vector<double>& b_to_f_scores(post_trans_scores[1][0]);
  vector<double>& b_to_b_scores(post_trans_scores[1][1]);

  vector<double> tmp_boundary_scores(f_to_b_scores.size());
  vector<int> transitions;
  for (size_t i = 0; i < tmp_classes.size(); ++i)
    if (i == 0 || tmp_classes[i] != tmp_classes[i - 1])
      transitions.push_back(i);

  transitions.push_back(tmp_classes.size());
  size_t j = 0;

  for (int i = 0; static_cast<size_t>(i) < f_to_b_scores.size(); ++i)
    {
      static const size_t sz = tmp_classes.size();

      if (abs(i - transitions[j]) > abs(i - transitions[j + 1]))
        ++j;
      const size_t trn = transitions[j];
      if (trn == 0 || trn == sz)
        tmp_boundary_scores[i] = 0;
      else
        tmp_boundary_scores[i] = (tmp_classes[trn]) ?
          b_to_f_scores[i] : f_to_b_scores[i];
    }

  vector<vector<double> > boundary_scores;
  expand_bins(tmp_boundary_scores, reset_points, boundary_scores);

  //// generate control sample
  const size_t rand_sample_size = std::min(size_t(150000), tmp_read_bins.size());
  vector<pair<double, double> > vals_control(rand_sample_size);

  for (size_t i = 0; i < rand_sample_size; ++i)
    vals_control[i] = make_pair(tmp_read_bins[i], scales[i]);

  vector<size_t> reset_points_control;
  for (size_t i = 0; i < reset_points.size(); ++i)
    if (reset_points[i] < rand_sample_size)
      reset_points_control.push_back(reset_points[i]);
    else break;
  reset_points_control.push_back(rand_sample_size);

  for (size_t i = 0; i < reset_points_control.size() - 1; ++i)
    std::random_shuffle(vals_control.begin() + reset_points_control[i],
                        vals_control.begin() + reset_points_control[i + 1]);

  vector<double> read_counts_control(rand_sample_size),
    scales_control(rand_sample_size);

  for (size_t i = 0; i < rand_sample_size; ++i) {
    read_counts_control[i] = vals_control[i].first;
    scales_control[i] = vals_control[i].second;
  }
  vals_control.clear();

  vector<vector<vector<double> > > post_trans_scores_control;
  hmm.TransitionPosteriors(read_counts_control, scales_control,
                           reset_points_control,
                           start_trans, trans, end_trans,
                           distros.front(), distros.back(),
                           post_trans_scores_control);

  vector<double>& f_to_b_scores_control = post_trans_scores_control[0][1];
  vector<double>& b_to_f_scores_control = post_trans_scores_control[1][0];

  vector<double> boundary_scores_control(f_to_b_scores_control.size(), 0);

  for (size_t i = 0; i < boundary_scores_control.size(); ++i)
    boundary_scores_control[i] = std::max(f_to_b_scores_control[i],
                                          b_to_f_scores_control[i]);

  std::sort(boundary_scores_control.begin(), boundary_scores_control.end());
  double fdr = 0.05;

  double cutoff = boundary_scores_control[static_cast<size_t>(boundary_scores_control.size() * (1 - fdr))];

  read_counts_control.clear();
  scales_control.clear();
  b_to_f_scores_control.clear();
  f_to_b_scores_control.clear();
  boundary_scores_control.clear();
  post_trans_scores_control.clear();

  //// finish generating control sample

  vector<GenomicRegion> boundaries;
  BoundEval be(1, 1);
  be.evaluate(bin_bounds, reset_points, tmp_classes, tmp_boundary_scores,
              f_to_f_scores, f_to_b_scores,
              b_to_f_scores, b_to_b_scores,
              cutoff, true, boundaries);

  // write result files
  WriteBEDFile(boundary_file, boundaries);
  if (VERBOSE)
    cout << "Boundary file: " + boundary_file << std::endl;

  if (!boundary_score_file.empty() && boundary_score_file != "None")
    {
      write_wigfile(boundary_scores, bin_bounds, boundary_score_file);
      if (VERBOSE)
        cout << "Boundary score file: " + boundary_score_file << std::endl;
    }
}

void
output_domains(
               const vector<double> &tmp_read_bins,
               const vector<double> &tmp_scales,
               const vector<size_t> &reset_points,
               const vector<vector<SimpleGenomicRegion> > &bin_bounds,
               const TwoStateScaleSplitHMM &hmm,
               const vector<SplitDistro> &distros,
               const vector<double> &start_trans,
               const vector<vector<double> > &trans,
               const vector<double> &end_trans,
               const vector<bool> &tmp_classes,
               const double posterior_cutoff,
               const size_t undef_region_cutoff,
               const double cdf_cutoff,
               const string &domain_file,
               const string &posterior_score_file,
               const bool VERBOSE)
{
  // Obtain the scores for the current domain class
  vector<double> tmp_scores;
  hmm.PosteriorScores(tmp_read_bins, tmp_scales, reset_points,
                      start_trans, trans, end_trans,
                      distros.front(), distros.back(), tmp_classes, tmp_scores);

  vector<vector<double> > read_bins, scores, scales;
  vector<vector<bool> > classes;
  expand_bins(tmp_read_bins, reset_points, read_bins);
  expand_bins(tmp_scores, reset_points, scores);
  expand_bins(tmp_scales, reset_points, scales);
  expand_bins(tmp_classes, reset_points, classes);

  vector<vector<GenomicRegion> > domains;
  build_domains(bin_bounds, classes, scores, posterior_cutoff,
                domains, undef_region_cutoff);
  pick_domains(bin_bounds, read_bins, scales, distros, domains, cdf_cutoff);

  // output domains
  write_bed_file(domains, domain_file);
  if (VERBOSE)
    cout << "Domains file: " + domain_file << std::endl;

  if (!posterior_score_file.empty() && posterior_score_file != "None")
    {
      size_t k = 0;
      for (size_t i = 0; i < scores.size(); ++i)
        for (size_t j = 0; j < scores[i].size(); ++j)
          {
            scores[i][j] == tmp_classes[k] ? scores[i][j] : 1 - scores[i][j];
            ++k;
          }
      write_wigfile(scores, bin_bounds, posterior_score_file);
      if (VERBOSE)
        cout << "Bin score file: " + posterior_score_file << std::endl;
    }
}

// end of functions for two-state modes


// for three-state mode
void
output_boundaries(const vector<double> &tmp_read_bins,
                  const vector<double> &scales,
                  const vector<size_t> &reset_points,
                  const vector<vector<SimpleGenomicRegion> > &bin_bounds,
                  const ThreeStateScaleSplitHMM &hmm,
                  const vector<SplitDistro> &distros,
                  const vector<double> &start_trans,
                  const vector<vector<double> > &trans,
                  const vector<double> &end_trans,
                  const vector<size_t> &classes,
                  const string &boundary_file,
                  const string &boundary_score_file,
                  const bool VERBOSE)
{

  const size_t NUM_OF_STATES = 3;

  vector<vector<vector<double> > > post_trans_scores;
  hmm.TransitionPosteriors(tmp_read_bins, scales, reset_points,
                           start_trans, trans, end_trans,
                           distros.front(), distros[1], distros.back(),
                           post_trans_scores);

  vector<size_t> change_points;
  for (size_t i = 0; i < classes.size(); ++i)
    if (i == 0 || classes[i] != classes[i - 1])
      change_points.push_back(i);
  change_points.push_back(classes.size());


  vector<double> tmp_boundary_scores(tmp_read_bins.size());
  size_t j = 0;

  for (int i = 0; static_cast<size_t>(i) < tmp_boundary_scores.size(); ++i)
    {
      static const size_t sz = tmp_boundary_scores.size();

      if (abs(i - change_points[j]) > abs(i - change_points[j + 1]))
        ++j;

      const size_t trn = change_points[j];

      if (trn == 0 || trn == sz) // indicates a starting domain
        tmp_boundary_scores[i] = 0;
      else
        tmp_boundary_scores[i] = post_trans_scores[classes[trn - 1]][classes[trn]][i];
    }

  vector< vector<double> > boundary_scores;
  expand_bins(tmp_boundary_scores, reset_points, boundary_scores);


  ///// genereate contronl sample
  const size_t rand_sample_size = std::min(size_t(150000), tmp_read_bins.size());
  vector<pair<double, double> > vals_control(rand_sample_size);

  for (size_t i = 0; i < rand_sample_size; ++i)
    vals_control[i] = make_pair(tmp_read_bins[i], scales[i]);

  vector<size_t> reset_points_control;
  for (size_t i = 0; i < reset_points.size(); ++i)
    if (reset_points[i] < rand_sample_size)
      reset_points_control.push_back(reset_points[i]);
    else
      break;
  reset_points_control.push_back(rand_sample_size);

  for (size_t i = 0; i < reset_points_control.size() - 1; ++i)
    std::random_shuffle(vals_control.begin() + reset_points_control[i],
                        vals_control.begin() + reset_points_control[i + 1]);

  vector<double> read_counts_control(rand_sample_size),
    scales_control(rand_sample_size);

  for (size_t i = 0; i < rand_sample_size; ++i)
    {
      read_counts_control[i] = vals_control[i].first;
      scales_control[i] = vals_control[i].second;
    }
  vals_control.clear();


  vector<vector<vector<double> > > post_trans_scores_control;
  hmm.TransitionPosteriors(read_counts_control, scales_control, reset_points_control,
                           start_trans, trans, end_trans,
                           distros.front(), distros[1], distros.back(),
                           post_trans_scores_control);

  vector<double> boundary_scores_control(read_counts_control.size(), 0);
  for (size_t i = 0; i < boundary_scores_control.size(); ++i)
    for (size_t j = 0; j < NUM_OF_STATES; ++j)
      for (size_t k = 0; k < NUM_OF_STATES; ++k)
        if (j != k &&
            post_trans_scores_control[j][k][i] > boundary_scores_control[i])
          boundary_scores_control[i] = post_trans_scores_control[j][k][i];

  std::sort(boundary_scores_control.begin(), boundary_scores_control.end());
  double fdr = 0.05;
  double cutoff = boundary_scores_control[
                                          static_cast<size_t>(boundary_scores_control.size() * (1 - fdr))];

  read_counts_control.clear();
  scales_control.clear();
  reset_points_control.clear();
  post_trans_scores_control.clear();
  boundary_scores_control.clear();
  //// finish generating control sample

  vector<GenomicRegion> boundaries;
  BoundEval be(1, 1);
  be.evaluate(bin_bounds, reset_points, classes, tmp_boundary_scores,
              post_trans_scores, cutoff, boundaries);

  // write result files
  WriteBEDFile(boundary_file, boundaries);
  if (VERBOSE)
    cout << "Boundary file: " + boundary_file << std::endl;

  if (!boundary_score_file.empty() && boundary_score_file != "None")
    {
      write_wigfile(boundary_scores, bin_bounds, boundary_score_file);
      if (VERBOSE)
        cout << "Boundary score file: " + boundary_score_file << std::endl;
    }
}

void
output_domains(
               const vector<double> &tmp_read_bins,
               const vector<double> &tmp_scales,
               const vector<size_t> &reset_points,
               const vector<vector<SimpleGenomicRegion> > &bin_bounds,
               const ThreeStateScaleSplitHMM &hmm,
               const vector<SplitDistro> &distros,
               const vector<double> &start_trans,
               const vector<vector<double> > &trans,
               const vector<double> &end_trans,
               const vector<size_t> &classes,
               const double posterior_cutoff,
               const size_t undef_region_cutoff,
               const double cdf_cutoff,
               const string &domain_file,
               const string &posterior_score_file,
               const bool VERBOSE)
{
  // Obtain the scores for the current domain class
  vector<double> tmp_scores;
  hmm.PosteriorScores(tmp_read_bins, tmp_scales, reset_points,
                      start_trans, trans, end_trans,
                      distros.front(), distros[1], distros.back(),
                      classes, tmp_scores);

  vector<vector<double> > read_bins, scores, scales;
  vector<vector<size_t> > expanded_classes;
  expand_bins(tmp_read_bins, reset_points, read_bins);
  expand_bins(tmp_scores, reset_points, scores);
  expand_bins(tmp_scales, reset_points, scales);
  expand_bins(classes, reset_points, expanded_classes);

  vector<vector<GenomicRegion> > domains;
  build_domains(bin_bounds, expanded_classes, scores, posterior_cutoff,
                domains, undef_region_cutoff);
  pick_domains_3s(bin_bounds, read_bins, scales,
                  distros, domains, cdf_cutoff);

  // output domains
  write_bed_file(domains, domain_file);
  if (VERBOSE)
    cout << "Domains file: " + domain_file << std::endl;

  if (!posterior_score_file.empty() && posterior_score_file != "None")
    {
      vector<double> fg_scores, bg_scores;
      hmm.PosteriorScores(tmp_read_bins, tmp_scales, reset_points,
                          start_trans, trans, end_trans,
                          distros.front(), distros[1], distros.back(),
                          fg_scores, bg_scores);
      write_wigfile(fg_scores, bg_scores, bin_bounds, posterior_score_file);
      if (VERBOSE)
        cout << "Bin score file: " + posterior_score_file << std::endl;
    }
}

// end of functions for three-state modes
int
main(int argc, const char **argv)  {
  try {

    string  deads_file, chroms_file;
    string in_param_file;
    string out_param_file;

    // expected size of a domain
    double fg_size = 20000;

    // flags
    bool USE_POSTERIOR = false;
    bool REMOVE_JACKPOT = true;
    bool VERBOSE = false;
    bool BAM_FORMAT = false;

    string domain_file = "/dev/stdout";
    string posterior_score_file = "";
    string boundary_file = "";
    string boundary_score_file = "";
    string read_counts_file = "";

    // mode
    int mode = 2;
    const int TEST_CONTROL_MODE = 2;
    const int TEST_TEST_MODE = 3;

    // name of emission distributions
    string fg_name("nbdiff"), bg_name("nbdiff");

    size_t desert_size = 20000;
    size_t bin_size_step = 50;
    size_t bin_size = 0;
    size_t FRAGMENT_LEN = 0;
    bool waterman = false;
    bool hideaki = false;
    bool hideaki_emp = false;
    bool smooth = true;
    size_t max_iterations = 20;
    size_t training_size = 0;
    double tolerance = 1e-20;
    double min_prob = 1e-20;

    double max_dead_proportion = 0.5;

    double posterior_cutoff = 0.5;
    size_t undef_region_cutoff = 3000;
    double cdf_cutoff = 0.1;

    // Determines how many iterations are used during the initialization
    // phase to find good starting values for the HMM
    const size_t MAX_INITIALIZATION_ITR = 5;

    ////////////////////// COMMAND LINE OPTIONS /////////////////////////
    OptionParser opt_parse(strip_path(argv[0]),
                           "segment the genome according to differential "
                           "mapped read density", "<mapped-read-locations-A> "
                           "<mapped-read-locations-B>");
    opt_parse.add_opt("out", 'o', "domain output file", false, domain_file);
    opt_parse.add_opt("score", '\0', "Posterior scores file",
                      false, posterior_score_file);
    opt_parse.add_opt("readcount", '\0', "readcounts file",
                      false, read_counts_file);
    opt_parse.add_opt("boundary", '\0', "domain boundary file",
                      false, boundary_file);
    opt_parse.add_opt("boundary-score", '\0', "boundary transition scores file",
                      false, boundary_score_file);
    opt_parse.add_opt("chrom", 'c', "file with chromosome sizes (BED format)",
                      true, chroms_file);
    opt_parse.add_opt("deadzones", 'd', "file of deadzones (BED format)",
                      false, deads_file);
    opt_parse.add_opt("bam", 'B', "Input reads file is BAM format",
                      false, BAM_FORMAT);
    opt_parse.add_opt("param-in", '\0', "Input parameters file",
                      false, in_param_file);
    opt_parse.add_opt("param-out", '\0', "Output parameters file",
                      false, out_param_file);
    opt_parse.add_opt("mode", 'm', "running mode 2:test-control; 3: test-test",
                      false, mode);
    opt_parse.add_opt("maxitr", 'i', "maximum iterations for training",
                      false, max_iterations);
    opt_parse.add_opt("bin-size", 'b', "bin size (default: based on data)",
                      false, bin_size);
    opt_parse.add_opt("bin-step", '\0',
                      "minimum bin size (default: " + toa(bin_size_step) + ")",
                      false, bin_size_step);
    opt_parse.add_opt("duplicates", '\0', "keep duplicate reads",
                      false, REMOVE_JACKPOT);
    opt_parse.add_opt("fragment_length", '\0',
                      "Extend reads to fragment length (default not to extend)",
                      false, FRAGMENT_LEN);
    opt_parse.add_opt("Waterman", '\0', "use Waterman's method for bin size",
                      false, waterman);
    opt_parse.add_opt("Hideaki", '\0', "use Hideaki's method for bin size",
                      false, hideaki);
    opt_parse.add_opt("Hideaki-emp", '\0', "use Hideaki's empirical method (default)",
                      false, hideaki_emp);
    opt_parse.add_opt("smooth", '\0', "Indicate whether the rate curve is assumed smooth",
                      false, smooth);
    opt_parse.add_opt("max-dead", '\0',
                      "max deadzone proportion for retained bins",
                      false, max_dead_proportion);
    opt_parse.add_opt("domain-size", 's', "expected domain size "
                      "(default: " + toa(fg_size) + ")",
                      false, fg_size);
    opt_parse.add_opt("desert", 'S', "desert size "
                      "(default: " + toa(desert_size) + ")",
                      false, desert_size);
    opt_parse.add_opt("fg", 'F', "foreground emission distribution", false, fg_name);
    opt_parse.add_opt("bg", 'B', "background emission distribution", false, bg_name);
    opt_parse.add_opt("training-size", '\0',
                      "Max number of data points for training (default: all)",
                      false, training_size);
    opt_parse.add_opt("posterior", 'P', "use posterior decoding "
                      "(default: Viterbi)", false, USE_POSTERIOR);
    opt_parse.add_opt("posterior-cutoff", '\0',
                      "Posterior threshold for signigicant bins",
                      false, posterior_cutoff);
    opt_parse.add_opt("undefined", '\0', "min size of unmappable region",
                      false, undef_region_cutoff);
    opt_parse.add_opt("cutoff", '\0', "cutoff in cdf for identified domains",
                      false, cdf_cutoff);
    opt_parse.add_opt("verbose", 'v', "print more run information",
                      false, VERBOSE);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);

    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
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

    if (leftover_args.size() < 2) {
      cerr << "Need two reads files" << endl;
      return EXIT_SUCCESS;
    }

    const string reads_file_a = leftover_args[0];
    const string reads_file_b = leftover_args[1];


    /**********************************************************************/

    /***********************************
     * STEP 1: READ IN THE DATA
     */
    vector<SimpleGenomicRegion> bin_boundaries;
    vector<double> read_bins_a;
    vector<double> read_bins_b;
    vector<double> scales;
    vector<size_t> reset_points;
    LoadReadsByRegion(VERBOSE, chroms_file, reads_file_a, reads_file_b,
                      deads_file, bin_size_step, bin_boundaries, read_bins_a,
                      read_bins_b, scales, reset_points,
                      FRAGMENT_LEN, BAM_FORMAT, REMOVE_JACKPOT);

    if (VERBOSE)
      cout << "[SELECTING BIN SIZE] ";
    if (bin_size == 0) {
      if (hideaki)
        bin_size = select_bin_size_hideaki(
                                           read_bins_a, scales, bin_size_step, smooth);
      else if (waterman)
        bin_size = select_bin_size_waterman(
                                            read_bins_a, scales, bin_size_step, smooth);
      else
        bin_size = select_bin_size_hideaki_emp(
                                               read_bins_a, scales, reset_points, bin_size_step,
                                               max_dead_proportion);
    }
    if (VERBOSE) cout << "bin size =  " << bin_size << endl;


    /***********************************
     * STEP 2: BIN THE READS
     */
    AdjustBinSize(bin_boundaries, read_bins_a, read_bins_b, scales,
                  reset_points, bin_size_step, bin_size);
    RemoveDeserts(bin_boundaries, read_bins_a, read_bins_b, scales,
                  reset_points, bin_size, desert_size, max_dead_proportion);

    vector<double> read_bins(read_bins_a.size());
    const double max_count = bin_size;
    for (size_t i = 0; i < read_bins.size(); ++i) {
      read_bins_a[i] = min(read_bins_a[i], max_count);
      read_bins_b[i] = min(read_bins_b[i], max_count);
      read_bins[i] = read_bins_a[i] - read_bins_b[i];
    }

    vector<vector<SimpleGenomicRegion> > bin_boundaries_folded;
    expand_bins(bin_boundaries, reset_points, bin_boundaries_folded);

    // mode specific code
    if (mode == TEST_CONTROL_MODE) {
      /***********************************
       * STEP 3: ESTIMATE EMISSION PARAMS
       */

      const TwoStateScaleSplitHMM hmm(min_prob, tolerance, max_iterations,
                                      VERBOSE);
      size_t state_num = 2;
      vector<SplitDistro> distros;
      vector<vector<double> > trans;
      vector<double> start_trans, end_trans;

      if (!in_param_file.empty())
        read_param_file(in_param_file, state_num,
                        start_trans, trans, end_trans, distros);
      else {
        if (VERBOSE)
          cout << "[ESTIMATING PARAMETERS]" << endl;

        fg_size = (fg_size > 0) ? fg_size : 20000;

        training_size =
          (training_size == 0) ? read_bins.size() : training_size;

        vector<double> read_bins_sample, read_bins_a_sample,
          read_bins_b_sample, scales_sample;
        vector<size_t> reset_points_sample;
        pick_training_sample(read_bins, read_bins_a, read_bins_b, scales,
                             reset_points,
                             training_size, read_bins_sample, read_bins_a_sample,
                             read_bins_b_sample, scales_sample, reset_points_sample);

        distros.push_back(SplitDistro(fg_name));
        distros.push_back(SplitDistro(bg_name));

        double mixing = 0;
        TwoStateSplitResolveMixture(read_bins_sample, read_bins_a_sample,
                                    read_bins_b_sample, scales_sample,
                                    MAX_INITIALIZATION_ITR, tolerance, VERBOSE,
                                    distros.front(), distros.back(), mixing);

        set_transitions(bin_size, fg_size, mixing, VERBOSE,
                        start_trans, trans, end_trans);

        hmm.BaumWelchTraining(read_bins_sample, read_bins_a_sample,
                              read_bins_b_sample, scales_sample, reset_points_sample,
                              start_trans, trans, end_trans, distros.front(), distros.back());

        clear_training_sample(read_bins_sample, read_bins_a_sample,
                              read_bins_b_sample, scales_sample, reset_points_sample);
      }

      if (!out_param_file.empty())
        write_param_file(out_param_file, state_num,
                         start_trans, trans, end_trans, distros);

      if (VERBOSE)
        report_final_values(distros, start_trans, trans, end_trans);

      /***********************************
       * STEP 5: DECODE THE DOMAINS
       */

      vector<bool> classes;
      vector<double> scores;
      if (USE_POSTERIOR)
        hmm.PosteriorDecoding(read_bins, scales, reset_points,
                              start_trans, trans, end_trans,
                              distros.front(), distros.back(), classes, scores);
      else
        hmm.ViterbiDecoding(read_bins, scales, reset_points,
                            start_trans, trans, end_trans,
                            distros.front(), distros.back(), classes);

      /***********************************
       * STEP 6: WRITE THE RESULTS
       */

      // make sure the output dir is valid
      output_domains(read_bins, scales, reset_points, bin_boundaries_folded,
                     hmm, distros, start_trans, trans, end_trans, classes,
                     posterior_cutoff, undef_region_cutoff, cdf_cutoff,
                     domain_file, posterior_score_file, VERBOSE);
      if (!boundary_file.empty() && boundary_file != "None")
        output_boundaries(read_bins, scales, reset_points,
                          bin_boundaries_folded,
                          hmm, distros, start_trans, trans, end_trans, classes,
                          boundary_file, boundary_score_file, VERBOSE);

      if (!read_counts_file.empty() && read_counts_file != "None") {
        write_read_counts_by_bin(bin_boundaries_folded,
                                 read_bins_a, read_bins_b, scales,
                                 classes, read_counts_file, VERBOSE);
      }
    }
    else if (mode == TEST_TEST_MODE) {
      /***********************************
       * STEP 3: ESTIMATE EMISSION PARAMS
       */
      const ThreeStateScaleSplitHMM hmm(min_prob, tolerance, max_iterations,
                                        VERBOSE);
      size_t state_num = 3;
      vector<SplitDistro> distros;
      vector<vector<double> > trans;
      vector<double> start_trans, end_trans;

      if (!in_param_file.empty())
        read_param_file(in_param_file, state_num,
                        start_trans, trans, end_trans, distros);
      else
        {
          if (VERBOSE)
            cout << "[ESTIMATING PARAMETERS]" << endl;

          fg_size = (fg_size > 0) ? fg_size : 6000;

          // All are using the fg name now
          training_size =
            (training_size == 0) ? read_bins.size() : training_size;

          vector<double> read_bins_sample, read_bins_a_sample,
            read_bins_b_sample, scales_sample;
          vector<size_t> reset_points_sample;
          pick_training_sample(
                               read_bins, read_bins_a, read_bins_b, scales, reset_points,
                               training_size, read_bins_sample, read_bins_a_sample,
                               read_bins_b_sample, scales_sample, reset_points_sample);

          distros.push_back(SplitDistro(fg_name));
          distros.push_back(SplitDistro(bg_name));
          distros.push_back(SplitDistro(fg_name));

          vector<double> mixing;
          ThreeStateScaleSplitResolveMixture(read_bins_sample, read_bins_a_sample,
                                             read_bins_b_sample, scales_sample,
                                             MAX_INITIALIZATION_ITR, tolerance, VERBOSE,
                                             distros.front(), distros[1], distros.back(), mixing);

          /***********************************
           * STEP 4: TRAIN THE HMM
           */
          set_transitions(bin_size, fg_size, mixing, VERBOSE,
                          start_trans, trans, end_trans);

          hmm.BaumWelchTraining(read_bins_sample, read_bins_a_sample,
                                read_bins_b_sample, scales_sample,
                                reset_points_sample,
                                start_trans, trans, end_trans,
                                distros.front(), distros[1], distros.back());

          clear_training_sample(read_bins_sample, read_bins_a_sample,
                                read_bins_b_sample, scales_sample, reset_points_sample);
        }

      if (!out_param_file.empty())
        write_param_file(out_param_file, state_num,
                         start_trans, trans, end_trans, distros);

      if (VERBOSE)
        report_final_values(distros, start_trans, trans, end_trans);

      /***********************************
       * STEP 5: DECODE THE DOMAINS
       */

      vector<size_t> classes;
      vector<double> scores;
      if (USE_POSTERIOR)
        hmm.PosteriorDecoding(read_bins, scales, reset_points,
                              start_trans, trans, end_trans,
                              distros.front(), distros[1], distros.back(),
                              classes, scores);
      else
        hmm.ViterbiDecoding(read_bins, scales, reset_points,
                            start_trans, trans, end_trans,
                            distros.front(), distros[1], distros.back(),
                            classes);

      /***********************************
       * STEP 6: WRITE THE RESULTS
       */

      // make sure the output dir is valid
      output_domains(read_bins, scales, reset_points, bin_boundaries_folded,
                     hmm, distros, start_trans, trans, end_trans,
                     classes, posterior_cutoff, undef_region_cutoff, cdf_cutoff,
                     domain_file, posterior_score_file, VERBOSE);
      if (!boundary_file.empty() && boundary_file != "None")
        output_boundaries(read_bins, scales, reset_points,
                          bin_boundaries_folded,
                          hmm, distros, start_trans, trans, end_trans, classes,
                          boundary_file, boundary_score_file, VERBOSE);


      if (!read_counts_file.empty() && read_counts_file != "None") {
        write_read_counts_by_bin(bin_boundaries_folded, read_bins_a, read_bins_b,
                                 scales, classes, read_counts_file, VERBOSE);
      }
    }
    else {
      cerr << "Please specifc the following value for mode: if you want to compare "
           << "a test sample and a control sample "
           << "and think there are two stats, use mode 2;"
           << "if you want to compare a test sample and another test sample "
           << "and think thre are three states, use mode 3" << std::endl;
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
