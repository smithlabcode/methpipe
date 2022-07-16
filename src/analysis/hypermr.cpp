/* Copyright (C) 2009 University of Southern California
 *                    Andrew D Smith
 * Author: Song Qiang, Andrew D. Smith
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
#include <stdexcept>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "zlib_wrapper.hpp"

#include "ThreeStateHMM.hpp"
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

using std::ostream_iterator;
using std::ofstream;


static GenomicRegion
as_gen_rgn(const MSite &s) {
  return GenomicRegion(s.chrom, s.pos, s.pos + 1);
}

static void
load_cpgs(const string &cpgs_file, vector<MSite> &cpgs,
          vector<pair<double, double> > &meth, vector<size_t> &reads) {

  igzfstream in(cpgs_file);
  if (!in)
    throw runtime_error("failed opening file: " + cpgs_file);

  MSite the_site;
  while (in >> the_site)
    cpgs.push_back(the_site);

  meth.resize(cpgs.size());
  reads.resize(cpgs.size());
  for (size_t i = 0; i < cpgs.size(); ++i) {
    meth[i] = make_pair(cpgs[i].n_meth(), cpgs[i].n_unmeth());
    reads[i] = cpgs[i].n_reads;
  }
}


template <class T, class U>
static void
separate_regions(const size_t desert_size,
                 vector<MSite> &cpgs,
                 vector<T> &meth, vector<U> &reads,
                 vector<size_t> &reset_points,
                 size_t &total_bases,
                 size_t &bases_in_deserts) {
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

  total_bases = 0;
  bases_in_deserts = 0;
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
    if (dist < numeric_limits<size_t>::max())
      total_bases += dist;

    prev_pos = cpgs[i].pos;
  }
  reset_points.push_back(cpgs.size());
}


static void
read_params_file(const string &params_file,
                 betabin &hypo_emission,
                 betabin &HYPER_emission,
                 betabin &HYPO_emission,
                 vector<vector<double> > &trans) {
  std::ifstream in(params_file);
  if (!in)
    throw runtime_error("failed to read param file: " + params_file);

  string hypo_emission_str;
  std::getline(in, hypo_emission_str);

  string HYPER_emission_str;
  std::getline(in, HYPER_emission_str);

  string HYPO_emission_str;
  std::getline(in, HYPO_emission_str);

  trans.resize(3, vector<double>(3, 0.0));
  for (size_t i = 0; i < trans.size(); ++i)
    for (size_t j = 0; j < trans[i].size(); ++j)
      in >> trans[i][j];
}

static void
write_params_file(const string &params_file,
                  const betabin &hypo_emission,
                  const betabin &HYPER_emission,
                  const betabin &HYPO_emission,
                  const vector<vector<double> > &trans) {
  ofstream out(params_file);
  out << hypo_emission.tostring() << endl
      << HYPER_emission.tostring() << endl
      << HYPO_emission.tostring() << endl;

  for (size_t i = 0; i < trans.size(); ++i) {
    for (size_t j = 0; j < trans[i].size(); ++j)
      out << trans[i][j] << "\t";
    out << endl;
  }
  out.close();
}


// static void
// assign_name(GenomicRegion &domain,
//             const STATE_LABELS &state, const size_t n_cpgs) {
//   domain.set_name((state == hypo ? "hypo:" : "hyper:") + to_string(n_cpgs));
// }


static void
build_domains(const bool VERBOSE,
              const vector<MSite> &cpgs,
              const vector<pair<double, double> > &meth,
              const vector<size_t> &reset_points,
              const vector<STATE_LABELS> &classes,
              vector<GenomicRegion> &domains) {
  domains.clear();

  for (size_t i = 0; i < reset_points.size() - 1; ++i)
    {
      const size_t start = reset_points[i];
      const size_t end = reset_points[i + 1];

      GenomicRegion domain(as_gen_rgn(cpgs[start]));
      STATE_LABELS prev_state = classes[start];
      size_t n = 1;
      // string hmrcpgs = cpgs[start].tostring();
      double meth_sum =
        meth[start].first / (meth[start].first + meth[start].second);

      for (size_t j = start + 1; j < end; ++j)
        {
          // if ((prev_state == hypo && classes[j] == HYPO)
          //     || (prev_state == HYPO && classes[j] == hypo))
          //     cerr << "WARNING: inconsist state sequences"
          //         " from posterior decoding" << endl;

          if ((prev_state == hypo && classes[j] == hypo)
              || (prev_state != hypo && classes[j] != hypo))
            {
              ++n;
              meth_sum += meth[j].first / (meth[j].first + meth[j].second);
              // hmrcpgs += ":" + cpgs[j].tostring();
            }
          else
            {
              domain.set_end(cpgs[j - 1].pos + 1);
              switch (prev_state) {
              case hypo:
                domain.set_name("hypo:" + smithlab::toa(n));
                break;
              case HYPER:
                domain.set_name("hyper:" + smithlab::toa(n));
                break;
              case HYPO:
                domain.set_name("hyper:" + smithlab::toa(n));
                break;
              }
              domain.set_score(meth_sum);
              domain.set_strand('+');
              if (prev_state == HYPER || prev_state == HYPO)
                domains.push_back(domain);

              domain = GenomicRegion(as_gen_rgn(cpgs[j]));
              n = 1;
              // hmrcpgs = cpgs[j].tostring();
              prev_state = classes[j];
              meth_sum = meth[j].first / (meth[j].first + meth[j].second);
            }
        }
      domain.set_end(cpgs[end - 1].pos + 1);
      switch (prev_state)
        {
        case hypo: domain.set_name("hypo:" + smithlab::toa(n)); break;
        case HYPER: domain.set_name("hyper:" + smithlab::toa(n)); break;
        case HYPO: domain.set_name("hyper:" + smithlab::toa(n)); break;
        }
      domain.set_score(meth_sum);
      domain.set_strand('+');
      if (prev_state == HYPER || prev_state == HYPO)
        domains.push_back(domain);
    }
}

// static void
// build_domains(const bool VERBOSE,
//               const vector<MSite> &cpgs,
//               const vector<size_t> &reset_points,
//               const vector<STATE_LABELS> &classes,
//               vector<GenomicRegion> &domains) {

//   STATE_LABELS prev_state = classes.front();
//   size_t n_cpgs = 0, reset_idx = 1, prev_pos = 0;
//   double total_meth = 0.0;
//   for (size_t i = 0; i < classes.size(); ++i) {
//     if (classes[i] != hypo)
//       throw runtime_error("found one!");
//     if (i == 0 || (reset_idx < reset_points.size() &&
//                    reset_points[reset_idx] == i) || classes[i] != prev_state) {
//       if (i > 0)  {
//         domains.back().set_end(prev_pos + 1);
//         assign_name(domains.back(), prev_state, n_cpgs);
//         domains.back().set_score(total_meth/n_cpgs);
//         n_cpgs = 0;
//         total_meth = 0.0;
//       }
//       domains.push_back(as_gen_rgn(cpgs[i]));
//       if (reset_points[reset_idx] == i)
//         ++reset_idx;
//     }
//     total_meth += cpgs[i].meth;
//     n_cpgs++;
//     prev_pos = cpgs[i].pos;
//     prev_state = classes[i];
//   }
//   domains.back().set_end(prev_pos + 1);
//   assign_name(domains.back(), prev_state, n_cpgs);
//   domains.back().set_score(total_meth/n_cpgs);
// }

static void
filter_domains(const bool VERBOSE, const double min_cumulative_meth,
               vector<GenomicRegion> &domains) {
  size_t j = 0;
  for (size_t i = 0; i < domains.size(); ++i)
    if (domains[i].get_score() >= min_cumulative_meth) {
      domains[j] = domains[i];
      ++j;
    }
  domains.erase(begin(domains) + j, end(domains));
}

static void
initialize_transitions(vector<vector<double> > &trans) {
  trans = vector<vector<double> >(3, vector<double>(3, 0.0));
  trans[hypo][hypo] = 0.99;
  trans[hypo][HYPER] = 1 - 0.99;
  trans[HYPER][hypo] = 1 - 0.95;
  trans[HYPER][HYPER] = 0.95 * 0.5;
  trans[HYPER][HYPO] = 0.95 * 0.5;
  trans[HYPO][HYPER] = 0.6;
  trans[HYPO][HYPO] = 0.3;
}

int
main_hypermr(int argc, const char **argv) {

  try {

    string outfile;
    string scores_file;
    string trans_file;

    size_t desert_size = 1000;
    size_t max_iterations = 10;

    // run mode flags
    bool VERBOSE = false;

    // corrections for small values (not parameters):
    double tolerance = 1e-10;

    double min_cumulative_meth = 4.0;
    bool USE_VITERBI_DECODING = false;

    string params_in_file;
    string params_out_file;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], "A program for segmenting DNA "
                           "methylation data"
                           "<cpg-BED-file>");
    opt_parse.add_opt("out", 'o', "output file (BED format)",
                      false, outfile);
    opt_parse.add_opt("scores", 's', "output file for posterior scores",
                      false, scores_file);
    opt_parse.add_opt("tolerance", 't', "Tolerance",
                      false, tolerance);
    opt_parse.add_opt("desert", 'd', "desert size", false, desert_size);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations);
    opt_parse.add_opt("viterbi", 'V', "Use Viterbi decoding",
                      false, USE_VITERBI_DECODING);
    opt_parse.add_opt("min-meth", 'M',
                      "min cumulative methylation level in HypeMR",
                      false, min_cumulative_meth);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("params-in", 'P', "HMM parameters file",
                      false, params_in_file);
    opt_parse.add_opt("params-out", 'p', "HMM parameters file",
                      false, params_out_file);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string cpgs_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "[loading_data]" << endl;
    vector<MSite> cpgs;
    vector<pair<double, double> > meth;
    vector<size_t> reads;
    load_cpgs(cpgs_file, cpgs, meth, reads);

    if (VERBOSE)
      cerr << "[n_sites=" << cpgs.size() << "]" << endl
           << "[mean_coverage="
           << accumulate(begin(reads), end(reads), 0.0)/reads.size()
           << "]" << endl;

    // separate the regions by chrom and by desert, and eliminate
    // those isolated CpGs
    vector<size_t> reset_points;
    size_t total_bases = 0, bases_in_deserts = 0;
    separate_regions(desert_size, cpgs, meth, reads, reset_points,
                     total_bases, bases_in_deserts);

    if (VERBOSE)
      cerr << "[n_sites_retained=" << cpgs.size() << "]" << endl
           << "[deserts_removed=" << reset_points.size() - 2 << "]" << endl
           << "[remaining_genome_fraction="
           << 1.0 - static_cast<double>(bases_in_deserts)/total_bases << "]"
           << endl;

    ThreeStateHMM hmm(meth, reset_points, tolerance, max_iterations, VERBOSE);

    vector<vector<double> > trans;
    initialize_transitions(trans);

    betabin hypo_emission, HYPER_emission, HYPO_emission;

    if (!params_in_file.empty())
      read_params_file(params_in_file,
                       hypo_emission, HYPER_emission, HYPO_emission, trans);
    else {
      const double n_reads =
        accumulate(begin(reads), end(reads), 0.0)/reads.size();
      const double fg_alpha = 0.33*n_reads;
      const double fg_beta = 0.67*n_reads;
      const double bg_alpha = 0.67*n_reads;
      const double bg_beta = 0.33*n_reads;

      hypo_emission = betabin(fg_alpha, fg_beta);
      HYPER_emission = betabin(bg_alpha, bg_beta);
      HYPO_emission = hypo_emission;
    }

    hmm.set_parameters(hypo_emission, HYPER_emission, HYPO_emission, trans);
    if (max_iterations > 0) hmm.BaumWelchTraining();
    hmm.get_parameters(hypo_emission, HYPER_emission, HYPO_emission, trans);

    if (!params_out_file.empty())
      write_params_file(params_out_file,
                        hypo_emission, HYPER_emission, HYPO_emission, trans);

    // DECODE THE STATES
    vector<STATE_LABELS> classes;
    if (USE_VITERBI_DECODING) hmm.ViterbiDecoding();
    else hmm.PosteriorDecoding();
    hmm.get_classes(classes);

    // DETERMINE THE DOMAINS
    vector<GenomicRegion> domains;
    build_domains(VERBOSE, cpgs, meth, reset_points, classes, domains);
    filter_domains(VERBOSE, min_cumulative_meth, domains);

    // WRITE THE RESULTS
    ofstream of;
    if (!outfile.empty()) of.open(outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    copy(begin(domains), end(domains), ostream_iterator<GenomicRegion>(out, "\n"));

    // IF REQUESTED, WRITE THE POSTERIOR SCORES
    if (!scores_file.empty()) {
      if (USE_VITERBI_DECODING)
        hmm.PosteriorDecoding();
      vector<Triplet> scores;
      hmm.get_state_posteriors(scores);
      ofstream score_out(scores_file.c_str());
      for (size_t i = 0; i < cpgs.size(); ++i) {
        score_out << cpgs[i] << "\t";
        switch (classes[i]) {
        case hypo:
          score_out << "hypo" << "\t" << scores[i].hypo << endl;
          break;
        case HYPER:
          score_out << "HYPER" << "\t" << scores[i].HYPER << endl;
          break;
        case HYPO:
          score_out << "HYPO" << "\t" << scores[i].HYPO << endl;
          break;
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
