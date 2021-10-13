/* Copyright (C) 2018 University of Southern California
 *                         Andrew D Smith
 * Authors: Andrew D. Smith, Song Qiang, Jenny Qu, Benjamin Decato
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

#include <stdexcept>
#include <numeric>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>

#include <unistd.h>
#include <gsl/gsl_sf.h>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "zlib_wrapper.hpp"
#include "bsutils.hpp"

#include "TwoStateHMM_PMD.hpp"
#include "MethpipeSite.hpp"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;
using std::max;
using std::min;
using std::make_pair;
using std::pair;
using std::runtime_error;

using std::ostream;
using std::ofstream;
using std::to_string;

template <class T> T &
methpipe_skip_header(T &in, string &line) {
  do {getline(in, line);}
  while (line[0] == '#');
  return in;
}

template <class T> T &
methpipe_read_site(const bool is_array_data,
                   T &in, string &chrom, size_t &pos,
                   string &strand, string &seq, double &meth,
                   size_t &coverage) {
  string line;
  methpipe_skip_header(in, line);

  /* GS: commenting out because it breaks when compiled with OSX clang
    std::istringstream iss;
    iss.rdbuf()->pubsetbuf(const_cast<char*>(line.c_str()), line.length()); */
  std::istringstream iss(line);
  iss >> chrom >> pos >> strand >> seq >> meth;
  if (!(iss >> chrom >> pos >> strand >> seq >> meth))
    throw std::runtime_error("bad methpipe site line: \"" + line + "\"");
  strand.resize(1);
  if (strand[0] != '-' && strand[0] != '+') {
    throw std::runtime_error("bad methpipe site line, strand incorrect: \"" +
                             line + "\"");
  }
  if (is_array_data)
    coverage = 1;
  else
    iss >> coverage;
  return in;
}


static void
get_adjacent_distances(const vector<GenomicRegion> &pmds,
                       vector<size_t> &dists) {
  for (size_t i = 1; i < pmds.size(); ++i)
    if (pmds[i].same_chrom(pmds[i - 1]))
      dists.push_back(pmds[i].get_start() - pmds[i - 1].get_end());
}


static bool
precedes(const string &chrom, const size_t position,
         const GenomicRegion &r) {
  return (chrom < r.get_chrom() ||
    (chrom == r.get_chrom() &&
     position < r.get_start()));
}


static bool
succeeds(const string &chrom, const size_t position,
         const GenomicRegion &r) {
  return (r.get_chrom() < chrom ||
    (chrom == r.get_chrom() &&
     r.get_end() <= position));
}


static void
merge_nearby_pmd(const size_t max_merge_dist,
     vector<GenomicRegion> &pmds) {
  size_t j = 0;
  for (size_t i = 1; i < pmds.size(); ++i) {
    if (pmds[j].same_chrom(pmds[i]) &&
        pmds[i].get_start() - pmds[j].get_end() <= max_merge_dist) {
      pmds[j].set_end(pmds[i].get_end());
      const string combined_name(pmds[j].get_name() + pmds[i].get_name());
      pmds[j].set_name(combined_name);
    }
    else pmds[++j] = pmds[i];
  }
  pmds.resize(j);
}


static size_t
find_best_bound(const bool IS_RIGHT_BOUNDARY,
                std::map<size_t, pair<size_t, size_t> > &pos_meth_tot,
                const vector<double> &fg_alpha,
                const vector<double> &fg_beta,
                const vector<double> &bg_alpha,
                const vector<double> &bg_beta) {

  vector<pair<size_t,size_t> > meth_tot;
  vector<size_t> positions;
  for (auto &&it : pos_meth_tot) {
    positions.push_back(it.first);
    meth_tot.push_back(it.second);
  }

  vector<size_t> cumu_left_meth(meth_tot.size(), 0);
  vector<size_t> cumu_left_tot(meth_tot.size(), 0);
  vector<size_t> cumu_right_meth(meth_tot.size(), 0);
  vector<size_t> cumu_right_tot(meth_tot.size(), 0);
  if (meth_tot.size() > 0)
    for (size_t i = 1; i < meth_tot.size()-1; ++i) {
      const size_t j = meth_tot.size() - 1 - i;
      cumu_left_meth[i] = cumu_left_meth[i-1] + meth_tot[i-1].first;
      cumu_left_tot[i] = cumu_left_tot[i-1] + meth_tot[i-1].second;
      cumu_right_meth[j] = cumu_right_meth[j+1] + meth_tot[j+1].first;
      cumu_right_tot[j] = cumu_right_tot[j+1] + meth_tot[j+1].second;
    }

  size_t best_idx = 0;
  double best_score = -std::numeric_limits<double>::max();
  if (meth_tot.size() > 0)
    for (size_t i = 1; i < meth_tot.size()-1; ++i) {
      size_t N_low, k_low, N_hi, k_hi;
      if (!IS_RIGHT_BOUNDARY) {
        N_low = cumu_right_tot[i] + meth_tot[i].second;
        k_low = cumu_right_meth[i] + meth_tot[i].first;
        N_hi = cumu_left_tot[i];
        k_hi = cumu_left_meth[i];
      }
      else {
        N_low = cumu_left_tot[i] + meth_tot[i].second;
        k_low = cumu_left_meth[i] + meth_tot[i].first;
        N_hi = cumu_right_tot[i];
        k_hi = cumu_right_meth[i];
      }
      if (N_hi > 0 && N_low > 0) {
        double score = 0;
        const double p_hi = static_cast<double>(k_hi)/N_hi;
        const double p_low = static_cast<double>(k_low)/N_low;

        for (size_t j = 0; j < fg_alpha.size(); ++j) {
          score += (((fg_alpha[j]-1.0)*log(p_low) +
                     ((fg_beta[j]-1.0)*log(1.0 - p_low)))
                    - gsl_sf_lnbeta(fg_alpha[j],fg_beta[j]))
            + (((bg_alpha[j]-1.0)*log(p_hi) +
                ((bg_beta[j]-1.0)*log(1.0 - p_hi)))
               - gsl_sf_lnbeta(bg_alpha[j],bg_beta[j]));
        } // beta max likelihood using learned emissions
        score /= fg_alpha.size();
        if (p_hi > p_low && score > best_score) {
          best_idx = i;
          best_score = score;
        }
      }
    }
  return (best_score > -std::numeric_limits<double>::max()) ?
    positions[best_idx] : std::numeric_limits<size_t>::max();
}


static void
get_boundary_positions(vector<GenomicRegion> &bounds,
           const vector<GenomicRegion> &pmds,
           const size_t &bin_size) {
  for (size_t i = 0; i < pmds.size(); ++i) {
    bounds.push_back(pmds[i]);
    pmds[i].get_start() > bin_size ?
      bounds.back().set_start(pmds[i].get_start() - bin_size) :
      bounds.back().set_start(0);
    bounds.back().set_end(pmds[i].get_start() + bin_size);

    bounds.push_back(pmds[i]);
    pmds[i].get_end() > bin_size ?
      bounds.back().set_start(pmds[i].get_end() - bin_size) :
      bounds.back().set_start(0);
    bounds.back().set_end(pmds[i].get_end() + bin_size);
  }
}


static void
get_optimized_boundary_likelihoods(const vector<string> &cpgs_file,
           vector<GenomicRegion> &bounds,
           const vector<bool> &array_status,
           const vector<double> &fg_alpha,
           const vector<double> &fg_beta,
           const vector<double> &bg_alpha,
           const vector<double> &bg_beta,
           vector<double> &boundary_scores,
           vector<size_t> &boundary_certainties) {
  // MAGIC NUMBER FOR WEIGHTING ARRAY
  // CONTRIBUTION TO BOUNDARY OBSERVATIONS
  static const double array_coverage_constant = 10;

  vector<igzfstream*> in(cpgs_file.size());
  for (size_t i = 0; i < cpgs_file.size(); ++i)
    in[i] = new igzfstream(cpgs_file[i]);

  std::map<size_t, pair<size_t, size_t> > pos_meth_tot;
  string chrom, site_name, strand;
  double meth_level;
  size_t position = 0ul, coverage = 0ul, n_meth = 0ul, n_unmeth = 0ul;
  size_t bound_idx = 0;
  for (; bound_idx < bounds.size(); ++bound_idx) { // for each boundary
    for (size_t i = 0; i < in.size(); ++i) {
      // get totals for all CpGs overlapping that boundary
      while (methpipe_read_site(array_status[i], *in[i], chrom, position,
                                strand, site_name, meth_level, coverage)
             && !succeeds(chrom, position, bounds[bound_idx])) {
        // check if CpG is inside boundary
        if (!precedes(chrom, position, bounds[bound_idx])) {
          if (array_status[i]) {
            if (meth_level != -1) {
              n_meth = round(meth_level*array_coverage_constant);
              n_unmeth = array_coverage_constant - n_meth;
            }
            else {
              n_meth = 0;
              n_unmeth = 0;
            }
          }
          else {
            n_meth = round(meth_level*coverage);
            n_unmeth = coverage - n_meth;
          }
          auto it(pos_meth_tot.find(position));
          if (it == end(pos_meth_tot)) {// does not exist in map
            pos_meth_tot[position] = make_pair(n_meth, n_meth+n_unmeth);
          }
          else { // add this file's contribution to the CpG's methylation
            pos_meth_tot[position].first += n_meth;
            pos_meth_tot[position].second += n_meth + n_unmeth;
          }
        }
      }
    }

    // Get the boundary position
    size_t boundary_position =
      (bounds[bound_idx].get_start() + bounds[bound_idx].get_end())/2;

    size_t N_low = 0, k_low = 0, N_hi = 0, k_hi = 0;
    for (auto it = begin(pos_meth_tot); it != end(pos_meth_tot); ++it) {
      if (it->first < boundary_position) {
        N_low += it->second.second;
        k_low += it->second.first;
      }
      else{
        N_hi += it->second.second;
        k_hi += it->second.first;
      }
    }

    double score = 0;
    const double p_hi = static_cast<double>(k_hi)/N_hi;
    const double p_low = static_cast<double>(k_low)/N_low;

    if (bound_idx %2 ) { // its a right boundary, p_low should go with fg
      for (size_t j = 0; j < fg_alpha.size(); ++j) {
        score += (((fg_alpha[j]-1.0)*log(p_low) +
                   ((fg_beta[j]-1.0)*log(1.0-p_low)))
                  - gsl_sf_lnbeta(fg_alpha[j],fg_beta[j]))
          + (((bg_alpha[j]-1.0)*log(p_hi) +
              ((bg_beta[j]-1.0)*log(1.0-p_hi)))
             - gsl_sf_lnbeta(bg_alpha[j],bg_beta[j]));
      }
    }
    else { // its a left boundary, p_low should go with bg
      for (size_t j = 0; j < fg_alpha.size(); ++j) {
        score +=  (((bg_alpha[j]-1.0)*log(p_low) +
                    ((bg_beta[j]-1.0)*log(1.0-p_low)))
                   - gsl_sf_lnbeta(bg_alpha[j],bg_beta[j]))
          + (((fg_alpha[j]-1.0)*log(p_hi) +
              ((fg_beta[j]-1.0)*log(1.0-p_hi)))
             - gsl_sf_lnbeta(fg_alpha[j],fg_beta[j]));
      }
    }
    boundary_certainties.push_back(std::min(N_low,N_hi));
    score /= fg_alpha.size();
    boundary_scores.push_back(exp(score));
    pos_meth_tot.clear();
  }

  for (size_t i = 0; i < cpgs_file.size(); ++i)
    delete in[i];
}


static void
find_exact_boundaries(const vector<string> &cpgs_file,
                      vector<GenomicRegion> &bounds,
                      const vector<bool> &array_status,
                      const vector<double> &fg_alpha,
                      const vector<double> &fg_beta,
                      const vector<double> &bg_alpha,
                      const vector<double> &bg_beta,
                      vector<size_t> &bound_site) {
  // MAGIC NUMBER FOR WEIGHTING ARRAY
  // CONTRIBUTION TO BOUNDARY OBSERVATIONS
  static const double array_coverage_constant = 10;

  vector<igzfstream*> in(cpgs_file.size());
  for (size_t i = 0; i < cpgs_file.size(); ++i)
    in[i] = new igzfstream(cpgs_file[i]);

  std::map<size_t, pair<size_t, size_t> > pos_meth_tot;
  string chrom, site_name, strand;
  double meth_level;
  size_t position = 0ul, coverage = 0ul, n_meth = 0ul, n_unmeth = 0ul;
  size_t bound_idx = 0;
  for (; bound_idx < bounds.size(); ++bound_idx) { // for each boundary
    for (size_t i = 0; i < in.size(); ++i) {
      // get totals for all CpGs overlapping that boundary
      while (methpipe_read_site(array_status[i], *in[i], chrom, position,
                                strand, site_name, meth_level, coverage)
             && !succeeds(chrom, position, bounds[bound_idx])) {
        // check if CpG is inside boundary
        if (!precedes(chrom, position, bounds[bound_idx])) {
          if (array_status[i]) {
            if (meth_level != -1) {
              n_meth = round(meth_level*array_coverage_constant);
              n_unmeth = array_coverage_constant - n_meth;
            }
            else {
              n_meth = 0;
              n_unmeth = 0;
            }
          }
          else {
            n_meth = round(meth_level*coverage);
            n_unmeth = coverage - n_meth;
          }
          auto it = pos_meth_tot.find(position);
          if (it == end(pos_meth_tot)) {// does not exist in map
            pos_meth_tot.emplace(position, make_pair(n_meth,n_meth+n_unmeth));
          }
          else { // add this file's contribution to the CpG's methylation
            pos_meth_tot[position].first += n_meth;
            pos_meth_tot[position].second += n_meth + n_unmeth;
          }
        }
      }
    }
    bound_site.push_back(find_best_bound(bound_idx % 2, pos_meth_tot,
                                         fg_alpha, fg_beta,
                                         bg_alpha, bg_beta));
    pos_meth_tot.clear();
  }
  for (size_t i = 0; i < in.size(); ++i)
    delete in[i];
}


static void
optimize_boundaries(const size_t bin_size,
                    const vector<string> &cpgs_file,
                    vector<GenomicRegion> &pmds,
                    const vector<bool> &array_status,
                    const vector<double> &fg_alpha,
                    const vector<double> &fg_beta,
                    const vector<double> &bg_alpha,
                    const vector<double> &bg_beta) {

  vector<GenomicRegion> bounds;
  get_boundary_positions(bounds, pmds, bin_size);
  vector<size_t> bound_site;
  find_exact_boundaries(cpgs_file, bounds, array_status, fg_alpha,
                        fg_beta, bg_alpha, bg_beta,
                        bound_site);

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ///// NOW RESET THE STARTS AND ENDS OF PMDs
  /////
  for (size_t i = 0; i < pmds.size(); ++i) {
    const size_t start_site = bound_site[2*i];
    if (start_site != std::numeric_limits<size_t>::max())
      pmds[i].set_start(start_site);
    const size_t end_site = bound_site[2*i + 1];
    if (end_site != std::numeric_limits<size_t>::max())
      pmds[i].set_end(end_site + 1);
  }

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ///// NOW MERGE THE PMDS THAT ARE TOO CLOSE
  /////
  vector<size_t> dists;
  get_adjacent_distances(pmds, dists);
  sort(begin(dists), end(dists));

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

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // LAST, GET THE CPG SITES WITHIN 1 BIN OF EACH BOUNDARY AND COMPUTE
  // THE LIKELIHOOD TO GET A "SCORE" ON THE BOUNDARY
  ///////
  bounds.clear(); // need updated boundaries after merging nearby PMDs
  get_boundary_positions(bounds, pmds, bin_size);

  vector<double> boundary_scores;
  vector<size_t> boundary_certainties;
  get_optimized_boundary_likelihoods(cpgs_file, bounds, array_status,
             fg_alpha, fg_beta, bg_alpha,
             bg_beta, boundary_scores, boundary_certainties);

  // Add the boundary scores to the PMD names
  for (size_t i = 0; i < pmds.size(); ++i)
    pmds[i].set_name(pmds[i].get_name()
                     + ":" + to_string(boundary_scores[2*i])
                     + ":" + to_string(boundary_certainties[2*i])
                     + ":" + to_string(boundary_scores[2*i+1])
                     + ":" + to_string(boundary_certainties[2*i+1]));
}


double
get_score_cutoff_for_fdr(const vector<double> &scores, const double fdr) {
  if (fdr <= 0)
    return numeric_limits<double>::max();
  else if (fdr > 1)
    return numeric_limits<double>::min();
  vector<double> local(scores);
  std::sort(begin(local), end(local));
  size_t i = 0;
  for (; i < local.size() - 1 &&
         local[i+1] < fdr*static_cast<double>(i+1)/local.size(); ++i);
  return local[i];
}


static inline double
score_contribution(const pair<double, double> &m) {
  const double denom = m.first + m.second;
  return (denom > 0) ? 1.0 - m.first/denom : 0.0;
}


static void
get_domain_scores(const vector<bool> &classes,
		  const vector<vector<pair<double, double> > > &meth,
		  const vector<size_t> &reset_points,
		  vector<double> &scores) {

  const size_t n_replicates = meth.size();
  size_t reset_idx = 1;
  bool in_domain = false;
  double score = 0;
  for (size_t i = 0; i < classes.size(); ++i) {
    if (reset_points[reset_idx] == i) {
      if (in_domain) {
        scores.push_back(score);
        score = 0;
        in_domain = false;
      }
      ++reset_idx;
    }
    if (classes[i]) {
      for (size_t r = 0; r < n_replicates ; ++r)
        score += score_contribution(meth[r][i]);
      in_domain = true;
    }
    else if (in_domain) {
      scores.push_back(score);
      score = 0;
      in_domain = false;
    }
  }
  if (in_domain)
    scores.push_back(score);
}


static void
build_domains(const bool VERBOSE,
              const vector<SimpleGenomicRegion> &cpgs,
              const vector<double> &post_scores,
              const vector<size_t> &reset_points,
              const vector<bool> &classes,
              const vector<size_t> &dists_btwn_cpgs,
              vector<GenomicRegion> &domains) {

  size_t n_cpgs = 0, reset_idx = 1, prev_end = 0;
  bool in_domain = false;
  double score = 0;
  for (size_t i = 0; i < classes.size(); ++i) {
    if (reset_points[reset_idx] == i) {
      if (in_domain) {
        domains.back().set_end(prev_end);
        domains.back().set_score(n_cpgs);
        n_cpgs = 0;
        score = 0;
        in_domain = false;
      }
      ++reset_idx;
    }
    if (classes[i]) {
      if (!in_domain) {
        domains.push_back(GenomicRegion(cpgs[i]));
        in_domain = true;
      }
      ++n_cpgs;
      score += post_scores[i];
    }
    else if (in_domain) {
      domains.back().set_end(prev_end);
      domains.back().set_score(n_cpgs);
      n_cpgs = 0;
      score = 0;
      in_domain = false;
    }
    prev_end = cpgs[i].get_end();
  }
}


//Modified to take multiple replicates
template <class T, class U> static void
separate_regions(const bool VERBOSE, const size_t desert_size,
                 vector<vector<SimpleGenomicRegion> > &cpgs,
                 vector<vector<T> > &meth, vector<vector<U> > &reads,
                 vector<size_t> &reset_points,
                 vector<size_t> &dists_btwn_cpgs) {
  if (VERBOSE)
    cerr << "[SEPARATING BY CPG DESERT]" << endl;
  const size_t n_replicates = cpgs.size();

  // eliminate the zero-read cpg sites if no coverage in any rep
  size_t end_coord_of_prev = 0;
  size_t j = 0;
  for (size_t i = 0; i < cpgs[0].size(); ++i) {
    bool all_empty = true;
    size_t rep_idx = 0;
    while (all_empty && rep_idx < n_replicates) {
      if (reads[rep_idx][i] == 0) ++rep_idx;
      else all_empty = false;
    }
    if (!all_empty) {
      dists_btwn_cpgs.push_back(cpgs[0][i].get_start() - end_coord_of_prev);
      end_coord_of_prev = cpgs[0][i].get_end();
      for (size_t r = 0; r < n_replicates; ++r) {
        cpgs[r][j] = cpgs[r][i];
        meth[r][j] = meth[r][i];
        reads[r][j] = reads[r][i];
      }
      ++j;
    }
  }

  for (size_t r = 0; r < n_replicates; ++r) {
    cpgs[r].resize(j);
    meth[r].resize(j);
    reads[r].resize(j);
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
         << "NUMBER OF DISTANCES BETWEEN: " << dists_btwn_cpgs.size() << endl
         << "DESERTS REMOVED: " << reset_points.size() - 2 << endl << endl;
}


static void
shuffle_cpgs(const TwoStateHMM &hmm,
	     vector<vector<pair<double, double> > > meth,
	     const vector<size_t> &reset_points,
	     const vector<double> &start_trans,
	     const vector<vector<double> > &trans,
	     const vector<double> &end_trans,
	     const vector<double> &fg_alpha, const vector<double> &fg_beta,
	     const vector<double> &bg_alpha, const vector<double> &bg_beta,
	     vector<double> &domain_scores,
	     vector<bool> &array_status) {

  size_t n_replicates = meth.size();

  for (size_t r =0 ; r < n_replicates; ++r)
    random_shuffle(begin(meth[r]), end(meth[r]));
  vector<bool> classes;
  vector<double> scores;
  hmm.PosteriorDecoding_rep(meth, reset_points, start_trans, trans,
                            end_trans, fg_alpha, fg_beta, bg_alpha,
                            bg_beta, classes, scores, array_status);
  get_domain_scores(classes, meth, reset_points, domain_scores);
  sort(begin(domain_scores), end(domain_scores));
}

static void
assign_p_values(const vector<double> &random_scores,
                const vector<double> &observed_scores,
                vector<double> &p_values) {
  const double n_randoms =
    random_scores.size() == 0 ? 1 : random_scores.size();
  for (size_t i = 0; i < observed_scores.size(); ++i)
    p_values.push_back((end(random_scores) -
                        upper_bound(begin(random_scores),
                                    end(random_scores),
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
write_posteriors_file(const string &posteriors_file,
                      const vector<vector<SimpleGenomicRegion> > &cpgs,
                      const vector<double> &scores) {
  static const size_t decimal_precision = 10;

  ofstream out(posteriors_file);
  out.precision(decimal_precision);
  for (size_t r = 0; r < scores.size(); ++r)
    out << cpgs[0][r] << '\t' << scores[r] << endl;
}


static void
write_params_file(const string &outfile,
                  const vector<double> &fg_alpha,
                  const vector<double> &fg_beta,
                  const vector<double> &bg_alpha,
                  const vector<double> &bg_beta,
                  const vector<double> &start_trans,
                  const vector<vector<double> > &trans,
                  const vector<double> &end_trans) {
  static const size_t decimal_precision = 30;
  ofstream out(outfile);
  out.precision(decimal_precision);
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


static bool
check_if_array_data(const string &infile) {
  igzfstream in(infile);
  if (!in)
    throw std::runtime_error("bad file: " + infile);

  std::string line;
  getline(in, line);
  std::istringstream iss(line);
  std::string chrom, pos, strand, seq, meth, cov;
  iss >> chrom >> pos >> strand >> seq >> meth;
  return (!(iss >> cov));
}


static void
load_array_data(const size_t bin_size,
                const string &cpgs_file,
                vector<SimpleGenomicRegion> &intervals,
                vector<pair<double, double> > &meth,
                vector<size_t> &reads) {
  // MAGIC. GS: minimum value for array?
  static const double meth_min = 1.0e-2;

  igzfstream in(cpgs_file);
  string curr_chrom;
  size_t prev_pos = 0ul, curr_pos = 0ul;
  double array_meth_bin = 0.0;
  double num_probes_in_bin = 0.0;
  string chrom, site_name, strand;
  double meth_level;
  size_t position = 0ul, coverage = 0ul;

  while (methpipe_read_site(true, in, chrom, position, strand, site_name,
                            meth_level, coverage)) {
    if (meth_level != -1 ) { // its covered by a probe
      ++num_probes_in_bin;
      if (meth_level < meth_min)
        array_meth_bin += meth_min;
      else if (meth_level > 1.0 - meth_min)
        array_meth_bin += (1.0 - meth_min);
      else
        array_meth_bin += meth_level;
    }

    if (curr_chrom != chrom) {
      if (!curr_chrom.empty()) {
        if (chrom < curr_chrom)
          throw runtime_error("CpGs not sorted in file \""
                              + cpgs_file + "\"");
        intervals.push_back(SimpleGenomicRegion(curr_chrom, curr_pos,
                                                prev_pos + 1));
        meth.push_back(make_pair(array_meth_bin, num_probes_in_bin));
        (num_probes_in_bin > 0) ? reads.push_back(1) : reads.push_back(0);
      }
      curr_chrom = chrom;
      curr_pos = position;
      array_meth_bin = 0.0;
      num_probes_in_bin = 0.0;
    }
    else if (curr_pos > position) {
      cerr << curr_pos << "\t" << position << endl;
      throw std::runtime_error("CpGs not sorted in file \"" + cpgs_file + "\"");
    }
    else if (position > curr_pos + bin_size) {
      intervals.push_back(SimpleGenomicRegion(curr_chrom, curr_pos,
                                              curr_pos + bin_size));
      meth.push_back(make_pair(array_meth_bin, num_probes_in_bin));
      (num_probes_in_bin > 0) ? reads.push_back(1) : reads.push_back(0);

      array_meth_bin = 0.0;
      num_probes_in_bin = 0.0;
      curr_pos += bin_size;
      while (curr_pos + bin_size < position) {
        intervals.push_back(SimpleGenomicRegion(curr_chrom, curr_pos,
                                                curr_pos + bin_size));
        reads.push_back(0);
        meth.push_back(make_pair(0.0, 0.0));
        curr_pos += bin_size;
      }
    }
  }

  if (meth_level != -1 ) { // its covered by a probe
    ++num_probes_in_bin;
    if (meth_level < meth_min)
      array_meth_bin += meth_min;
    else if (meth_level > 1.0 - meth_min)
      array_meth_bin += (1.0 - meth_min);
    else
      array_meth_bin += meth_level;
  }
  prev_pos = position;
  if (!curr_chrom.empty()) {
    intervals.push_back(SimpleGenomicRegion(curr_chrom, curr_pos,
                                            prev_pos + 1));
    meth.push_back(make_pair(array_meth_bin, num_probes_in_bin));
    if (num_probes_in_bin > 0)
      reads.push_back(1);
    else reads.push_back(0);
  }
}


static void
load_wgbs_data(const size_t bin_size,
               const string &cpgs_file,
               vector<SimpleGenomicRegion> &intervals,
               vector<pair<double, double> > &meth,
               vector<size_t> &reads) {

  igzfstream in(cpgs_file);

  string curr_chrom;
  size_t prev_pos = 0ul, curr_pos = 0ul;
  size_t n_meth_bin = 0ul, n_unmeth_bin = 0ul;
  string chrom, site_name, strand;
  double meth_level;
  size_t position = 0ul, coverage = 0ul, n_meth = 0ul, n_unmeth = 0ul;

  while (methpipe_read_site(false, in, chrom, position,
                            strand, site_name, meth_level, coverage)) {

    n_meth = round(meth_level * coverage);
    n_unmeth = coverage - n_meth;

    // no range, or no chrom
    if (curr_chrom != chrom) {
      if (!curr_chrom.empty()) {
        if (chrom < curr_chrom)
          throw std::runtime_error("CpGs not sorted in file \""
                                   + cpgs_file + "\"");
        intervals.push_back(SimpleGenomicRegion(curr_chrom, curr_pos,
                                                prev_pos + 1));
        reads.push_back(n_meth_bin + n_unmeth_bin);
        meth.push_back(make_pair(n_meth_bin, n_unmeth_bin));
      }
      curr_chrom = chrom;
      curr_pos = position;
      n_meth_bin = 0ul;
      n_unmeth_bin = 0ul;
    }
    else if (curr_pos > position) {
      throw std::runtime_error("CpGs not sorted in file \"" + cpgs_file + "\"");
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
    intervals.push_back(SimpleGenomicRegion(curr_chrom,
                                            curr_pos, prev_pos + 1));
    reads.push_back(n_meth_bin + n_unmeth_bin);
    meth.push_back(make_pair(n_meth_bin, n_unmeth_bin));
  }
}

static void
apply_ci_criterion(const double conf_level,
		   const size_t &n_reads,
		   size_t &min_cov_to_pass,
		   size_t &total_passed) {
  static const double fixed_phat = 0.5; // magic
  double lower = 0.0, upper = 0.0;
  wilson_ci_for_binomial(1.0 - conf_level, n_reads,
			 fixed_phat, lower, upper);
  if ((upper - lower) < (1.0 - conf_level)) {
    if (n_reads < min_cov_to_pass)
      min_cov_to_pass = n_reads;
    ++total_passed;
  }
}


/* Bin first chromosome into bins of increasing size, starting
 * at @param bin_size, until a sufficient percentage of bins
 * contain confident binomial estimates of the true methylation
 * level in the bin.
 */
static double
binsize_selection(const bool &VERBOSE,
                  const size_t initial_bin_size, const string &cpgs_file,
                  const double conf_level, const double ACCEPT_THRESHOLD) {
  static const size_t bin_size_step = 500; // MAGIC (ADS: must be fixed)
  static const size_t max_bin_size = 50000; // MORE MAGIC...

  size_t bin_size = initial_bin_size;

  for (; bin_size < max_bin_size; bin_size += bin_size_step) {

    if (VERBOSE)
      cerr << "evaluating bin size " << bin_size << endl;

    // ADS: loading data each iteration should be put outside loop
    igzfstream in(cpgs_file);
    size_t curr_pos = 0ul;
    size_t n_reads_bin = 0ul;
    string chrom, site_name, strand;
    double meth_level;
    size_t position = 0ul, coverage = 0;

    vector<size_t> reads;
    string first_chrom;

    while (methpipe_read_site(false, in, chrom, position, strand, site_name,
            meth_level, coverage) &&
     (first_chrom.empty() || first_chrom == chrom)) {
      if (curr_pos > position)
        throw std::runtime_error("CpGs not sorted in file: " + cpgs_file);
      else if (position > curr_pos + bin_size) {
        reads.push_back(n_reads_bin);
        n_reads_bin = 0ul;
        curr_pos += bin_size;
        while (curr_pos + bin_size < position) {
          reads.push_back(0);
          curr_pos += bin_size;
        }
      }
      n_reads_bin += coverage;
      swap(chrom, first_chrom);
    }

    size_t total_passed = 0, total_covered = 0;
    size_t total = 0;

    // min_cov_to_pass will change during iteration
    size_t min_cov_to_pass = numeric_limits<size_t>::max();

    double prev_frac_passed = 0.0;
    for (size_t i = 0; i < reads.size(); ++i) {
      if (reads[i] > 0) {
        ++total_covered;
	apply_ci_criterion(conf_level, reads[i],
			   min_cov_to_pass, total_passed);
      }
      ++total;
    }

    // fails if minimum coverage could not be calculated
    if (min_cov_to_pass == numeric_limits<size_t>::max()) {
      throw runtime_error("Insufficient coverage in chrom " +
                          first_chrom + " to estimate a bin size");
    }
    if (VERBOSE)
      cerr << "Min. cov. to pass: " << min_cov_to_pass << endl;

    const double frac_passed = static_cast<double>(total_passed)/total;
    if (frac_passed > ACCEPT_THRESHOLD)
      return bin_size;
    else {
      if (VERBOSE)
        cerr << "% bins passed: " << (double)total_passed/total << endl;

      if (frac_passed < prev_frac_passed)
        throw runtime_error("Insufficient coverage in chrom " +
                            first_chrom + " to estimate a bin size");

      prev_frac_passed = frac_passed;
    }
  }
  return bin_size;
}


static void
load_intervals(const size_t bin_size,
               const string &cpgs_file,
               vector<SimpleGenomicRegion> &intervals,
               vector<pair<double, double> > &meth,
               vector<size_t> &reads, vector<bool> &array_status) {

  const bool is_array_data = check_if_array_data(cpgs_file);

  array_status.push_back(is_array_data);

  if (is_array_data)
    load_array_data(bin_size, cpgs_file, intervals, meth, reads);
  else
    load_wgbs_data(bin_size, cpgs_file, intervals, meth, reads);
}


int
main(int argc, const char **argv) {
  try {

    const char* sep = ",";
    string outfile;

    size_t rng_seed = 408;

    bool DEBUG = false;
    size_t desert_size = 5000;
    size_t bin_size = 1000;
    size_t max_iterations = 10;
    // run mode flags
    bool VERBOSE = false;
    bool ARRAY_MODE = false;
    bool fixed_bin_size = false;

    // MAGIC: corrections for small values (not parameters):
    static const double tolerance = 1e-5;
    static const double min_prob  = 1e-10;

    string params_in_files;
    string params_out_file;
    string posteriors_out_prefix;


    const string description =
      "Identify PMDs in methylomes. Methylation must be provided in the \
      methcounts file format (chrom, position, strand, context, \
      methylation, reads). See the methcounts documentation for \
      details. This program assumes only data at CpG sites and that \
      strands are collapsed so only the positive site appears in the \
      file, but reads counts are from both strands.";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<methcount-files>");
    opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("desert", 'd', "max dist between bins with data in PMD",
                      false, desert_size);
    opt_parse.add_opt("fixedbin", 'f', "Fixed bin size", false, fixed_bin_size);
    opt_parse.add_opt("bin", 'b', "Starting bin size", false, bin_size);
    opt_parse.add_opt("arraymode",'a', "All samples are array", false, ARRAY_MODE);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("debug", 'D', "print more run info", false, DEBUG);
    opt_parse.add_opt("params-in", 'P', "HMM parameter files for "
                      "individual methylomes (separated with comma)",
                      false, params_in_files);
    opt_parse.add_opt("posteriors-out", 'r',
                      "write out posterior probabilities in methcounts format",
                      false, posteriors_out_prefix);
    opt_parse.add_opt("params-out", 'p', "write HMM parameters to this file",
                      false, params_out_file);
    opt_parse.add_opt("seed", 's', "specify random seed", false, rng_seed);

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
    /****************** END COMMAND LINE OPTIONS *****************/

    srand(rng_seed);

    vector<string> cpgs_file = leftover_args;
    vector<string> params_in_file;
    if (!params_in_files.empty()) {
      params_in_file = smithlab::split(params_in_files, sep, false);
      assert(cpgs_file.size() == params_in_file.size());
    }

    size_t n_replicates = cpgs_file.size();

    // Sanity checks input file format and dynamically selects bin
    // size from WGBS samples.
    if (!fixed_bin_size && !ARRAY_MODE) {
      if (VERBOSE)
        cerr << "[DYNAMICALLY SELECTING BIN SIZE]" << endl;
      double confidence_interval = 0.80;
      double prop_accept = 0.80;
      for (size_t i = 0; i < n_replicates; ++i) {
        const bool arrayData = check_if_array_data(cpgs_file[i]);
        if (!arrayData) {
          bin_size = binsize_selection(VERBOSE, bin_size, cpgs_file[i],
                                       confidence_interval, prop_accept);
          desert_size = 5*bin_size; // TODO: explore extrapolation number
        }
        else {
          // same as the parameters below
          bin_size = 1000;
          desert_size = 200000;
        }
      }
    }

    if (ARRAY_MODE) {
      bin_size = 1000;    // MAGIC NUMBERS FROM PAPER
      desert_size = 200000; // PERFORM WITH HIGHEST JACCARD INDEX TO WGBS
    }
    if (VERBOSE)
      cerr << "[READING IN AT BIN SIZE " << bin_size << "]" << endl;

    // separate the regions by chrom and by desert
    vector<vector<SimpleGenomicRegion> > cpgs(n_replicates);
    vector<vector<pair<double, double> > > meth(n_replicates);
    vector<vector<size_t> > reads(n_replicates);
    vector<bool> array_status;

    for (size_t i = 0; i < n_replicates; ++i) {
      if (VERBOSE)
        cerr << "[READING CPGS AND METH PROPS] from " << cpgs_file[i] << endl;

      load_intervals(bin_size, cpgs_file[i], cpgs[i], meth[i],
                     reads[i], array_status);
      if (VERBOSE)
        cerr << "TOTAL CPGS: " << cpgs[i].size() << endl
             << "MEAN COVERAGE: "
             << accumulate(begin(reads[i]), end(reads[i]), 0.0)/reads[i].size()
             << endl;
    }

    // ensure replicates have same number of sites
    const size_t n_sites = cpgs[0].size();
    for (size_t i = 0; i < n_replicates; ++i)
      if (cpgs[i].size() != n_sites)
        throw runtime_error("inputs contain different number of sites");

    // separate regions by chrom and desert; eliminate isolated CpGs
    vector<size_t> reset_points;
    vector<size_t> dists_btwn_cpgs;
    separate_regions(VERBOSE, desert_size, cpgs, meth, reads,
                     reset_points, dists_btwn_cpgs);

    /****************** Read in params *****************/
    vector<double> start_trans(2, 0.5), end_trans(2, 1e-10);
    vector<vector<double> > trans(2, vector<double>(2, 0.01));
    trans[0][0] = trans[1][1] = 0.99;
    const TwoStateHMM hmm(min_prob, tolerance, max_iterations, VERBOSE, DEBUG);
    vector<double> reps_fg_alpha(n_replicates, 0.05);
    vector<double> reps_fg_beta(n_replicates, 0.95);
    vector<double> reps_bg_alpha(n_replicates, 0.95);
    vector<double> reps_bg_beta(n_replicates, 0.05);
    double score_cutoff_for_fdr = std::numeric_limits<double>::max();

    if (!params_in_file.empty()) {
      // read parameters files
      for (size_t i= 0; i< n_replicates; ++i)
        read_params_file(VERBOSE, params_in_file[i], reps_fg_alpha[i],
                         reps_fg_beta[i], reps_bg_alpha[i], reps_bg_beta[i],
                         start_trans, trans, end_trans, score_cutoff_for_fdr);
    }

    // train model (default behavior; not done when params supplied)
    if (max_iterations > 0)
      hmm.BaumWelchTraining_rep(meth, reset_points,
                                start_trans, trans, end_trans,
                                reps_fg_alpha, reps_fg_beta,
                                reps_bg_alpha, reps_bg_beta, array_status);

    if (!params_out_file.empty()) {
      // write all the HMM parameters
      write_params_file(params_out_file,
                        reps_fg_alpha, reps_fg_beta,
                        reps_bg_alpha, reps_bg_beta,
                        start_trans, trans, end_trans);
    }

    /***********************************/

    if (!posteriors_out_prefix.empty()) {
      vector<double> into_scores;
      hmm.TransitionPosteriors_rep(meth, reset_points, start_trans, trans,
                                   end_trans, reps_fg_alpha, reps_fg_beta,
                                   reps_bg_alpha, reps_bg_beta, array_status,
                                   2, into_scores);
      write_posteriors_file(posteriors_out_prefix + ".intoTrans", cpgs,
                            into_scores);
      vector<double> outof_scores;
      hmm.TransitionPosteriors_rep(meth, reset_points, start_trans, trans,
                                   end_trans, reps_fg_alpha, reps_fg_beta,
                                   reps_bg_alpha, reps_bg_beta, array_status,
                                   1, outof_scores);
      write_posteriors_file(posteriors_out_prefix + ".outofTrans", cpgs,
                            outof_scores);
    }

    vector<bool> classes;
    vector<double> scores;
    hmm.PosteriorDecoding_rep(meth, reset_points, start_trans, trans, end_trans,
                              reps_fg_alpha, reps_fg_beta,
                              reps_bg_alpha, reps_bg_beta, classes,
                              scores, array_status);

    if (!posteriors_out_prefix.empty()) {
      write_posteriors_file(posteriors_out_prefix + ".emissions", cpgs, scores);
    }

    vector<double> domain_scores;
    get_domain_scores(classes, meth, reset_points, domain_scores);

    if (VERBOSE)
      cerr << "[RANDOMIZING SCORES FOR FDR]" << endl;

    vector<double> random_scores;
    shuffle_cpgs(hmm, meth, reset_points, start_trans, trans, end_trans,
		 reps_fg_alpha, reps_fg_beta, reps_bg_alpha, reps_bg_beta,
		 random_scores, array_status);

    vector<double> p_values;
    assign_p_values(random_scores, domain_scores, p_values);

    if (score_cutoff_for_fdr == numeric_limits<double>::max() &&
        !p_values.empty())
      score_cutoff_for_fdr = get_score_cutoff_for_fdr(p_values, 0.01);

    if (!params_out_file.empty()) {
      ofstream out(params_out_file, std::ios::app);
      out << "FDR_CUTOFF\t"
          << std::setprecision(30) << score_cutoff_for_fdr << endl;
    }

    vector<GenomicRegion> domains;
    build_domains(VERBOSE, cpgs[0], scores, reset_points, classes,
                  dists_btwn_cpgs, domains);

    size_t good_pmd_count = 0;
    vector<GenomicRegion> good_domains;
    for (size_t i = 0; i < domains.size(); ++i)
      if (p_values[i] < score_cutoff_for_fdr) {
        good_domains.push_back(domains[i]);
        good_domains.back().set_name("PMD" + to_string(good_pmd_count++));
      }

    optimize_boundaries(bin_size, cpgs_file, good_domains, array_status,
      reps_fg_alpha, reps_fg_beta, reps_bg_alpha, reps_bg_beta);

    ofstream of;
    if (!outfile.empty()) of.open(outfile);
    ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    copy(begin(good_domains), end(good_domains),
         std::ostream_iterator<GenomicRegion>(out, "\n"));
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
