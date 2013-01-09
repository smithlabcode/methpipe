/*    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith and Fang Fang
 *
 *    Authors: Fang Fang and Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "EpireadStats.hpp"

#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include <numeric>
#include<limits>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_cdf.h>

using std::string;
using std::vector;


static const double EPIREAD_STATS_TOLERANCE = 1e-20;

double
log_likelihood(const epiread &r, const double z, const vector<double> &a) {
  double ll = 0.0;
  for (size_t i = 0; i < r.seq.length(); ++i) {
    if (r.seq[i] == 'C' || r.seq[i] == 'T') {
      const double val = (r.seq[i] == 'C') ? a[r.pos + i] : (1.0 - a[r.pos + i]);
      assert(finite(val) && finite(log(val)));
      ll += z*log(val);
    }
  }
  return ll;
}

double
log_likelihood(const epiread &r, const double z,
 	       const vector<double> &a1, const vector<double> &a2) {
  return log_likelihood(r, z, a1) + log_likelihood(r, 1.0 - z, a2);
}

double
log_likelihood(const epiread &r, const vector<double> &a1, const vector<double> &a2) {
  return log(exp(log_likelihood(r, 1.0, a1)) + exp(log_likelihood(r, 1.0, a2)));
}

double
log_likelihood(const vector<epiread> &reads, const vector<double> &indicators,
	       const vector<double> &a1, const vector<double> &a2) {
  double ll = 0.0;
  for (size_t i = 0; i < reads.size(); ++i)
    ll += log_likelihood(reads[i], indicators[i], a1, a2);
  return ll;
}


static double
expectation_step(const vector<epiread> &reads, const double mixing,
		 const vector<double> &a1, const vector<double> &a2, 
		 vector<double> &indicators) {
  const double log_mixing1 = log(mixing);
  const double log_mixing2 = log(1.0 - mixing);
  assert(finite(log_mixing1) && finite(log_mixing2));
  
  double score = 0;
  for (size_t i = 0; i < reads.size(); ++i) {
    const double ll1 = log_mixing1 + log_likelihood(reads[i], 1.0, a1);
    const double ll2 = log_mixing2 + log_likelihood(reads[i], 1.0, a2);
    assert(finite(ll1) && finite(ll2));
    const double log_denom = log(exp(ll1) + exp(ll2));
    score += log_denom;
    indicators[i] = exp(ll1 - log_denom);
    assert(finite(log_denom) && finite(indicators[i]));
  }
  return score;
}


void
fit_epiallele(const vector<epiread> &reads, 
	      const vector<double> &indicators, vector<double> &a) {
  const size_t n_cpgs = a.size();
  vector<double> meth(n_cpgs, 0.0), total(n_cpgs, 0.0);
  for (size_t i = 0; i < reads.size(); ++i) {
    const size_t start = reads[i].pos;
    const double weight = indicators[i];
    for (size_t j = 0; j < reads[i].seq.length(); ++j)
      if (reads[i].seq[j] == 'C' || reads[i].seq[j] == 'T') {
	meth[start + j] += weight*(reads[i].seq[j] == 'C');
	total[start + j] += weight;
      }
  }
  for (size_t i = 0; i < n_cpgs; ++i)
    a[i] = (meth[i] + 0.5)/(total[i] + 1.0);
}


static void
maximization_step(const vector<epiread> &reads, const vector<double> &indicators,
		  double &mixing, vector<double> &a1, vector<double> &a2) {
  
  vector<double> inverted_indicators(indicators);
  for (size_t i = 0; i < inverted_indicators.size(); ++i)
    inverted_indicators[i] = 1.0 - inverted_indicators[i];
  
  // Fit the regular model parameters
  fit_epiallele(reads, indicators, a1);
  fit_epiallele(reads, inverted_indicators, a2);

  // Fit the mixing parameters
  //!!!! NO NEED BECAUSE THESE ARE ALWAYS 0.5!!!
  mixing = 0.5;
}


static double
expectation_maximization(const size_t max_itr, const vector<epiread> &reads, 
			 double &mixing, vector<double> &indicators, 
			 vector<double> &a1, vector<double> &a2) {
  //!!!! MIXING ALWAYS 0.5 IN CURRENT IMPLEMENTATION !!
  mixing = 0.5;
  double prev_score = -std::numeric_limits<double>::max();
  for (size_t i = 0; i < max_itr; ++i) {
    
    const double score = expectation_step(reads, mixing, a1, a2, indicators);
    maximization_step(reads, indicators, mixing, a1, a2);
    
    if ((prev_score - score)/prev_score < EPIREAD_STATS_TOLERANCE)
      break;
    prev_score = score;
  }
  return prev_score;
}


double
resolve_epialleles(const size_t max_itr, const vector<epiread> &reads, 
		   vector<double> &indicators, 
		   vector<double> &a1, vector<double> &a2) {
  // In current implementation mixing will never change
  /*static const*/ 
  double MIXING_PARAMETER = 0.5; 
  indicators.clear();
  indicators.resize(reads.size(), 0.0);
  for (size_t i = 0; i < reads.size(); ++i) {
    const double l1 = log_likelihood(reads[i], 1.0, a1);
    const double l2 = log_likelihood(reads[i], 1.0, a2);
    indicators[i] = exp(l1 - log(exp(l1) + exp(l2)));
  }
  
  return expectation_maximization(max_itr, reads, MIXING_PARAMETER, 
				  indicators, a1, a2);
}

double
fit_single_epiallele(const vector<epiread> &reads, vector<double> &a) {
  assert(reads.size() > 0);
  vector<double> indicators(reads.size(), 1.0);
  fit_epiallele(reads, indicators, a);
  
  double score = 0.0;
  for (size_t i = 0; i < reads.size(); ++i) {
    score += log_likelihood(reads[i], 1.0, a);
    assert(std::isfinite(score));
  }  
  return score;
}



double
fit_single_epiallele(const size_t read_start, const size_t read_end,
		     const size_t start, const size_t end,
		     const vector<epiread> &reads, vector<double> &a) {
  
  assert(reads.size() > 0);
  vector<double> indicators(read_end - read_start, 1.0);
  fit_epiallele(read_start, read_end, start, end, reads, indicators, a);
  
  double score = 0.0;
  for (size_t i = 0; i < reads.size(); ++i) {
    score += log_likelihood(start, end, reads[i], 1.0, a);
    assert(std::isfinite(score));
  }  
  return score;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////
//////   FUNCTIONS BELOW ARE FOR THE CASE WHEN NOT ALL PARTS OF SOME
//////   READS ARE TO BE INCLUDED, AND WHEN ONLY A SUBSET OF THE CPGS
//////   AVAILABLE WILL BE USED
//////
//////
//////
//////


double
log_likelihood(const size_t start, const size_t end,
	       const epiread &r, const double z,
	       const vector<double> &a) {
  double ll = 0.0;
  const size_t p = r.pos;
  string::const_iterator s(r.seq.begin());
  const size_t lim = r.seq.length();
  for (size_t i = 0; i < lim; ++i) {
    if (p + i >= start && p + i < end) {
      const char curr = *(s + i);
      if (curr == 'C' || curr == 'T') {
	const double val = (curr == 'C') ? 
	  a[p + i - start] : (1.0 - a[p + i - start]);
	ll += log(val);
	assert(finite(ll));
      }
    }
  }
  return z*ll;
}


static double
log_likelihood(const size_t start, const size_t end,
	       const epiread &r, const double z,
	       const vector<double> &a1, const vector<double> &a2) {
  return (log_likelihood(start, end, r, z, a1) + 
	  log_likelihood(start, end, r, 1.0 - z, a2));
}


double
log_likelihood(const size_t read_start, const size_t read_end,
	       const size_t start, const size_t end,
	       const vector<epiread> &reads, const vector<double> &indicators,
	       const vector<double> &a1, const vector<double> &a2) {
  double ll = 0.0;
  for (size_t i = read_start; i < read_end; ++i)
    ll += log_likelihood(start, end, reads[i], indicators[i - read_start], a1, a2);
  return ll;
}

static double
expectation_step(const size_t read_start, const size_t read_end,
		 const size_t start, const size_t end,
		 const vector<epiread> &reads, const double mixing,
		 const vector<double> &a1, const vector<double> &a2,
		 vector<double> &indicators) {

  const double log_mixing1 = log(mixing);
  const double log_mixing2 = log(1.0 - mixing);
  assert(finite(log_mixing1) && finite(log_mixing2));
  double score = 0;
  for (size_t i = read_start; i < read_end; ++i) {
    if (reads[i].pos + reads[i].seq.length() > start && reads[i].pos < end) {
      const double ll1 = log_mixing1 + log_likelihood(start, end, reads[i], 1.0, a1);
      const double ll2 = log_mixing2 + log_likelihood(start, end, reads[i], 1.0, a2);
      assert(finite(ll1) && finite(ll2));
      const double log_denom = log(exp(ll1) + exp(ll2));
      score += log_denom;
      indicators[i - read_start] = exp(ll1 - log_denom);
      assert(finite(log_denom) && finite(indicators[i - read_start]));
    }
  }
  return score;
}


void
fit_epiallele(const size_t read_start, const size_t read_end,
	      const size_t start, const size_t end,
	      const vector<epiread> &reads, const vector<double> &indicators,
	      vector<double> &a) {
  
  vector<double> meth(end - start, 0.0), total(end - start, 0.0);
  for (size_t i = read_start; i < read_end; ++i) {
    const size_t st = reads[i].pos;
    const double weight = indicators[i - read_start];
    const size_t lim = reads[i].seq.length();
    string::const_iterator s(reads[i].seq.begin());
    for (size_t j = 0; j < lim; ++j) {
      const size_t curr_pos = st + j;
      const char curr_base = *(s + j);
      if (curr_pos >= start && curr_pos < end && 
	  (curr_base == 'C' || curr_base == 'T')) {
	const size_t idx = curr_pos - start;
	meth[idx] += weight*(curr_base == 'C');
	total[idx] += weight;
      }
    }
  }
  const size_t lim = end - start;
  for (size_t i = 0; i < lim; ++i)
    a[i] = (meth[i] + 0.5)/(total[i] + 1.0);
}




static void
maximization_step(const size_t read_start, const size_t read_end,
		  const size_t start, const size_t end,
		  const vector<epiread> &reads, const vector<double> &indicators, 
		  const double mixing, vector<double> &a1, vector<double> &a2) {
  
  vector<double> inverted_indicators(indicators);
  for (size_t i = 0; i < inverted_indicators.size(); ++i)
    inverted_indicators[i] = 1.0 - inverted_indicators[i];
  
  // Fit the regular model parameters
  fit_epiallele(read_start, read_end, start, end, reads, indicators, a1);
  fit_epiallele(read_start, read_end, start, end, reads, inverted_indicators, a2);
  // Fit the mixing parameters
  // NO NEED BECAUSE THESE ARE ALWAYS 0.5!!!
  // mixing = MIXING_PARAMETER;
}


static double
expectation_maximization(const size_t read_start, const size_t read_end,
			 const size_t start, const size_t end,
			 const size_t max_itr, const vector<epiread> &reads, 
			 const double mixing, vector<double> &indicators, 
			 vector<double> &a1, vector<double> &a2) {
  double prev_score = -std::numeric_limits<double>::max();
  for (size_t i = 0; i < max_itr; ++i) {
    const double score = expectation_step(read_start, read_end, start, end, 
					  reads, mixing, a1, a2, indicators);
    maximization_step(read_start, read_end, start, end, 
		      reads, indicators, mixing, a1, a2);
    if ((prev_score - score)/prev_score < EPIREAD_STATS_TOLERANCE)
      break;
    prev_score = score;
  }
  return prev_score;
}


double
resolve_epialleles(const size_t read_start, const size_t read_end,
		   const size_t start, const size_t end,
		   const size_t max_itr, const vector<epiread> &reads,
		   vector<double> &indicators, 
		   vector<double> &a1, vector<double> &a2) {
  
  static const double MIXING_PARAMETER = 0.5;
  vector<double>(reads.size(), 0.0).swap(indicators);
  
  for (size_t i = read_start; i < read_end; ++i) {
    const double l1 = log_likelihood(start, end, reads[i], 1.0, a1);
    const double l2 = log_likelihood(start, end, reads[i], 1.0, a2);
    indicators[i - read_start] = exp(l1 - log(exp(l1) + exp(l2)));
  }
  return expectation_maximization(read_start, read_end, start, end, max_itr,
				  reads, MIXING_PARAMETER, indicators, a1, a2);
}



double
test_asm_lrt(const size_t max_itr, const double low_prob, const double high_prob, 
	     vector<epiread> reads) {
  
  adjust_read_offsets(reads);
  const size_t n_cpgs = get_n_cpgs(reads);
  
  // try a single epi-allele
  vector<double> a0(n_cpgs, 0.5);
  const double single_score = fit_single_epiallele(reads, a0);
  
  // initialize the pair epi-alleles and indicators, and do the actual
  // computation to infer alleles
  vector<double> a1(n_cpgs, low_prob), a2(n_cpgs, high_prob), indicators;
  resolve_epialleles(max_itr, reads, indicators, a1, a2);
  
  double log_likelihood_pair = 0.0;
  for (size_t i = 0; i < reads.size(); ++i)
    log_likelihood_pair += log_likelihood(reads[i], a1, a2);
  log_likelihood_pair += reads.size()*log(0.5);
  // static const double correction = 1.5;
  const size_t df = n_cpgs;
  
  const double llr_stat = -2*(single_score - log_likelihood_pair);
  const double p_value = 1.0 - gsl_cdf_chisq_P(llr_stat, df);
  return p_value;
}

double
test_asm_lrt2(const size_t max_itr, const double low_prob, const double high_prob, 
	      vector<epiread> reads) {
  
  adjust_read_offsets(reads);
  const size_t n_cpgs = get_n_cpgs(reads);
  
  // try a single epi-allele
  vector<double> a0(n_cpgs, 0.5);
  const double single_score = fit_single_epiallele(reads, a0);
  
  // initialize the pair epi-alleles and indicators, and do the actual
  // computation to infer alleles
  vector<double> a1(n_cpgs, low_prob), a2(n_cpgs, high_prob), indicators;
  resolve_epialleles(max_itr, reads, indicators, a1, a2);
  
  const double indic_sum = 
    std::accumulate(indicators.begin(), indicators.end(), 0.0);
  
  const double log_likelihood_pair = log_likelihood(reads, indicators, a1, a2)
    - gsl_sf_lnbeta(indic_sum, reads.size() - indic_sum)
    + reads.size()*log(0.5);
  const size_t df = n_cpgs + reads.size();
  
  const double llr_stat = -2*(single_score - log_likelihood_pair);
  const double p_value = 1.0 - gsl_cdf_chisq_P(llr_stat, df);
  return p_value;
}

double
test_asm_bic(const size_t max_itr, const double low_prob, const double high_prob,
	     vector<epiread> reads) {
  
  adjust_read_offsets(reads);
  const size_t n_cpgs = get_n_cpgs(reads);
  
  // try a single epi-allele
  vector<double> a0(n_cpgs, 0.5);
  const double single_score = fit_single_epiallele(reads, a0);
  
  // initialize the pair epi-alleles and indicators, and do the actual
  // computation to infer alleles
  vector<double> a1(n_cpgs, low_prob), a2(n_cpgs, high_prob), indicators;
  resolve_epialleles(max_itr, reads, indicators, a1, a2);
  
  const double indic_sum = 
    std::accumulate(indicators.begin(), indicators.end(), 0.0);
  
  const double pair_score = log_likelihood(reads, indicators, a1, a2)
    - gsl_sf_lnbeta(indic_sum, reads.size() - indic_sum)
    + gsl_sf_lnbeta(reads.size()/2.0, reads.size()/2.0)
    + reads.size()*log(0.5);
  
  const double bic_single = n_cpgs*log(reads.size()) - 2*single_score;
  const double bic_pair = 2*n_cpgs*log(reads.size()) - 2*pair_score;
  
  return bic_pair - bic_single;
}


