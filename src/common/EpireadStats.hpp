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

#ifndef EPIREAD_STATS
#define EPIREAD_STATS

#include "Epiread.hpp"
#include <vector>


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////
//////  FUNCTIONS FOR A SINGLE EPITYPE
//////

double
log_likelihood(const epiread &r, const std::vector<double> &a);
void
fit_epiallele(const std::vector<epiread> &reads, 
	      const std::vector<double> &indicators, std::vector<double> &a);
double
fit_single_epiallele(const std::vector<epiread> &reads, std::vector<double> &a);

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////
//////  FUNCTIONS FOR TWO EPITYPES
//////

double
log_likelihood(const epiread &r, const double z,
 	       const std::vector<double> &a1, const std::vector<double> &a2);
double
log_likelihood(const epiread &r, const std::vector<double> &a1, 
	       const std::vector<double> &a2);
double
log_likelihood(const std::vector<epiread> &reads, const std::vector<double> &indicators,
 	       const std::vector<double> &a1, const std::vector<double> &a2);

double
resolve_epialleles(const size_t max_itr,
		   const std::vector<epiread> &reads, 
		   std::vector<double> &indicators, 
		   std::vector<double> &a1, std::vector<double> &a2);

double
test_asm_lrt(const size_t max_itr, const double low_prob,
	     const double high_prob, std::vector<epiread> reads);

double
test_asm_lrt2(const size_t max_itr, const double low_prob,
	     const double high_prob, std::vector<epiread> reads);

double
test_asm_bic(const size_t max_itr, const double low_prob,
	     const double high_prob, std::vector<epiread> reads);


class EpireadStats {
public:
  EpireadStats(const double lp,
	       const double hp,
	       const double cv,
	       const size_t mi,
	       const bool UB) :
    low_prob(lp), high_prob(hp), 
    critical_value(cv), max_itr(mi),
    USE_BIC(UB) {}

  double 
  test_asm(const std::vector<epiread> &reads, bool &is_significant) const {
    const double score = (USE_BIC) ?
      test_asm_bic(max_itr, low_prob, high_prob, reads) :
      test_asm_lrt(max_itr, low_prob, high_prob, reads);
    is_significant = (score < critical_value || (USE_BIC && score < 0.0));
    return score;
  }
  
private:
  double low_prob;
  double high_prob;
  double critical_value;
  size_t max_itr;
  bool USE_BIC;
};

#endif
