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
log_likelihood(const epiread &r, const double z, const std::vector<double> &a);
double
log_likelihood(const size_t start, const size_t end, 
	       const epiread &r, const double z, const std::vector<double> &a);
void
fit_epiallele(const std::vector<epiread> &reads, 
	      const std::vector<double> &indicators, std::vector<double> &a);
void
fit_epiallele(const size_t read_start, const size_t read_end,
	      const size_t start, const size_t end,
	      const std::vector<epiread> &reads, 
	      const std::vector<double> &indicators, std::vector<double> &a);
double
fit_single_epiallele(const std::vector<epiread> &reads, std::vector<double> &a);
double
fit_single_epiallele(const size_t read_start, const size_t read_end,
		     const size_t start, const size_t end,
		     const std::vector<epiread> &reads, std::vector<double> &a);

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
log_likelihood(const size_t read_start, const size_t read_end,
	       const size_t start, const size_t end,
	       const std::vector<epiread> &reads, const std::vector<double> &indicators,
	       const std::vector<double> &a1, const std::vector<double> &a2);

double
resolve_epialleles(const size_t max_itr,
		   const std::vector<epiread> &reads, 
		   std::vector<double> &indicators, 
		   std::vector<double> &a1, std::vector<double> &a2);
double
resolve_epialleles(const size_t read_start, const size_t read_end,
		   const size_t start, const size_t end,
		   const size_t max_itr, const std::vector<epiread> &reads,
		   std::vector<double> &indicators, 
		   std::vector<double> &a1, std::vector<double> &a2);

double
test_asm_lrt(const size_t max_itr, const double low_prob, const double high_prob, 
	     std::vector<epiread> reads);

double
test_asm_lrt2(const size_t max_itr, const double low_prob, const double high_prob, 
	      std::vector<epiread> reads);

double
test_asm_bic(const size_t max_itr, const double low_prob, const double high_prob,
	     std::vector<epiread> reads);


#endif
