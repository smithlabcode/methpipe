/*
 *    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith, Song Qiang
 *
 *    Authors: Andrew D. Smith, Song Qiang
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

#include <cmath>

#include <gsl/gsl_sf.h>

#include "contingency-table.hpp"
#include "numerical_utils.hpp"

using std::min;

static inline double
log_prob_hypergeo(const size_t meth_a, const size_t unmeth_a, 
                  const size_t meth_b, const size_t unmeth_b,
                  const size_t k)
{
    return  gsl_sf_lnchoose(meth_b + unmeth_b - 1, k) + 
        gsl_sf_lnchoose(meth_a + unmeth_a - 1, meth_a + meth_b - 1 - k) -
        gsl_sf_lnchoose(meth_a + unmeth_a + meth_b + unmeth_b - 2, 
                        meth_a + meth_b -  1);
}

double
ContingencyTable::beta_population_greater(
    const size_t meth_a, const size_t unmeth_a, 
    const size_t meth_b, const size_t unmeth_b) 
{
    double p = 0;
    
    for (size_t k = meth_b > unmeth_a ? meth_b - unmeth_a : 0;
         k < meth_b; ++k)
        p = log_sum_log(p, log_prob_hypergeo(
                            meth_a, unmeth_a, meth_b, unmeth_b, k));
    return exp(p);
}

