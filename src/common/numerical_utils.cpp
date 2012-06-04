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

#include "numerical_utils.hpp"

#include <cmath>
#include <vector>
#include <algorithm>

using std::vector;

double
log_sum_log_vec(const std::vector<double> &vals, const size_t limit) 
{
    const std::vector<double>::const_iterator x = 
        std::max_element(vals.begin(), vals.begin() + limit);
    const double max_val = *x;
    const size_t max_idx = x - vals.begin();
    double sum = 1.0;
    for (size_t i = 0; i < limit; ++i) 
    {
        if (i != max_idx) 
        {
            sum += exp(vals[i] - max_val);
        }
    }
    return max_val + log(sum);
}

double
log_sum_log(const std::vector<double>::const_iterator &begin,
            const std::vector<double>::const_iterator &end)
{
    const std::vector<double>::const_iterator max_itr = 
        std::max_element(begin, end);
    const double max_val = *max_itr;

    double sum = 1.0;
    for (std::vector<double>::const_iterator itr = begin; itr < end; ++itr) 
        if (itr != max_itr) sum += exp(*itr - max_val);
    
    return max_val + log(sum);
}


