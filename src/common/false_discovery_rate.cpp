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

#include <cassert>

#include <vector>
#include <algorithm>
#include <functional>

#include "false_discovery_rate.hpp"

using std::vector;
using std::upper_bound;

double
FDR::get_empirical_p_value(const vector<double> &random_scores, 
                           const double &observed_score)
{
    return random_scores.size() == 0 ? 0 :
        (random_scores.end() - 
         upper_bound(random_scores.begin(), random_scores.end(), observed_score))
        / static_cast<double>(random_scores.size());
}


void
FDR::assign_empirical_p_values(
    const vector<double> &random_scores, 
    const vector<double> &observed_scores, 
    vector<double> &p_values) 
{
    // make sure random_scores are sorted
    assert(std::adjacent_find(random_scores.begin(), random_scores.end(),
                              std::greater<double>())
           == random_scores.end()); 

    // get p_values
    p_values.resize(observed_scores.size());
    for (size_t i = 0; i < observed_scores.size(); ++i)
        p_values[i] = get_empirical_p_value(random_scores, observed_scores[i]);
    
    // std::transform(observed_scores.begin(), observed_scores.end(),
    //                p_values.begin(),
    //                std::bind1st(std::ptr_fun(get_empirical_p_value),
    //                             random_scores)); 
}

double
FDR::get_fdr_cutoff(const vector<double> &p_values, const double fdr) 
{
    if (fdr < 0) return 0;
    else if (fdr > 1) return 1;

    vector<double> local(p_values);
    std::sort(local.begin(), local.end());
    assert(local.size() > 0);
    size_t i = 0;
    for (; i < local.size() - 1 && 
             local[i+1] < fdr*static_cast<double>(i+1)/local.size(); ++i);
    assert(i < local.size());
    return local[i];
}

