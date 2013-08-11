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

#ifndef FALSE_DISCOVERY_RATE
#define FALSE_DISCOVERY_RATE

#include <vector>

namespace FDR
{
    double
    get_empirical_p_value(const std::vector<double> &random_scores, 
                          const double &observed_score);
    void
    assign_empirical_p_values(const std::vector<double> &random_scores, 
                              const std::vector<double> &observed_scores, 
                              std::vector<double> &p_values);
    double
    get_fdr_cutoff(const std::vector<double> &p_values, const double fdr);
};

#endif
