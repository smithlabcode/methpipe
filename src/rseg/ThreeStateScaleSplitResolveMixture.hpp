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

#ifndef THREE_STATE_SCALE_SPLIT_RESOLVE_MIXTURE_HPP
#define THREE_STATE_SCALE_SPLIT_RESOLVE_MIXTURE_HPP

#include "smithlab_utils.hpp"
#include "SplitDistro.hpp"

void
ThreeStateScaleSplitResolveMixture(const std::vector<double> &values,
			      const std::vector<double> &vals_a,
			      const std::vector<double> &vals_b,
                  const std::vector<double> &scales,
			      const size_t max_iterations,
			      const double tolerance,
			      int VERBOSE,
			      SplitDistro &fg_distro,
			      SplitDistro &mid_distro,
			      SplitDistro &bg_distro,
			      std::vector<double> &mixing);
#endif
