/*
  Copyright (C) 2008 Cold Spring Harbor Laboratory
  Authors: Andrew D. Smith

  This file is part of rmap.

  rmap is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  rmap is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with rmap; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef SMOOTHING_HPP
#define SMOOTHING_HPP

#include <vector>

void
KernelSmoothing(const double bandwidth,
		const std::vector<double> &x_values, 
		const std::vector<double> &y_values,
		const std::vector<double> &x_target, 
		std::vector<double> &y_target);

void
LocalLinearRegression(const double bandwidth,
		      const std::vector<double> &x_values, 
		      const std::vector<double> &y_values,
		      const std::vector<double> &x_target, 
		      std::vector<double> &y_target);


void
KernelSmoothing(const double bandwidth,
		const std::vector<double> &y_vals,
		std::vector<double> &y_target);


#endif
