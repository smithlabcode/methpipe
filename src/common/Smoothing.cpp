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

#include "Smoothing.hpp"

#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iostream>
#include <stdexcept>

#include "smithlab_utils.hpp"

using std::vector;
using std::transform;
using std::divides;
using std::runtime_error;

static double
Epanechnikov_kernel(double i, double j, double bandwidth) {
  const double u = (j - i)/bandwidth;
  return 0.75*(1.0 - u*u);
}

void
KernelSmoothing(const double bandwidth,
                const vector<double> &x_vals,
                const vector<double> &y_vals,
                const vector<double> &x_target,
                vector<double> &y_target) {
  assert(x_vals.size() == y_vals.size());

  // allocate the space for the new y vals
  y_target.resize(x_target.size(), 0);

  // set the index into the current point to use for smoothing
  size_t x_start = 0;
  size_t x_end = 0;

  // iterate over the x target vals
  for (size_t i = 0; i < x_target.size(); ++i) {

    // calculate the x starting point
    while (x_start < x_vals.size() && x_vals[x_start] < x_target[i] - bandwidth)
      ++x_start;

    // calculate the x ending point
    while (x_end < x_vals.size() && x_vals[x_end] < x_target[i] + bandwidth)
      ++x_end;

    if (x_start >= x_end)
      throw runtime_error("smoothing using an interval of size 0");

    // set the number of points used for smoothing current value
    const size_t lim = x_end - x_start;

    // calculate the weights
    vector<double> weights(lim);
    for (size_t j = 0; j < lim; ++j)
      weights[j] = Epanechnikov_kernel(x_target[i], x_vals[x_start+j], bandwidth);
    const double weight_sum = accumulate(weights.begin(), weights.end(), 0.0);
    transform(weights.begin(), weights.end(), weights.begin(),
              [weight_sum] (const double w) {return w / weight_sum;});

    // apply the weights
    y_target[i] = 0;
    for (size_t j = 0; j < lim; ++j)
      y_target[i] += y_vals[x_start + j]*weights[j];
  }
}



void
KernelSmoothing(const double bandwidth, const vector<double> &y_vals,
                vector<double> &y_target) {

  // allocate the space for the new y vals
  y_target.resize(y_vals.size(), 0);

  // set the index into the current point to use for smoothing
  size_t x_start = 0;
  size_t x_end = 0;

  // iterate over the x target vals
  for (size_t i = 0; i < y_vals.size(); ++i) {

    // calculate the x starting point
    while (x_start < y_vals.size() && x_start < i - bandwidth)
      ++x_start;

    // calculate the x ending point
    while (x_end < y_vals.size() && x_end < i + bandwidth)
      ++x_end;

    if (x_start >= x_end)
      throw runtime_error("smoothing using an interval of size 0");

    // set the number of points used for smoothing current value
    const size_t lim = x_end - x_start;

    // calculate the weights
    vector<double> weights(lim);
    for (size_t j = 0; j < lim; ++j)
      weights[j] = Epanechnikov_kernel(i, x_start + j, bandwidth);

    const double weight_sum = accumulate(weights.begin(), weights.end(), 0.0);
    transform(weights.begin(), weights.end(), weights.begin(),
              [weight_sum] (const double w) {return w / weight_sum;});

    // apply the weights
    y_target[i] = 0;
    for (size_t j = 0; j < lim; ++j)
      y_target[i] += y_vals[x_start + j]*weights[j];
  }
}

#include <gsl/gsl_fit.h>

void
LocalLinearRegression(const double bandwidth,
                      const vector<double> &x_vals,
                      const vector<double> &y_vals,
                      const vector<double> &x_target,
                      vector<double> &y_target) {

  // Make sure the x and y vectors are of the same length
  assert(x_vals.size() == y_vals.size());

  // allocate the space for the new y vals
  y_target.resize(x_target.size(), 0);

  // set the index into the current point to use for smoothing
  size_t x_start = 0, x_end = 0;

  // iterate over the x target vals
  for (size_t i = 0; i < x_target.size(); ++i) {

    // calculate the x starting point
    while (x_start < x_vals.size() && x_vals[x_start] < x_target[i] - bandwidth)
      ++x_start;

    // calculate the x ending point
    while (x_end < x_vals.size() && x_vals[x_end] < x_target[i] + bandwidth)
      ++x_end;

    if (x_start >= x_end)
      throw runtime_error("smoothing using an interval of size 0");

    // set the number of points used for smoothing current value
    const size_t lim = x_end - x_start;

    // calculate the weights
    vector<double> weights(lim);
    for (size_t j = 0; j < lim; ++j)
      weights[j] = Epanechnikov_kernel(x_target[i], x_vals[x_start+j], bandwidth);
    const double weight_sum = accumulate(weights.begin(), weights.end(), 0.0);
    transform(weights.begin(), weights.end(), weights.begin(),
              [weight_sum] (const double w) {return w / weight_sum;});

    double intercept = 0, slope = 0;
    double c00 = 0, c10 = 0, c11 = 0;
    double ssq = 0;
    gsl_fit_wlinear(&x_vals[x_start], 1, &weights[0], 1, &y_vals[x_start], 1, lim,
                    &intercept, &slope, &c00, &c10, &c11, &ssq);
    y_target[i] = intercept + slope*x_target[i];
  }
}
