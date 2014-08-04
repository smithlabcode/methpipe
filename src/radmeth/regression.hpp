/*    Copyright (C) 2013 University of Southern California and
 *                       Egor Dolzhenko
 *                       Andrew D Smith
 *
 *    Authors: Andrew D. Smith and Egor Dolzhenko
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
 */

// Declares the Regression class which implements beta-binomial regression.

#ifndef REGRESSION_HPP_
#define REGRESSION_HPP_

// Std headers.
#include <vector>

// GSL headers.
#include <gsl/gsl_matrix_double.h>

// Local headers.
#include "design.hpp"

class Regression {
public:
  Regression() : num_parameters_(0), maximum_likelihood_(0) {}
  Regression(const Design &design) 
    : design_(design), num_parameters_(design_.num_factors() + 1), 
      maximum_likelihood_(0) {}
  Design design() const { return design_; }
  void set_response(const std::vector<size_t> &response_total, 
                    const std::vector<size_t> &response_meth);
  std::vector<size_t> response_total() const { return response_total_; }
  std::vector<size_t> response_meth() const {return response_meth_; }
  double p(size_t sample, const gsl_vector *v) const;
  double loglik(const gsl_vector *parameters) const;
  void gradient(const gsl_vector *parameters, gsl_vector *output) const;
  std::vector<double> fitted_parameters() { return fitted_parameters_; }
  std::vector<double> fitted_distribution_parameters();
  double maximum_likelihood() { return maximum_likelihood_; }
  double min_methdiff(size_t factor);
  double log_fold_change(size_t factor);
  friend bool gsl_fitter(Regression &r, std::vector<double> initial_parameters);
private:
  Design design_;
  std::vector<size_t> response_total_;
  std::vector<size_t> response_meth_;
  const size_t num_parameters_;
  std::vector<double> fitted_parameters_;
  double maximum_likelihood_;
};

#endif //REGRESSION_HPP_
