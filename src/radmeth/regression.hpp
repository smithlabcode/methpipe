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
#include <sstream>
#include <string>

// GSL headers.
#include <gsl/gsl_matrix_double.h>

// Objects of this class represent design matrices.
//
// A Design object can be created like this:
//
//  string encoding = "f1\tf2\ns1\t1\t1\ns2\t1\t0";
//  istringstream iss(encoding);
//  Design design(iss);
//
// To drop the second factor from the desin:
//
//  design.remove_factor(1);

class Design {
public:
  Design() {}
  Design(std::istream &is);
  size_t num_factors() const { return factor_names_.size(); }
  size_t num_samples() const { return sample_names_.size(); }
  std::vector<std::string> factor_names() const { return factor_names_; }
  std::vector<std::string> sample_names() const { return sample_names_; }
  std::vector<std::vector<double> > matrix() const { return matrix_; }
  double operator() (size_t sample, size_t factor) const;
  void remove_factor(size_t factor);
  friend std::ostream& operator<<(std::ostream& os, const Design &design);
private:
  std::vector<std::string> factor_names_;
  std::vector<std::string> sample_names_;
  std::vector<std::vector<double> > matrix_;
};

class Regression {
public:
  Regression() : num_parameters_(0), maximum_likelihood_(0) {}
  Regression(const Design &design)
    : design_(design), num_parameters_(design_.num_factors() + 1),
      maximum_likelihood_(0) {}
  void set_response(const std::vector<size_t> &response_total,
                    const std::vector<size_t> &response_meth);
  double loglik(const gsl_vector *parameters) const;
  std::vector<size_t> response_total() const { return response_total_; }
  std::vector<size_t> response_meth() const { return response_meth_; }
  double p(size_t sample, const gsl_vector *v) const;
  std::vector<double> fitted_parameters() { return fitted_parameters_; }
  std::vector<double> fitted_distribution_parameters();
  double maximum_likelihood() { return maximum_likelihood_; }
  double min_methdiff(size_t factor);
  double log_fold_change(size_t factor);
  friend bool gsl_fitter(Regression &r, std::vector<double> initial_parameters);
  void gradient(const gsl_vector *parameters, gsl_vector *output) const;
private:
  Design design_;
  std::vector<size_t> response_total_;
  std::vector<size_t> response_meth_;
  const size_t num_parameters_;
  std::vector<double> fitted_parameters_;
  double maximum_likelihood_;
};

bool gsl_fitter(Regression &r,
              std::vector<double> initial_parameters = std::vector<double> ());

#endif //REGRESSION_HPP_
