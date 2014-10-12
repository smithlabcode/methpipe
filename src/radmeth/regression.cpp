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

#include <vector>
#include <algorithm>

#include "smithlab_utils.hpp"

#include "regression.hpp"

#include <gsl/gsl_multimin.h>

using std::istream; using std::string;
using std::vector; using std::istringstream;

void
Regression::set_response(const vector<size_t> &response_total,
                          const vector<size_t> &response_meth) {
  response_total_ = response_total;
  response_meth_ = response_meth;
}

double
Regression::p(size_t sample, const gsl_vector *parameters) const {
  double dot_prod = 0;

  for(size_t factor = 0; factor < design_.num_factors(); ++factor)
    dot_prod += design_(sample, factor)*gsl_vector_get(parameters, factor);

  double p = exp(dot_prod)/(1 + exp(dot_prod));

  return p;
}

double
Regression::loglik(const gsl_vector *parameters) const {
  double log_lik = 0;

  //dispersion parameter phi is the last element of parameter vector
  const double dispersion_param = gsl_vector_get(parameters,
                                                  num_parameters_ - 1);
  const double phi = exp(dispersion_param)/(1 + exp(dispersion_param));

  for(size_t s = 0; s < design_.num_samples(); ++s) {
    const double n_s = response_total_[s];
    const double y_s = response_meth_[s];
    const double p_s = p(s, parameters);

    for(int k = 0; k < y_s; ++k) {
      log_lik += log((1 - phi)*p_s + phi*k);
    }

    for(int k = 0; k < n_s - y_s; ++k) {
      log_lik += log((1 - phi)*(1 - p_s) + phi*k);
    }

    for(int k = 0; k < n_s; ++k) {
      log_lik -= log(1 + phi*(k - 1));
    }
  }

  return log_lik;
}

void
Regression::gradient(const gsl_vector *parameters, gsl_vector *output) const {

  const double dispersion_param = gsl_vector_get(parameters,
                                                  num_parameters_ - 1);

  const double phi = exp(dispersion_param)/(1 + exp(dispersion_param));

  for(size_t f = 0; f < num_parameters_; ++f) {

    double deriv = 0;

    for(size_t s = 0; s < design_.num_samples(); ++s) {
      int n_s = response_total_[s];
      int y_s = response_meth_[s];
      double p_s = p(s, parameters);

      double term = 0;

      //a parameter linked to p
      if(f < design_.num_factors()) {
        double factor = (1 - phi)*p_s*(1 - p_s)*design_(s, f);
        if (factor == 0) continue;

        for(int k = 0; k < y_s; ++k)
          term += 1/((1 - phi)*p_s + phi*k);

        for(int k = 0; k < n_s - y_s; ++k)
          term -= 1/((1 - phi)*(1 - p_s) + phi*k);

        deriv += term*factor;
      } else { // the parameter linked to phi
        for(int k = 0; k < y_s; ++k)
          term += (k - p_s)/((1 - phi)*p_s + phi*k);

        for(int k = 0; k < n_s - y_s; ++k)
          term += (k - (1 - p_s))/((1 - phi)*(1 - p_s) + phi*k);

        for(int k = 0; k < n_s; ++k) {
          term -= (k - 1)/(1 + phi*(k - 1));
        }

        deriv += term * phi * (1 - phi);
      }
    }

    gsl_vector_set(output, f, deriv);
  }
}

std::vector<double>
Regression::fitted_distribution_parameters() {
  vector<double> distribution_parameters;

	for (size_t sample = 0; sample < design_.num_samples(); ++sample) {
    double sample_parameter = 0;
    for (size_t factor = 0; factor < design_.num_factors(); ++factor) {
      sample_parameter += fitted_parameters_[factor]*design_(sample, factor);
    }
    distribution_parameters.push_back(
            exp(sample_parameter)/(1 + exp(sample_parameter)) );
	}

	double phi_param = fitted_parameters_.back();
	double phi = exp(phi_param)/(1 + exp(phi_param));
	distribution_parameters.push_back(phi);

  return distribution_parameters;
}

double
Regression::min_methdiff(size_t test_factor) {

  vector<double> methdiffs;

  for (size_t sample = 0; sample < design_.num_samples(); ++sample) {
    if (design_(sample, test_factor) == 1) {
      double full_parameter_sum = 0;
      double reduced_parameter_sum = 0;
      for (size_t factor = 0; factor < design_.num_factors(); ++factor) {
        full_parameter_sum +=
                          fitted_parameters_[factor]*design_(sample, factor);
      }
      reduced_parameter_sum = full_parameter_sum -
                              fitted_parameters_[test_factor];
      const double full_pi =
                          exp(full_parameter_sum)/(1 + exp(full_parameter_sum));
      const double reduced_pi =
                    exp(reduced_parameter_sum)/(1 + exp(reduced_parameter_sum));
      methdiffs.push_back(fabs(full_pi - reduced_pi));
    }
  }

  return *std::min_element(methdiffs.begin(), methdiffs.end());
}

double
Regression::log_fold_change(size_t factor) {
  return fitted_parameters_[factor];
}



Design::Design(istream &is) {
  string header_encoding;
  getline(is, header_encoding);

  istringstream header_is(header_encoding);
  string header_name;
  while (header_is >> header_name)
    factor_names_.push_back(header_name);

  string row;
  while (getline(is, row)) {

    if (row.empty())
      continue;

    istringstream row_is(row);
    string token;
    row_is >> token;
    sample_names_.push_back(token);

    vector<double> matrix_row;
    while (row_is >> token) {
      if (token.length() == 1 && (token == "0" || token == "1"))
        matrix_row.push_back(token == "1");
      else
        throw SMITHLABException("only binary factor levels are allowed:\n"
                                + row);
    }

    if (matrix_row.size() != num_factors())
      throw SMITHLABException("each row must have as many columns as "
            "factors:\n" + row);

    matrix_.push_back(vector<double>());
    swap(matrix_.back(), matrix_row);
  }
}

double
Design::operator() (size_t sample, size_t factor) const {
  return matrix_[sample][factor];
}

std::ostream&
operator<<(std::ostream& os, const Design &design) {

  for(size_t factor = 0; factor < design.num_factors(); ++factor) {
    os << design.factor_names_[factor];
    if (factor + 1 != design.num_factors())
      os << "\t";
  }
  os << std::endl;

  for(size_t sample = 0; sample < design.num_samples(); ++sample) {
    os << design.sample_names_[sample] << "\t";
    for(size_t factor = 0; factor < design.num_factors(); ++factor) {
      os << design.matrix_[sample][factor];
      if (factor + 1 != design.num_factors())
        os << "\t";
    }
      os << "\n";
  }
  return os;
}

void
Design::remove_factor(size_t factor) {
  factor_names_.erase(factor_names_.begin() + factor);

  for (size_t sample = 0; sample < num_samples(); ++sample)
    matrix_[sample].erase(matrix_[sample].begin() + factor);
}

void
Design::remove_factor_name(std::string name) {
  vector<string>::const_iterator name_it =
                std::find(factor_names_.begin(), factor_names_.end(), name);

  if (name_it == factor_names_.end())
    throw SMITHLABException(name + " is not one of the factor names");

  size_t factor =  name_it - factor_names_.begin();

  remove_factor(factor);
}

bool
Design::operator== (const Design &other_design) const {
  return (factor_names_ == other_design.factor_names_) &&
          (sample_names_ == other_design.sample_names_) &&
          (matrix_ == other_design.matrix_);
}

bool
Design::operator!= (const Design &other_design) const {
  return !(*this == other_design);
}
