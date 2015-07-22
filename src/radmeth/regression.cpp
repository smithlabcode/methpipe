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
#include <stdexcept>

// GSL headers.
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_multimin.h>

#include "smithlab_utils.hpp"

#include "regression.hpp"

using std::istream; using std::string;
using std::vector; using std::istringstream;

std::istream&
operator>> (std::istream &is, Design &design) {
  string header_encoding;
  getline(is, header_encoding);

  istringstream header_is(header_encoding);
  string header_name;
  while (header_is >> header_name)
    design.factor_names.push_back(header_name);

  string row;
  while (getline(is, row)) {

    if (row.empty())
      continue;

    istringstream row_is(row);
    string token;
    row_is >> token;
    design.sample_names.push_back(token);

    vector<double> matrix_row;
    while (row_is >> token) {
      if (token.length() == 1 && (token == "0" || token == "1"))
        matrix_row.push_back(token == "1");
      else
        throw SMITHLABException("only binary factor levels are allowed:\n"
                                + row);
    }

    if (matrix_row.size() != design.factor_names.size())
      throw SMITHLABException("each row must have as many columns as "
            "factors:\n" + row);

    design.matrix.push_back(vector<double>());
    swap(design.matrix.back(), matrix_row);
  }
  return is;
}

std::ostream&
operator<< (std::ostream &os, const Design &design) {
  for(size_t factor = 0; factor < design.factor_names.size(); ++factor) {
    os << design.factor_names[factor];
    if (factor + 1 != design.factor_names.size())
      os << "\t";
  }
  os << std::endl;

  for(size_t sample = 0; sample < design.sample_names.size(); ++sample) {
    os << design.sample_names[sample] << "\t";
    for(size_t factor = 0; factor < design.factor_names.size(); ++factor) {
      os << design.matrix[sample][factor];
      if (factor + 1 != design.factor_names.size())
        os << "\t";
    }
      os << "\n";
  }
  return os;
}

void
remove_factor(Design &design, size_t factor) {
  design.factor_names.erase(design.factor_names.begin() + factor);
  for (size_t sample = 0; sample < design.sample_names.size(); ++sample)
    design.matrix[sample].erase(design.matrix[sample].begin() + factor);
}

// Parses a natural number from its string representation. Throws exception if
// the string does not encode one.
static size_t
parse_natural_number(string encoding) {
  istringstream iss(encoding);
  size_t number;
  iss >> number;
  if (!iss)
    throw SMITHLABException("The token \"" +encoding + "\" "
                            "does not encode a natural number");
  return number;
}

std::istream&
operator>>(std::istream &table_encoding, SiteProportions &props) {
  props.chrom.clear();
  props.position = 0;
  props.strand.clear();
  props.context.clear();
  props.meth.clear();
  props.total.clear();

  string row_encoding;
  getline(table_encoding, row_encoding);

  // Skip lines contining only the newline character (e.g. the last line of the
  // proportion table).
  if(row_encoding.empty())
    return table_encoding;

  // Get the row name (which must be specified like this: "chr:position") and
  // parse it.
  istringstream row_stream(row_encoding);
  string row_name_encoding;
  row_stream >> row_name_encoding;

  // Every row must start an identifier consisiting of genomic loci of the
  // corresponding site. Here we check this identifier has the correct number
  // of colons.
  const size_t num_colon =
            std::count(row_name_encoding.begin(), row_name_encoding.end(), ':');

  if (num_colon != 3)
    throw SMITHLABException("Each row in the count table must start with "
                            "a line chromosome:position:strand:context."
                            "Got \"" + row_name_encoding + "\" instead." );

  // First parse the row identifier.
  istringstream name_stream(row_name_encoding);
  getline(name_stream, props.chrom, ':');

  if (props.chrom.empty())
    throw SMITHLABException("Error parsing " + row_name_encoding +
                            ": chromosome name is missing.");

  string position_encoding;

  getline(name_stream, position_encoding, ':');
  props.position = parse_natural_number(position_encoding);
  getline(name_stream, props.strand, ':');
  getline(name_stream, props.context, ':');

  // After parsing the row identifier, parse count proportions.
  size_t total_count, meth_count;

  while (row_stream >> total_count >> meth_count) {
    props.total.push_back(total_count);
    props.meth.push_back(meth_count);
  }

  if (!row_stream.eof())
    throw SMITHLABException("Some row entries are not natural numbers: " +
                            row_stream.str());

  if (props.total.size() != props.meth.size())
    throw SMITHLABException("This row does not encode proportions"
                            "correctly:\n" + row_encoding);
  return table_encoding;
}

static double
pi(Regression *reg, size_t sample, const gsl_vector *parameters) {
  double dot_prod = 0;

  for(size_t factor = 0; factor < reg->design.factor_names.size(); ++factor)
    dot_prod +=
          reg->design.matrix[sample][factor]*gsl_vector_get(parameters, factor);

  double p = exp(dot_prod)/(1 + exp(dot_prod));

  return p;
}

static double
neg_loglik(const gsl_vector *parameters, void *object) {
  Regression *reg = (Regression *)(object);
  const size_t num_parameters = reg->design.factor_names.size() + 1;

  double log_lik = 0;

  //dispersion parameter phi is the last element of parameter vector
  const double dispersion_param = gsl_vector_get(parameters,
                                                  num_parameters - 1);
  const double phi = exp(dispersion_param)/(1 + exp(dispersion_param));

  for(size_t s = 0; s < reg->design.sample_names.size(); ++s) {
    const double n_s = reg->props.total[s];
    const double y_s = reg->props.meth[s];
    const double p_s = pi(reg, s, parameters);

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

  return (-1)*log_lik;
}

static void
neg_gradient(const gsl_vector *parameters, void *object,
                      gsl_vector *output) {

  Regression *reg = (Regression *)(object);
  const size_t num_parameters = reg->design.factor_names.size() + 1;

  const double dispersion_param = gsl_vector_get(parameters,
                                                  num_parameters - 1);

  const double phi = exp(dispersion_param)/(1 + exp(dispersion_param));

  for(size_t f = 0; f < num_parameters; ++f) {

    double deriv = 0;

    for(size_t s = 0; s < reg->design.sample_names.size(); ++s) {
      int n_s = reg->props.total[s];
      int y_s = reg->props.meth[s];
      double p_s = pi(reg, s, parameters);

      double term = 0;

      //a parameter linked to p
      if(f < reg->design.factor_names.size()) {
        double factor = (1 - phi)*p_s*(1 - p_s)*reg->design.matrix[s][f];
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

  gsl_vector_scale(output, -1.0);
}

static void
neg_loglik_and_grad(const gsl_vector *parameters,
                                void *object,
                                double *loglik_val,
                                gsl_vector *d_loglik_val) {

  *loglik_val = neg_loglik(parameters, object);
  neg_gradient(parameters, object, d_loglik_val);
}

bool
fit(Regression &r, vector<double> initial_parameters) {
  const size_t num_parameters = r.design.factor_names.size() + 1;

  if (initial_parameters.empty()) {
    for(size_t ind = 0; ind < num_parameters - 1; ++ind)
      initial_parameters.push_back(0.0);
    initial_parameters.push_back(-2.5);
  }

  if (initial_parameters.size() != num_parameters)
    throw std::runtime_error("Wrong number of initial parameters.");

  int status = 0;

  size_t iter = 0;

  gsl_multimin_function_fdf loglik_bundle;

  loglik_bundle.f = &neg_loglik;
  loglik_bundle.df = &neg_gradient;
  loglik_bundle.fdf = &neg_loglik_and_grad;
  loglik_bundle.n = num_parameters;
  loglik_bundle.params = (void *)&r;

  gsl_vector *parameters = gsl_vector_alloc(num_parameters);

  for (size_t parameter = 0; parameter < initial_parameters.size();
        ++parameter) {
    gsl_vector_set(parameters, parameter, initial_parameters[parameter]);
  }

  const gsl_multimin_fdfminimizer_type *T;

  //can also try gsl_multimin_fdfminimizer_conjugate_pr;
  T = gsl_multimin_fdfminimizer_conjugate_fr;

  gsl_multimin_fdfminimizer *s;
  s = gsl_multimin_fdfminimizer_alloc (T, num_parameters);

  gsl_multimin_fdfminimizer_set (s, &loglik_bundle, parameters, 0.001, 1e-4);

  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);

    if (status)
      break;

    status = gsl_multimin_test_gradient (s->gradient, 1e-4);
  }
  while (status == GSL_CONTINUE && iter < 700);
  //It it reasonable to reduce the number of iterations to 500?

  r.max_loglik = (-1)*neg_loglik(s->x, &r);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(parameters);

  return status == GSL_SUCCESS;
}
