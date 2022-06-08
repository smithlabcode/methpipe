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

#include "Distro.hpp"
#include "smithlab_utils.hpp"

#include <gsl/gsl_sf_gamma.h>

#include <cmath>
#include <algorithm>
#include <sstream>
#include <iterator>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <iostream>
#include <limits>

using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::max;
using std::min;
using std::bind2nd;
using std::divides;
using std::numeric_limits;
using std::pair;
using std::make_pair;
using std::accumulate;
using std::max_element;


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// DISTRO__

double
Distro_::log_sum_log_vec(const vector<double> &vals, size_t limit) {
  const vector<double>::const_iterator x =
    max_element(vals.begin(), vals.begin() + limit);
  const double max_val = *x;
  const size_t max_idx = x - vals.begin();
  double sum = 1.0;
  for (size_t i = 0; i < limit; ++i) {
    if (i != max_idx) {
      sum += exp(vals[i] - max_val);
#ifdef DEBUG
      assert(std::isfinite(sum));
#endif
    }
  }
  return max_val + log(sum);
}

// Distro_::Distro_(const Distro_ &other) :
//   params(other.params),
//   workspace_vals(other.workspace_vals),
//   workspace_probs(other.workspace_probs) {
//   rng = gsl_rng_clone(other.rng);
// }

Distro_::Distro_() {
}

Distro_::Distro_(const std::vector<double> p) : params(p) {
}

Distro_::~Distro_() {}


string
Distro_::tostring() const {
  std::ostringstream os;
  if (!params.empty())
    os << std::setprecision(4) << params.front();
  for (size_t i = 1; i < params.size(); ++i)
    os << " " << std::setprecision(4) << params[i];
  return os.str();
}


double
Distro_::log_likelihood(std::vector<double>::const_iterator a,
                        std::vector<double>::const_iterator b) const {
  double l = 0;
  for (; a < b; ++a)
    l += this->log_likelihood(*a);
  return l;
}

double
Distro_::operator()(const double val) const {
  return exp(log_likelihood(val));
}


double
Distro_::operator()(const vector<double> &vals) const {
  const size_t lim = vals.size();
  double l = 1;
  for (size_t i = 0; i < lim; ++i)
    l *= this->operator()(vals[i]);
  return l;
}


double
Distro_::log_likelihood(const std::vector<double> &vals) const {
  double l = 0;
  const size_t lim = vals.size();
  for (size_t i = 0; i < lim; ++i)
    l += this->log_likelihood(vals[i]);
  return l;
}

double
Distro_::log_likelihood(const std::vector<double> &vals,
                        const std::vector<double> &scales) const {
  double l = 0;
  const size_t lim = vals.size();
  for (size_t i = 0; i < lim; ++i)
    l += this->log_likelihood(vals[i], scales[i]);
  return l;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// DISTRO

bool
Distro::has_params(const string &name) {
  return (smithlab::split(name, ",").size() >= 2);
}

Distro::Distro(const std::string &n, const std::string &params) :
  name(n), d(distro_factory(n, params)) {}

Distro::Distro(const std::string &n, const std::vector<double> &params) :
  name(n), d(distro_factory(name)) {d->set_params(params);}

Distro::Distro(const std::string &s) : d(distro_factory(s)) {
  vector<string> name_split;
  if (s.find(",") == string::npos)  // whitespaces seperated
    name_split = smithlab::split_whitespace_quoted(s);
  else // comma seperated
    name_split = smithlab::split(s, ",");
  name = name_split[0];
  name = name_split[0];
}

Distro::Distro(const Distro &rhs) : name(rhs.name), d(distro_factory(rhs.name)) {
  vector<double> tmp_params(rhs.get_params());
  d->set_params(tmp_params);
}

Distro&
Distro::operator=(const Distro &rhs) {
  if (this != &rhs) {
    name = rhs.name;
    d = distro_factory(rhs.name);
    vector<double> tmp_params(rhs.get_params());
    d->set_params(tmp_params);
  }
  return *this;
}

Distro::~Distro() {if (d) delete d;}

double
Distro::operator()(double val) const {return (*d)(val);}

double
Distro::operator()(const std::vector<double> &vals) const {
  return (*d)(vals);
}

void
Distro::estimate_params_ml(const std::vector<double> &vals) {
  d->estimate_params_ml(vals);
}

void
Distro::estimate_params_ml(const std::vector<double> &vals,
                           const std::vector<double> &scales,
                           const std::vector<double> &probs) {
  d->estimate_params_ml(vals, scales, probs);
}

void
Distro::estimate_params_ml(const std::vector<double> &vals,
                           const std::vector<double> &weights) {
  d->estimate_params_ml(vals, weights);
}

double
Distro::log_likelihood(const std::vector<double> &vals) const {
  return d->log_likelihood(vals);
}

double
Distro::log_likelihood(const std::vector<double> &vals,
                       const std::vector<double> &scales) const {
  return d->log_likelihood(vals, scales);
}



double
Distro::log_likelihood(std::vector<double>::const_iterator a,
                       std::vector<double>::const_iterator b) const {
  return d->log_likelihood(a, b);
}

double
Distro::log_likelihood(double val) const {return d->log_likelihood(val);}

double
Distro::log_likelihood(const double &val,
                       const double &scale) const
{return d->log_likelihood(val, scale);}


std::string
Distro::tostring() const {return name + string(" ") + d->tostring();}

std::ostream&
operator<<(std::ostream& s, const Distro& distro) {
  return s << distro.tostring();
}

double
Distro::log_sum_log_vec(const vector<double> &vals, size_t limit) {
  return Distro_::log_sum_log_vec(vals, limit);
}

////////////////////////////////////////////////////////////////////////
//
// DISTRO FACTORY


Distro_ *
distro_factory(string name, string params) {
  Distro_ *distro;
  if (name == "exp")
    distro = new ExpDistro();
  else if (name == "pois")
    distro = new PoisDistro();
  else if (name == "nbd")
    distro = new NegBinomDistro();
  else throw SMITHLABException("bad distribution name \"" + name + "\"");

  vector<string> params_split = smithlab::split(params, ",");
  if (params_split.size() != distro->required_params())
    throw SMITHLABException("bad number of params: " + smithlab::toa(params_split.size()) +
                            " for distro " + name + "\"");
  else {
    vector<double> params_vec;
    for (size_t i = 0; i < params_split.size(); ++i)
      params_vec.push_back(atof(params_split[i].c_str()));
    distro->set_params(params_vec);
  }
  return distro;
}


Distro_ *
distro_factory(string name_arg) {
  vector<string> name_split;
  if (name_arg.find(",") == string::npos)  // whitespaces seperated
    name_split = smithlab::split_whitespace_quoted(name_arg);
  else // comma seperated
    name_split = smithlab::split(name_arg, ",");

  const string name = name_split.front();

  Distro_ *distro;
  if (name == "exp")
    distro = new ExpDistro();
  else if (name == "pois")
    distro = new PoisDistro();
  else if (name == "nbd")
    distro = new NegBinomDistro();
  else throw SMITHLABException("bad distribution name \"" + name + "\"");

  if (name_split.size() > 1) {

    vector<string> params_split(vector<string>(name_split.begin() + 1,
                                               name_split.end()));
    if (params_split.size() != distro->required_params())
      throw SMITHLABException("bad number of params: " + smithlab::toa(params_split.size()) +
                              " for distro " + name + "\"");
    else {
      vector<double> params_vec;
      for (size_t i = 0; i < params_split.size(); ++i)
        params_vec.push_back(atof(params_split[i].c_str()));
      distro->set_params(params_vec);
    }
  }
  return distro;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



double
ExpDistro::log_likelihood(const double val) const {
  return -log(params[0]) - val/params[0];
}

double
ExpDistro::log_likelihood(const double &val,
                          const double &scale) const
{
  //// TEST NEEDED
  return -log(params[0] * scale) - val / params[0] / scale;
}


ExpDistro::ExpDistro(const ExpDistro &rhs) :
  Distro_(rhs.params) {}

ExpDistro&
ExpDistro::operator=(const ExpDistro &rhs) {
  if (this != &rhs) {
    this->Distro_::params = rhs.Distro_::params;
  }
  return *this;
}


void
ExpDistro::estimate_params_ml(const vector<double> &vals) {
  params.front() = accumulate(vals.begin(), vals.end(), 0.0)/vals.size();
}

void
ExpDistro::estimate_params_ml(const vector<double> &vals,
                              const vector<double> &scales,
                              const vector<double> &probs) {
  cerr << "ExpDistro::estimate_params_ml:  "
       << "Using an untested function" << endl;

  vector<double> values(vals.begin(), vals.end());
  for (size_t i = 0; i < values.size(); ++i)
    values[i] /= scales[i];
  if (probs.size() == 0)
    estimate_params_ml(vals);
  else
    estimate_params_ml(vals, probs);
}


void
ExpDistro::estimate_params_ml(const vector<double> &vals,
                              const vector<double> &probs) {
  const size_t lim = vals.size();
  if (workspace_vals.size() < lim) {
    workspace_vals.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = log(probs[i]);
    workspace_vals[i] = log(vals[i]) + log(probs[i]);
  }
  const double prob_sum = exp(log_sum_log_vec(workspace_probs, lim));
  params.front() = exp(log_sum_log_vec(workspace_vals, lim))/prob_sum;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////




double
PoisDistro::log_likelihood(const double val) const {
  return -params.front() + val*log(params.front()) -
    gsl_sf_lnfact(static_cast<size_t>(val));
}

double
PoisDistro::log_likelihood(const double &val,
                           const double &scale) const
{
  const double lambda = params[0] * scale;
  return -lambda + val*log(lambda) -
    gsl_sf_lnfact(static_cast<size_t>(val));
}




PoisDistro::PoisDistro(const PoisDistro &rhs) : Distro_(rhs.params) {}


PoisDistro&
PoisDistro::operator=(const PoisDistro &rhs) {
  if (this != &rhs) {
    this->Distro_::params = rhs.Distro_::params;
  }
  return *this;
}


void
PoisDistro::estimate_params_ml(const vector<double> &vals) {
  params.front() = accumulate(vals.begin(), vals.end(), 0.0)/vals.size();
}


void
PoisDistro::estimate_params_ml(const vector<double> &vals,
                               const vector<double> &probs) {
  const size_t lim = vals.size();
  if (workspace_vals.size() < lim) {
    workspace_vals.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = log(probs[i]);
    workspace_vals[i] = log(vals[i]) + log(probs[i]);
  }
  const double prob_sum = exp(log_sum_log_vec(workspace_probs, lim));
  params.front() = exp(log_sum_log_vec(workspace_vals, lim))/prob_sum;
}

void
PoisDistro::estimate_params_ml(const vector<double> &vals,
                               const vector<double> &scales,
                               const vector<double> &probs) {
  const size_t lim = vals.size();
  if (workspace_vals.size() < lim) {
    workspace_vals.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i)
    {
      workspace_probs[i] = log(probs[i]) + log(scales[i]);
      workspace_vals[i] = log(vals[i]) + log(probs[i]);
    }
  const double prob_sum = exp(log_sum_log_vec(workspace_probs, lim));
  params.front() = exp(log_sum_log_vec(workspace_vals, lim))/prob_sum;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

const double
NegBinomDistro::max_allowed_alpha = 100;

const double
NegBinomDistro::min_allowed_alpha = 1e-20;

const double
NegBinomDistro::alpha_allowed_error = 1e-10;


void
NegBinomDistro::set_helpers() {
  n_helper = 1/params[1];
  p_helper = n_helper/(n_helper + params[0]);
  n_log_p_minus_lngamma_n_helper = n_helper*log(p_helper) -
    gsl_sf_lngamma(n_helper);
  log_q_helper = log(1 - p_helper);
  //TODO: should check that these are valid!!!
}


void
NegBinomDistro::set_params(const std::vector<double> &p) {
  Distro_::set_params(p);
  set_helpers();
}




double
NegBinomDistro::log_likelihood(const double val) const {
  const double P = (gsl_sf_lngamma(val + n_helper) -
                    gsl_sf_lnfact(static_cast<size_t>(val))) +
    n_log_p_minus_lngamma_n_helper + val*log_q_helper;
  if (!std::isfinite(P))
    return -40;
  return P;
}

double
NegBinomDistro::log_likelihood(const double &val,
                               const double &scale) const
{
  //// TEST NEEDED
  // alpha is not scaled
  const double scaled_p_helper = n_helper / (n_helper + params[0] * scale);
  const double scaled_n_log_p_minus_lngamma_n_helper =
    n_helper*log(scaled_p_helper) - gsl_sf_lngamma(n_helper);
  const double scaled_log_q_helper = log(1 - scaled_p_helper);


  const double P = (gsl_sf_lngamma(val + n_helper) -
                    gsl_sf_lnfact(static_cast<size_t>(val))) +
    scaled_n_log_p_minus_lngamma_n_helper + val * scaled_log_q_helper;
  if (!std::isfinite(P))
    return -40;
  return P;
}



NegBinomDistro::NegBinomDistro(const NegBinomDistro &rhs) :
  Distro_(rhs.params) {
  set_helpers();
}


NegBinomDistro&
NegBinomDistro::operator=(const NegBinomDistro &rhs) {
  if (this != &rhs) {
    this->Distro_::params = rhs.Distro_::params;
    set_helpers();
  }
  return *this;
}


static inline double
score_fun_first_term(const vector<double> &vals_hist,
                     const double mu, const double alpha) {
  double sum = 0;
  for (size_t i = 0; i < vals_hist.size(); ++i)
    if (vals_hist[i] > 0) {
      double inner_sum = 0;
      for (size_t j = 0; j < i; ++j)
        inner_sum += j/(1 + alpha*j);
      sum += vals_hist[i]*inner_sum;
    }
  return sum;
}


static inline double
alpha_score_function(const vector<double> &vals_hist, const double mu,
                     const double alpha, const double vals_count) {
  const double one_plus_alpha_mu = 1 + alpha*mu;
  return (score_fun_first_term(vals_hist, mu, alpha)/vals_count +
          (log(one_plus_alpha_mu)/alpha - mu)/alpha);
}


void
NegBinomDistro::estimate_params_ml(const vector<double> &vals) {
  // This is the mu
  params.front() = std::accumulate(vals.begin(), vals.end(), 0.0)/vals.size();

  // Now for the alpha
  const double max_value = *std::max_element(vals.begin(), vals.end());
  vector<double> vals_hist(static_cast<size_t>(max_value) + 1, 0.0);
  for (size_t i = 0; i < vals.size(); ++i)
    ++vals_hist[static_cast<size_t>(vals[i])];

  const double vals_count = vals.size();

  const double mu = params.front();
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;

  double a_mid = max_allowed_alpha;
  double diff = std::numeric_limits<double>::max();
  double prev_val = std::numeric_limits<double>::max();
  while (diff > alpha_allowed_error && fabs((a_high - a_low)/max(a_high, a_low)) > alpha_allowed_error) {
    a_mid = (a_low + a_high)/2;
    const double mid_val = alpha_score_function(vals_hist, mu, a_mid, vals_count);
    if (mid_val < 0)
      a_high = a_mid;
    else
      a_low = a_mid;
    diff = std::fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  params[1] = a_mid;

  set_helpers();
  //   const size_t lim = vals.size();

  //   // This is the mu
  //   params.front() = std::accumulate(vals.begin(), vals.begin() + lim, 0.0)/lim;
  //   const double mu = params.front();
  //   const double var = gsl_stats_variance_m(&vals.front(), 1, vals.size(), mu);
  //   const double r = (mu*mu)/(var - mu);
  //   // const double p = r/(r + params[0]);
  //   params[1] = max(0.01, 1/r);

  //   set_helpers();
}

void
NegBinomDistro::estimate_params_ml(const vector<double> &vals,
                                   const vector<double> &probs) {
  //   const size_t lim = vals.size();
  //   if (workspace_vals.size() < lim) {
  //     workspace_vals.resize(lim);
  //     workspace_probs.resize(lim);
  //   }
  //   for (size_t i = 0; i < lim; ++i) {
  //     workspace_probs[i] = log(probs[i]);
  //     workspace_vals[i] = log(vals[i]) + log(probs[i]);
  //   }
  //   const double vals_count = exp(log_sum_log_vec(workspace_probs, lim));
  //   const double mu = exp(log_sum_log_vec(workspace_vals, lim))/vals_count;
  //   const double var = gsl_stats_wvariance_m(&vals.front(), 1, &probs.front(), 1,
  //                                       vals.size(), mu);
  //   const double r = (mu*mu)/(var - mu);
  //   // const double p = r/(r + params[0]);
  //   params[0] = mu;
  //   params[1] = max(0.01, 1/r);
  //   set_helpers();
  const size_t lim = vals.size();
  if (workspace_vals.size() < lim) {
    workspace_vals.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = log(probs[i]);// - centering_value;
    workspace_vals[i] = log(vals[i]) + log(probs[i]);// - centering_value;
  }

  const double vals_count = exp(log_sum_log_vec(workspace_probs, lim));
  params.front() = exp(log_sum_log_vec(workspace_vals, lim))/vals_count;

  // Now for the alpha
  const double max_value = *std::max_element(vals.begin(), vals.begin() + lim);
  vector<double> vals_hist(static_cast<size_t>(max_value) + 1, 0.0);
  for (size_t i = 0; i < lim; ++i)
    vals_hist[static_cast<size_t>(vals[i])] += probs[i];

  const double mu = params.front();
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;

  double a_mid = max_allowed_alpha;
  double diff = std::numeric_limits<double>::max();
  double prev_val = std::numeric_limits<double>::max();
  while (diff > alpha_allowed_error &&
         fabs((a_high - a_low)/max(a_high, a_low)) > alpha_allowed_error) {
    a_mid = (a_low + a_high)/2;
    const double mid_val = alpha_score_function(vals_hist, mu, a_mid, vals_count);
    if (mid_val < 0)
      a_high = a_mid;
    else
      a_low = a_mid;
    //     cerr << diff << "\t"
    //   << prev_val << "\t"
    //   << mid_val << "\t"
    //   << a_low << "\t"
    //   << a_mid << "\t"
    //   << a_high << endl;
    diff = std::fabs((prev_val - mid_val)/std::max(mid_val, prev_val));
    prev_val = mid_val;
  }
  params[1] = a_mid;

  set_helpers();
}

static double
llh_derivative_rt_alpha(const vector<double> &vals,
                        const vector<double> &scales,
                        const vector<double> &probs,
                        const vector<double> &vals_hist,
                        const double mu,
                        const double alpha)
{
  const double first_term = score_fun_first_term(vals_hist, mu, alpha);

  const double mu_times_alpha = mu * alpha;
  const double alpha_inverse = 1 / alpha;
  const double alpha_square_inverse = pow(alpha_inverse, 2.0);

  double second_term = 0;
  for (size_t i = 0; i < vals.size(); ++i)
    {
      const double one_plus_extra = 1 + scales[i] * mu_times_alpha;
      second_term += probs[i] *
        (alpha_square_inverse * log(one_plus_extra) -
         scales[i] * mu * (alpha_inverse + vals[i]) / one_plus_extra);
    }

  return first_term + second_term;
}

void
NegBinomDistro::estimate_params_ml(const std::vector<double> &vals,
                                   const std::vector<double> &scales,
                                   const std::vector<double> &probs)
{
  const size_t lim = vals.size();
  if (workspace_vals.size() < lim)
    {
      workspace_vals.resize(lim);
      workspace_probs.resize(lim);
    }
  for (size_t i = 0; i < lim; ++i)
    {
      workspace_probs[i] = log(probs[i]) + log(scales[i]);
      workspace_vals[i] = log(vals[i]) + log(probs[i]);
    }

  // this is mu
  const double mu = exp(log_sum_log_vec(workspace_vals, lim)) /
    exp(log_sum_log_vec(workspace_probs, lim));

  // Now for the alpha
  const double max_value = *std::max_element(vals.begin(), vals.end());
  vector<double> vals_hist(static_cast<size_t>(max_value) + 1, 0.0);
  for (size_t i = 0; i < vals.size(); ++i)
    vals_hist[static_cast<size_t>(vals[i])] += probs[i];

  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;

  double a_mid = max_allowed_alpha;
  double diff = std::numeric_limits<double>::max();
  double prev_val = std::numeric_limits<double>::max();
  while (diff > alpha_allowed_error &&
         fabs((a_high - a_low) / max(a_high, a_low)) > alpha_allowed_error)
    {
      a_mid = (a_low + a_high)/2;
      const double mid_val =
        llh_derivative_rt_alpha(vals, scales, probs, vals_hist, mu, a_mid);

      if (mid_val < 0)
        a_high = a_mid;
      else
        a_low = a_mid;

      diff = std::fabs((prev_val - mid_val)/prev_val);
      prev_val = mid_val;
    }

  params[0] = mu;
  params[1] = a_mid;

  set_helpers();
}
