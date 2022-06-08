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

#include "SplitDistro.hpp"
#include "smithlab_utils.hpp"

#include <cmath>
#include <algorithm>
#include <sstream>
#include <iterator>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <iostream>
#include <limits>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_bessel.h>

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
using std::tr1::unordered_map;

static double
split_distro_log_sum_log_vec(const vector<double> &vals, size_t limit) {
  const vector<double>::const_iterator x =
    std::max_element(vals.begin(), vals.begin() + limit);
  const double max_val = *x;
  const size_t max_idx = x - vals.begin();
  double sum = 1.0;
  for (size_t i = 0; i < limit; ++i) {
    if (i != max_idx) {
      sum += exp(vals[i] - max_val);
#ifdef DEBUG
      assert(finite(sum));
#endif
    }
  }
  return max_val + log(sum);
}

double
SplitDistro_::log_likelihood(std::vector<double>::const_iterator a,
                             std::vector<double>::const_iterator b) const {
  double l = 0;
  for (; a < b; ++a)
    l += this->log_likelihood(*a);
  return l;
}

double
SplitDistro_::operator()(const double val) const {
  return exp(log_likelihood(val));
}

double
SplitDistro_::operator()(const vector<double> &vals) const {
  const size_t lim = vals.size();
  double l = 1;
  for (size_t i = 0; i < lim; ++i)
    l *= this->operator()(vals[i]);
  return l;
}

double
SplitDistro_::log_likelihood(const std::vector<double> &vals) const {
  double l = 0;
  const size_t lim = vals.size();
  for (size_t i = 0; i < lim; ++i)
    l += this->log_likelihood(vals[i]);
  return l;
}

double
SplitDistro_::log_likelihood(const std::vector<double> &vals,
                             const std::vector<double> &scales) const {
  double l = 0;
  const size_t lim = vals.size();
  for (size_t i = 0; i < lim; ++i)
    l += this->log_likelihood(vals[i], scales[i]);
  return l;
}

SplitDistro::SplitDistro(const std::string &n, const std::string &params) :
  name(n), d(split_distro_factory(n, params)) {}

SplitDistro::SplitDistro(const std::string &n, const std::vector<double> &params) :
  name(n), d(split_distro_factory(name)) {d->set_params(params);}

SplitDistro::SplitDistro(const std::string &s) : d(split_distro_factory(s)) {
  vector<string> name_split;
  if (s.find(",") == string::npos)  // whitespaces seperated
    name_split = smithlab::split_whitespace_quoted(s);
  else // comma seperated
    name_split = smithlab::split(s, ",");
  name = name_split[0];
}

SplitDistro::SplitDistro(const SplitDistro &rhs) : name(rhs.name),
                                                   d(split_distro_factory(rhs.name)) {
  vector<double> tmp_params(rhs.get_params());
  d->set_params(tmp_params);
}

SplitDistro&
SplitDistro::operator=(const SplitDistro &rhs) {
  if (this != &rhs) {
    name = rhs.name;
    d = split_distro_factory(rhs.name);
    vector<double> tmp_params(rhs.get_params());
    d->set_params(tmp_params);
  }
  return *this;
}

SplitDistro::~SplitDistro() {if (d) delete d;}

double
SplitDistro::operator()(double val) const {return (*d)(val);}

double
SplitDistro::operator()(const std::vector<double> &vals) const {
  return (*d)(vals);
}

void
SplitDistro::estimate_params_ml(const std::vector<double> &vals_a,
                                const std::vector<double> &vals_b) {
  d->estimate_params_ml(vals_a, vals_b);
}

void
SplitDistro::estimate_params_ml(const std::vector<double> &vals_a,
                                const std::vector<double> &vals_b,
                                const std::vector<double> &weights) {
  d->estimate_params_ml(vals_a, vals_b, weights);
}

void
SplitDistro::estimate_params_ml(const std::vector<double> &vals_a,
                                const std::vector<double> &vals_b,
                                const std::vector<double> &scales,
                                const std::vector<double> &weights) {
  d->estimate_params_ml(vals_a, vals_b, scales, weights);
}

double
SplitDistro::log_likelihood(const std::vector<double> &vals) const {
  return d->log_likelihood(vals);
}

double
SplitDistro::log_likelihood(const std::vector<double> &vals,
                            const std::vector<double> &scales) const {
  return d->log_likelihood(vals, scales);
}


double
SplitDistro::log_likelihood(std::vector<double>::const_iterator a,
                            std::vector<double>::const_iterator b) const {
  return d->log_likelihood(a, b);
}


double
SplitDistro::log_likelihood(double val) const {
  return d->log_likelihood(val);
}

double
SplitDistro::log_likelihood(const double val, const double scale) const {
  return d->log_likelihood(val, scale);
}

std::string
SplitDistro::tostring() const {
  return name + string(" ") + d->tostring();
}



std::ostream&
operator<<(std::ostream& s, const SplitDistro& distro) {
  return s << distro.tostring();
}


double
SplitDistro::log_sum_log_vec(const vector<double> &vals, size_t limit) {
  return split_distro_log_sum_log_vec(vals, limit);
}

////////////////////////////////////////////////////////////////////////
//
// DISTRO FACTORY


SplitDistro_ *
split_distro_factory(string name, string params) {
  SplitDistro_ *distro;
  if (name == "skel")
    distro = new SkellamDistro();
  else if (name == "gauss")
    distro = new GaussianDistro();
  else if (name == "nbdiff")
    distro = new NegBinomDiffDistro();
  else {
    cerr << "bad distribution name \"" << name << "\"" << endl;
    exit(EXIT_FAILURE);
  }

  vector<string> params_split = smithlab::split(params, ",");
  if (params_split.size() != distro->required_params())
    throw SMITHLABException("bad number of params: " + smithlab::toa(params_split.size()) +
                            " for split distro " + name + "\"");
  else {
    vector<double> params_vec;
    for (size_t i = 0; i < params_split.size(); ++i) {
      params_vec.push_back(atof(params_split[i].c_str()));
    }
    distro->set_params(params_vec);
  }
  return distro;
}


SplitDistro_ *
split_distro_factory(string name_arg) {

  vector<string> name_split;
  if (name_arg.find(",") == string::npos)  // whitespaces seperated
    name_split = smithlab::split_whitespace_quoted(name_arg);
  else // comma seperated
    name_split = smithlab::split(name_arg, ",");

  const string name = name_split.front();

  SplitDistro_ *distro;
  if (name == "skel")
    distro = new SkellamDistro();
  else if (name == "nbdiff")
    distro = new NegBinomDiffDistro();
  else if (name == "gauss")
    distro = new GaussianDistro();
  else
    throw SMITHLABException("bad split distribution name \"" + name + "\"");

  if (name_split.size() > 1) {
    vector<string> params_split(vector<string>(name_split.begin() + 1,
                                               name_split.end()));
    if (params_split.size() != distro->required_params())
      throw SMITHLABException("bad number of params: " + smithlab::toa(params_split.size()) +
                              " for split distro " + name + "\"");
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
//
// DISTRO

SplitDistro_::SplitDistro_() {
}

SplitDistro_::SplitDistro_(const std::vector<double> p) : params(p) {
}

SplitDistro_::~SplitDistro_() {
}


string
SplitDistro_::tostring() const {
  const std::ios_base::fmtflags split_distro_opts = std::ios_base::showpoint;
  std::ostringstream os;
  os.flags(split_distro_opts);
  if (!params.empty())
    os << std::setprecision(3) << params.front();
  for (size_t i = 1; i < params.size(); ++i)
    os << " " << std::setprecision(3) << params[i];
  return os.str();
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


double
SkellamDistro::log_likelihood(const double val) const {
  if (std::fabs(val) > 120)
    return -40;

  //   const double r = -(params.front() + params.back()) +
  //     (val/2.0)*(log(params.front()) - log(params.back())) +
  //     log(gsl_sf_bessel_In(static_cast<int>(val),
  //                     2*sqrt(params.front()*params.back())));

  const double r = -(params.front() + params.back()) +
    (val/2.0)*(log(params.front()) - log(params.back())) +
    log(gsl_sf_bessel_In(static_cast<int>(fabs(val)),
                         2*sqrt(params.front()*params.back())));
  return r;
}

double
SkellamDistro::log_likelihood(const double val,
                              const double scale) const {
  if (std::fabs(val) > 120)
    return -40;

  const double lambda1 = params[0] * scale;
  const double lambda2 = params[1] * scale;

  const double r = -(lambda1 + lambda2) +
    (val/2.0)*(log(lambda1) - log(lambda2)) +
    log(gsl_sf_bessel_In(static_cast<int>(fabs(val)),
                         2*sqrt(lambda1 * lambda2)));
  return r;
}

SkellamDistro::SkellamDistro(const SkellamDistro &rhs) :
  SplitDistro_(rhs.params) {}

SkellamDistro&
SkellamDistro::operator=(const SkellamDistro &rhs) {
  if (this != &rhs) {
    this->SplitDistro_::params = rhs.SplitDistro_::params;
  }
  return *this;
}

void
SkellamDistro::estimate_params_ml(const vector<double> &vals_a,
                                  const vector<double> &vals_b) {
  assert(vals_a.size() == vals_b.size());
  const size_t lim = vals_a.size();
  params.front() = accumulate(vals_a.begin(), vals_a.end(), 0.0)/lim;
  params.back() = accumulate(vals_b.begin(), vals_b.end(), 0.0)/lim;
}

void
SkellamDistro::estimate_params_ml(const vector<double> &vals_a,
                                  const vector<double> &vals_b,
                                  const vector<double> &probs) {
  const size_t lim = vals_a.size();
  if (workspace_vals_a.size() < lim) {
    workspace_vals_a.resize(lim);
    workspace_vals_b.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = log(probs[i]);
    workspace_vals_a[i] = log(vals_a[i]) + log(probs[i]);
    workspace_vals_b[i] = log(vals_b[i]) + log(probs[i]);
  }
  const double prob_sum = exp(split_distro_log_sum_log_vec(workspace_probs, lim));
  params.front() = exp(split_distro_log_sum_log_vec(workspace_vals_a, lim))/prob_sum;
  params.back() = exp(split_distro_log_sum_log_vec(workspace_vals_b, lim))/prob_sum;
}

void
SkellamDistro::estimate_params_ml(const vector<double> &vals_a,
                                  const vector<double> &vals_b,
                                  const vector<double> &scales,
                                  const vector<double> &probs) {

  const size_t lim = vals_a.size();
  if (workspace_vals_a.size() < lim) {
    workspace_vals_a.resize(lim);
    workspace_vals_b.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i)
    {
      workspace_probs[i] = log(probs[i]) + log(scales[i]);
      workspace_vals_a[i] = log(vals_a[i]) + log(probs[i]);
      workspace_vals_b[i] = log(vals_b[i]) + log(probs[i]);
    }
  const double prob_sum = exp(split_distro_log_sum_log_vec(workspace_probs, lim));
  params.front() = exp(split_distro_log_sum_log_vec(workspace_vals_a, lim))/prob_sum;
  params.back() = exp(split_distro_log_sum_log_vec(workspace_vals_b, lim))/prob_sum;
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

double
GaussianDistro::log_likelihood(const double val) const {
  //   if (std::fabs(val) > 120)
  //     return -40;
  //   const double r = -(params.front() + params.back()) +
  //     (val/2.0)*(log(params.front()) - log(params.back())) +
  //     log(gsl_sf_bessel_In(static_cast<int>(val),
  //                     2*sqrt(params.front()*params.back())));
  if (gsl_ran_gaussian_pdf(val - params.front(), params.back()) > 1e-30)
    return log(gsl_ran_gaussian_pdf(val - params.front(), params.back()));
  else
    return -30.0;
}

double
GaussianDistro::log_likelihood(const double val,
                               const double scale) const {

  const double mean = params[0] * scale;
  const double var = params[1] * pow(scale, 2.0);
  const double p = gsl_ran_gaussian_pdf(val - mean, var);
  return (p > 1e-30) ? p : -30.0;
}

GaussianDistro::GaussianDistro(const GaussianDistro &rhs) :
  SplitDistro_(rhs.params) {}

GaussianDistro&
GaussianDistro::operator=(const GaussianDistro &rhs) {
  if (this != &rhs) {
    this->SplitDistro_::params = rhs.SplitDistro_::params;
  }
  return *this;
}

void
GaussianDistro::estimate_params_ml(const vector<double> &vals_a,
                                   const vector<double> &vals_b) {
  assert(vals_a.size() == vals_b.size());
  const size_t lim = vals_a.size();

  vector<double> vals_c(vals_a);
  for (size_t i = 0; i < vals_c.size(); ++i)
    vals_c[i] -= vals_b[i];
  params.front() = accumulate(vals_c.begin(), vals_c.end(), 0.0)/lim;
  params.back() = gsl_stats_sd_m(&vals_c.front(), 1, vals_c.size(), params.front());
  //  cerr << params.front() << "\t" << params.back() << endl;
}

void
GaussianDistro::estimate_params_ml(const vector<double> &vals_a,
                                   const vector<double> &vals_b,
                                   const vector<double> &probs) {
  const size_t lim = vals_a.size();
  if (workspace_vals_a.size() < lim) {
    workspace_vals_a.resize(lim);
    workspace_vals_b.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = log(probs[i]);
    workspace_vals_a[i] = log(vals_a[i]) + log(probs[i]);
    workspace_vals_b[i] = log(vals_b[i]) + log(probs[i]);
  }

  vector<double> vals_c(vals_a);
  for (size_t i = 0; i < vals_c.size(); ++i)
    vals_c[i] -= vals_b[i];
  params.front() = gsl_stats_wmean(&probs.front(), 1,
                                   &vals_c.front(), 1, lim);
  params.back() = gsl_stats_wsd_m(&probs.front(), 1,
                                  &vals_c.front(), 1, lim, params.front());
  //  cerr << params.front() << "\t" << params.back() << endl;
}

void
GaussianDistro::estimate_params_ml(const vector<double> &vals_a,
                                   const vector<double> &vals_b,
                                   const vector<double> &scales,
                                   const vector<double> &probs) {

  cerr << "GaussianDistro::estimate_params_ml:  "
       << "Using an untested function" << endl;

  const size_t lim = vals_a.size();
  if (workspace_vals_a.size() < lim) {
    workspace_vals_a.resize(lim);
    workspace_vals_b.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = log(probs[i]);
    workspace_vals_a[i] = log(vals_a[i] / scales[i]) + log(probs[i]);
    workspace_vals_b[i] = log(vals_b[i] / scales[i]) + log(probs[i]);
  }

  vector<double> vals_c(vals_a);
  for (size_t i = 0; i < vals_c.size(); ++i)
    vals_c[i] = (vals_a[i] - vals_b[i]) / scales[i];
  params.front() = gsl_stats_wmean(&probs.front(), 1,
                                   &vals_c.front(), 1, lim);
  params.back() = gsl_stats_wsd_m(&probs.front(), 1,
                                  &vals_c.front(), 1, lim, params.front());
  //  cerr << params.front() << "\t" << params.back() << endl;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void
NegBinomDiffDistro::set_helpers() {
  r1 = 1/params[1];
  p1 = r1/(r1 + params[0]);

  r2 = 1/params[3];
  p2 = r2/(r2 + params[2]);

  q1 = 1 - p1;
  q2 = 1 - p2;
  c = log(q1) + log(q2);

  r1_log_p1_M_lnG_r1_P_r2_log_p2_M_lnG_r2 =
    (r1*log(p1) - gsl_sf_lngamma(r1) + r2*log(p2) - gsl_sf_lngamma(r2));
  //TODO: should check that these are valid!!!
}

const double NegBinomDiffDistro::tolerance = 1e-160;

static inline double
log_sum_log(const double p, const double q) {
  // Never pass a '0' into this function!!!!
  return ((p > q) ? p + log(1.0 + exp(q - p)) : q + log(1.0 + exp(p - q)));
}

/// (Andrew) Function below not used so commented out
/*
inline static double
log_nbd_pdf(const unsigned int k, const double p, double n) {
  const double f = gsl_sf_lngamma(k + n);
  const double a = gsl_sf_lngamma(n);
  const double b = gsl_sf_lngamma(k + 1.0);
  return (f-a-b) + n*log(p) + k*log(1 - p);
}
*/

double
NegBinomDiffDistro::log_likelihood(const double k) const {
  const int k_int = static_cast<int>(k);
  if (ll_hash.find(k_int) == ll_hash.end()) {

    double sum = 0;
    double prev = -1;
    double n = std::max(0.0, -k);

    while (fabs((prev - sum)/sum) > tolerance) {
      prev = sum;
      const double log_form =
        (gsl_sf_lngamma(k + n + r1) -
         gsl_sf_lnfact(static_cast<size_t>(k + n)) +
         gsl_sf_lngamma(n + r2) -
         gsl_sf_lnfact(static_cast<size_t>(n))) + n*c;
      sum = ((sum == 0) ? log_form : log_sum_log(sum, log_form));
      ++n;
    }
    sum += k*log(q1) + r1_log_p1_M_lnG_r1_P_r2_log_p2_M_lnG_r2;
    ll_hash[k_int] = sum;
  }
  return ll_hash[k_int];
}

double
NegBinomDiffDistro::log_likelihood(const double k,
                                   const double scale) const {

  const double scaled_p1 = r1 / (r1 + params[0] * scale);
  const double scaled_p2 = r2 / (r2 + params[2] * scale);
  const double scaled_log_q1_log_q2 = log(1 - scaled_p1)
    + log(1 - scaled_p2);
  const double scaled_r1_log_p1_M_lnG_r1_P_r2_log_p2_M_lnG_r2 =
    (r1 * log(scaled_p1) - gsl_sf_lngamma(r1)
     + r2 * log(scaled_p2) - gsl_sf_lngamma(r2));

  double sum = 0;
  double prev = -1;
  double n = std::max(0.0, -k);

  while (fabs((prev - sum)/sum) > tolerance)
    {
      prev = sum;
      const double log_form =
        (gsl_sf_lngamma(k + n + r1) -
         gsl_sf_lnfact(static_cast<size_t>(k + n)) +
         gsl_sf_lngamma(n + r2) -
         gsl_sf_lnfact(static_cast<size_t>(n))) + n * scaled_log_q1_log_q2;
      sum = ((sum == 0) ? log_form : log_sum_log(sum, log_form));
      ++n;
    }

  sum += k * log(1 - scaled_p1) + scaled_r1_log_p1_M_lnG_r1_P_r2_log_p2_M_lnG_r2;
  return sum;
}


NegBinomDiffDistro::NegBinomDiffDistro(const NegBinomDiffDistro &rhs) :
  SplitDistro_(rhs.params) {
  set_helpers();
}

NegBinomDiffDistro&
NegBinomDiffDistro::operator=(const NegBinomDiffDistro &rhs) {
  if (this != &rhs) {
    this->SplitDistro_::params = rhs.SplitDistro_::params;
    set_helpers();
  }
  return *this;
}


static inline double
score_first_term(const vector<double> &v_hist, const double mu, const double alpha) {
  double sum = 0;
  const size_t lim = v_hist.size();
  for (size_t i = 0; i < lim; ++i)
    if (v_hist[i] > 0) {
      double inner_sum = 0;
      for (size_t j = 0; j < i; ++j)
        inner_sum += j/(1.0 + alpha*j);
      sum += v_hist[i]*inner_sum;
    }
  return sum;
}

static inline double
alpha_score_function(const vector<double> &vals_hist, const double mu,
                     const double alpha, const double vals_count) {
  const double one_plus_alpha_mu = 1 + alpha*mu;
  return (score_first_term(vals_hist, mu, alpha)/vals_count +
          (log(one_plus_alpha_mu)/alpha - mu)/alpha);
}

static const double max_allowed_alpha = 10;
static const double min_allowed_alpha = 1e-5;
static const double alpha_allowed_error = 1e-10;

static void
estimate_params_pair_ml(const vector<double> &vals, double &mu, double &alpha) {
  const size_t lim = vals.size();
  // This is the mu
  mu = std::accumulate(vals.begin(), vals.end(), 0.0)/lim;

  // Now for the alpha
  const double max_value = *std::max_element(vals.begin(), vals.end());
  vector<double> vals_hist(static_cast<size_t>(max_value) + 1, 0.0);
  for (size_t i = 0; i < lim; ++i)
    ++vals_hist[static_cast<size_t>(vals[i])];

  const double vals_count = lim;

  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  double a_mid = max_allowed_alpha;

  double diff = std::numeric_limits<double>::max();
  double prev_val = std::numeric_limits<double>::max();
  while (diff > alpha_allowed_error) {
    a_mid = (a_low + a_high)/2;
    const double mid_val = alpha_score_function(vals_hist, mu, a_mid, vals_count);
    if (mid_val < 0)
      a_high = a_mid;
    else
      a_low = a_mid;
    diff = std::fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  alpha = a_mid;
  // METHOD OF MOMENTS:
  //   const size_t lim = vals.size();
  //   mu = std::accumulate(vals.begin(), vals.begin() + lim, 0.0)/lim;
  //   const double var = gsl_stats_variance_m(&vals.front(), 1, vals.size(), mu);
  //   const double r = (mu*mu)/(var - mu);
  //   alpha = max(min_alpha, 1/r);
}

static void
estimate_params_pair_ml(const vector<double> &vals, const vector<double> &probs,
                        vector<double> &workspace_vals, vector<double> &workspace_probs,
                        double &mu, double &alpha) {
  const size_t lim = vals.size();
  if (workspace_vals.size() < lim) {
    workspace_vals.resize(lim);
    workspace_probs.resize(lim);
  }

  for (size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = log(probs[i]);
    workspace_vals[i] = log(vals[i]) + log(probs[i]);
  }

  const double vals_count = exp(SplitDistro::log_sum_log_vec(workspace_probs, lim));
  mu = exp(SplitDistro::log_sum_log_vec(workspace_vals, lim))/vals_count;

  const double max_value = *std::max_element(vals.begin(), vals.end());
  vector<double> vals_hist(static_cast<size_t>(max_value) + 1, 0.0);
  for (size_t i = 0; i < lim; ++i)
    vals_hist[static_cast<size_t>(vals[i])] += probs[i];

  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  double a_mid = max_allowed_alpha;

  double diff = std::numeric_limits<double>::max();
  double prev_val = std::numeric_limits<double>::max();
  while (diff > alpha_allowed_error) {
    a_mid = (a_low + a_high)/2;
    const double mid_val = alpha_score_function(vals_hist, mu, a_mid, vals_count);
    if (mid_val < 0)
      a_high = a_mid;
    else
      a_low = a_mid;
    diff = std::fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  alpha = a_mid;
  // METHOD OF MOMENTS:
  //   mu = exp(split_distro_log_sum_log_vec(workspace_vals, lim))/vals_count;
  //   const double var = gsl_stats_wvariance_m(&vals.front(), 1, &probs.front(), 1,
  //                                       vals.size(), mu);
  //   const double r = (mu*mu)/(var - mu);
  //   alpha = max(min_alpha, 1/r);
}

void
NegBinomDiffDistro::estimate_params_ml(const vector<double> &vals_a,
                                       const vector<double> &vals_b) {
  hq_estimate_params_ml(vals_a, vals_b);
}

void
NegBinomDiffDistro::estimate_params_ml(const vector<double> &vals_a,
                                       const vector<double> &vals_b,
                                       const vector<double> &probs) {
  hq_estimate_params_ml(vals_a, vals_b, probs);
}


////// andrew's orignial version
void
NegBinomDiffDistro::andrew_estimate_params_ml(const vector<double> &vals_a,
                                              const vector<double> &vals_b) {
  estimate_params_pair_ml(vals_a, params[0], params[1]);
  estimate_params_pair_ml(vals_b, params[2], params[3]);

  ll_hash.clear();
  set_helpers();
}


void
NegBinomDiffDistro::andrew_estimate_params_ml(const vector<double> &vals_a,
                                              const vector<double> &vals_b,
                                              const vector<double> &probs) {
  estimate_params_pair_ml(vals_a, probs, workspace_vals_a, workspace_probs,
                          params[0], params[1]);
  estimate_params_pair_ml(vals_b, probs, workspace_vals_b, workspace_probs,
                          params[2], params[3]);
  ll_hash.clear();
  set_helpers();
}
///////

/////// experimental code for gradient descent method
inline void
my_dummy_func(const double &r)
{
  // used to make sure a tempory variable are not optimized out by compiler
  return;
}

double
local_log_likelihood(const vector<double> &vals, const vector<double> &probs,
                     const double &mu1, const double &alpha1,
                     const double &mu2, const double &alpha2)
{
  static const double tolerance = 1e-160;
  const double r1 = 1 / alpha1;
  const double p1 = r1 / (r1 + mu1);
  const double log_q1 = log(1 - p1);

  const double r2 = 1 / alpha2;
  const double p2 = r2 / (r2 + mu2);

  const double fixed_term = (r1*log(p1) - gsl_sf_lngamma(r1)
                             + r2*log(p2) - gsl_sf_lngamma(r2));

  const double log_q1_log_q2 = log(1 - p1) + log(1 - p2);

  unordered_map<int, double> ll_hash;

  double llh = 0;
  for (size_t i = 0; i < vals.size(); ++i)
    {
      const int val_int = static_cast<int>(vals[i]);
      if (ll_hash.find(val_int) == ll_hash.end())
        {
          double sum = 0;
          double prev = -1;
          double n = std::max(0.0, - vals[i]);

          while (fabs((prev - sum)/sum) > tolerance) {
            prev = sum;
            const double log_form =
              (gsl_sf_lngamma(val_int + n + r1) -
               gsl_sf_lnfact(static_cast<size_t>(val_int + n)) +
               gsl_sf_lngamma(n + r2) -
               gsl_sf_lnfact(static_cast<size_t>(n))) + n * log_q1_log_q2;
            sum = ((sum == 0) ? log_form : log_sum_log(sum, log_form));
            ++n;
          }
          sum += val_int * log_q1 + fixed_term;
          ll_hash[val_int] = sum;
        }
      llh += probs[i] * ll_hash[val_int];
    }

  return llh;
}

double
local_log_likelihood(const vector<double> &vals,
                     const double &mu1, const double &alpha1,
                     const double &mu2, const double &alpha2)
{
  const vector<double> probs(vals.size(), 1);
  return local_log_likelihood(vals, probs, mu1, alpha1,
                              mu2, alpha2);
}


void
NegBinomDiffDistro::hq_estimate_params_ml(const vector<double> &vals_a,
                                          const vector<double> &vals_b)
{
  const vector<double> probs(vals_a.size(), 1);
  hq_estimate_params_ml(vals_a, vals_b, probs);
}

void
NegBinomDiffDistro::hq_estimate_params_ml(const std::vector<double> &vals_a,
                                          const std::vector<double> &vals_b,
                                          const std::vector<double> &probs)
{

  //  static const double max_allow_error = 1e-10;
  static const double max_allow_llh_diff = 1e-12;
  //    static const double max_delta_norm = 1e-20;
  static const double sqrt_error_machine = 1e-10;

  double global_scale = 1.0 / 16;

  vector<double> vals(vals_a.size());
  for (size_t i = 0; i < vals.size(); ++i)
    vals[i] = vals_a[i] - vals_b[i];

  double mu1, alpha1, mu2, alpha2;

  // get initial value
  estimate_params_pair_ml(vals_a, probs,
                          workspace_vals_a, workspace_probs,
                          mu1, alpha1);
  estimate_params_pair_ml(vals_b, probs,
                          workspace_vals_b, workspace_probs,
                          mu2, alpha2);

  // newton gradient
  double temp(0);
  double h(0);

  double mu1_grad(0), alpha1_grad(0), mu2_grad(0), alpha2_grad(0);
  double prev_llh = local_log_likelihood(vals, probs, mu1, alpha1, mu2, alpha2);
  double llh_diff = std::numeric_limits<double>::max();
  double curr_llh = 0;

  size_t iter(0);
  while (fabs(llh_diff / curr_llh) > max_allow_llh_diff)
    {
      ++iter;
      mu1 += mu1_grad;
      alpha1 += alpha1_grad;
      mu2 += mu2_grad;
      alpha2 += alpha2_grad;

      // compute gradient
      h = sqrt_error_machine * mu1; // stepwise to be determined wisely
      temp = h + mu1;
      my_dummy_func(temp);
      h = temp - mu1;
      mu1_grad  = (local_log_likelihood(vals, probs, mu1 + h, alpha1, mu2, alpha2)
                   -  local_log_likelihood(vals, probs, mu1 - h, alpha1, mu2, alpha2))
        / h / 2;

      h = sqrt_error_machine * alpha1; // stepwise to be determined wisely
      temp = h + alpha1;
      my_dummy_func(temp);
      h = temp - alpha1;
      alpha1_grad  = (local_log_likelihood(vals, probs, mu1, alpha1 + h, mu2, alpha2)
                      - local_log_likelihood(vals, probs, mu1, alpha1 - h, mu2, alpha2))
        / h / 2;

      h = sqrt_error_machine * mu2; // stepwise to be determined wisely
      temp = h + mu2;
      my_dummy_func(temp);
      h = temp - mu2;
      mu2_grad  = (local_log_likelihood(vals, probs, mu1, alpha1, mu2 + h, alpha2)
                   - local_log_likelihood(vals, probs, mu1, alpha1, mu2 - h, alpha2))
        / h / 2;

      h = sqrt_error_machine * alpha2; // stepwise to be determined wisely
      temp = h + alpha2;
      my_dummy_func(temp);
      h = temp - alpha2;
      alpha2_grad = (local_log_likelihood(vals, probs, mu1, alpha1, mu2, alpha2 + h)
                     - local_log_likelihood(vals, probs, mu1, alpha1, mu2, alpha2 - h))
        / h / 2;

      double local_scale = min(min(mu1 / fabs(mu1_grad + 1e-20),
                                   alpha1 / fabs(alpha1_grad + 1e-20)),
                               min(mu2 / fabs(mu2_grad + 1e-20),
                                   alpha2 / fabs(alpha2_grad + 1e-20)));
      mu1_grad *= local_scale * global_scale;
      alpha1_grad *= local_scale * global_scale;
      mu2_grad *= local_scale * global_scale;
      alpha2_grad *= local_scale * global_scale;

      curr_llh = local_log_likelihood(vals, probs,
                                      mu1 + mu1_grad, alpha1 + alpha1_grad,
                                      mu2 + mu2_grad, alpha2 + alpha2_grad);
      llh_diff = curr_llh - prev_llh;

      while (llh_diff < - 1e-100) // if the move can not lead to greater likelihood
        // cancel this remove and try a smaller step
        {
          mu1_grad /= 2;
          alpha1_grad /= 2;
          mu2_grad /= 2;
          alpha2_grad /= 2;

          curr_llh = local_log_likelihood(vals, probs,
                                          mu1 + mu1_grad, alpha1 + alpha1_grad,
                                          mu2 + mu2_grad, alpha2 + alpha2_grad);
          llh_diff = curr_llh - prev_llh;
        }
      prev_llh = curr_llh;
    }

  params[0] = mu1;
  params[1] = alpha1;
  params[2] = mu2;
  params[3] = alpha2;

  ll_hash.clear();
  set_helpers();
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
estimate_params_pair_ml(const std::vector<double> &vals,
                        const std::vector<double> &scales,
                        const std::vector<double> &probs,
                        vector<double> &workspace_vals,
                        vector<double> &workspace_probs,
                        double &mu, double &alpha)
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
  mu = exp(split_distro_log_sum_log_vec(workspace_vals, lim)) /
    exp(split_distro_log_sum_log_vec(workspace_probs, lim));

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
  alpha = a_mid;
}

double
local_log_likelihood(const vector<double> &vals,
                     const vector<double> &scales,
                     const vector<double> &probs,
                     const double &mu1, const double &alpha1,
                     const double &mu2, const double &alpha2)
{
  static const double tolerance = 1e-160;

  double llh = 0;
  for (size_t i = 0; i < vals.size(); ++i)
    {
      const double r1 = 1 / alpha1;
      const double p1 = r1 / (r1 + mu1 * scales[i]);
      const double log_q1 = log(1 - p1);

      const double r2 = 1 / alpha2;
      const double p2 = r2 / (r2 + mu2 * scales[i]);

      const double fixed_term = (r1*log(p1) - gsl_sf_lngamma(r1)
                                 + r2*log(p2) - gsl_sf_lngamma(r2));

      const double log_q1_log_q2 = log(1 - p1) + log(1 - p2);

      const int val_int = static_cast<int>(vals[i]);

      double sum = 0;
      double prev = -1;
      double n = std::max(0.0, - vals[i]);

      while (fabs((prev - sum)/sum) > tolerance)
        {
          prev = sum;
          const double log_form =
            (gsl_sf_lngamma(val_int + n + r1) -
             gsl_sf_lnfact(static_cast<size_t>(val_int + n)) +
             gsl_sf_lngamma(n + r2) -
             gsl_sf_lnfact(static_cast<size_t>(n))) + n * log_q1_log_q2;
          sum = ((sum == 0) ? log_form : log_sum_log(sum, log_form));
          ++n;
        }
      sum += val_int * log_q1 + fixed_term;

      llh += probs[i] * sum;
    }

  return llh;
}


void
NegBinomDiffDistro::estimate_params_ml(const std::vector<double> &vals_a,
                                       const std::vector<double> &vals_b,
                                       const std::vector<double> &scales,
                                       const std::vector<double> &probs) {
  estimate_params_pair_ml(vals_a, scales, probs,
                          workspace_vals_a, workspace_probs,
                          params[0], params[1]);
  estimate_params_pair_ml(vals_b, scales, probs,
                          workspace_vals_b, workspace_probs,
                          params[2], params[3]);
}
