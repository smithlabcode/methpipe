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

#ifndef SPLIT_DISTRO_HPP
#define SPLIT_DISTRO_HPP

#include <algorithm>
#include <vector>
#include <iostream>
#include <string>

#include <tr1/unordered_map>

class SplitDistro_ {
public:
  SplitDistro_();
  SplitDistro_(const std::vector<double> p);
  SplitDistro_(const SplitDistro_ &);
  SplitDistro_& operator=(const SplitDistro_ &);
  virtual ~SplitDistro_();

  virtual size_t required_params() const = 0;

  virtual double get_mean() const = 0;
  double operator()(double val) const;
  double operator()(const std::vector<double> &) const;

  virtual double log_likelihood(double val) const = 0;
  virtual double log_likelihood(const double val, const double scale) const = 0;
  double log_likelihood(const std::vector<double> &vals) const;
  double log_likelihood(const std::vector<double> &vals,
			const std::vector<double> &scales) const;
  double log_likelihood(std::vector<double>::const_iterator a,
			std::vector<double>::const_iterator b) const;
  virtual void
  estimate_params_ml(const std::vector<double> &,
		     const std::vector<double> &) = 0;
  virtual void
  estimate_params_ml(const std::vector<double> &,
		     const std::vector<double> &,
		     const std::vector<double> &) = 0;
  virtual void
  estimate_params_ml(const std::vector<double> &vals_a,
                     const std::vector<double> &vals_b,
                     const std::vector<double> &scales,
                     const std::vector<double> &probs) = 0;

  std::string tostring() const;
  std::vector<double> get_params() const {return params;}

  virtual void
  set_params(const std::vector<double> &p) {params = p;}
  static double
  log_sum_log_vec(const std::vector<double> &vals, size_t limit);

protected:

  std::vector<double> params;
  std::vector<double> workspace_vals_a;
  std::vector<double> workspace_vals_b;
  std::vector<double> workspace_probs;
};

SplitDistro_ *
split_distro_factory(std::string name, std::string params);

SplitDistro_ *
split_distro_factory(std::string name);

class SplitDistro {
public:

  SplitDistro() : d(0) {}
  SplitDistro(const std::string &name, const std::string &params);
  SplitDistro(const std::string &name, const std::vector<double> &params);
  SplitDistro(const std::string &s);
  SplitDistro(const SplitDistro &);
  SplitDistro& operator=(const SplitDistro &);
  ~SplitDistro();

  double get_mean() const {return d->get_mean();}


  double operator()(double val) const;
  double operator()(const std::vector<double> &) const;
  void estimate_params_ml(const std::vector<double> &,
			  const std::vector<double> &);
  void estimate_params_ml(const std::vector<double> &,
			  const std::vector<double> &,
			  const std::vector<double> &);
  void estimate_params_ml(const std::vector<double> &vals_a,
			  const std::vector<double> &vals_b,
			  const std::vector<double> &scales,
			  const std::vector<double> &probs);

  double
  log_likelihood(const std::vector<double> &) const;
  double
  log_likelihood(const std::vector<double> &vals,
                 const std::vector<double> &scales) const;
  double
  log_likelihood(std::vector<double>::const_iterator a,
		 std::vector<double>::const_iterator b) const;
  double log_likelihood(double val) const;
  double log_likelihood(const double val, const double scale) const;
  std::string tostring() const;
  std::vector<double> get_params() const {return d->get_params();}
  void set_params(const std::vector<double> &p) {d->set_params(p);}

  static double
  log_sum_log_vec(const std::vector<double> &vals, size_t limit);


private:
  std::string name;
  SplitDistro_ *d;
};

std::ostream&
operator<<(std::ostream& s, const SplitDistro& distro);

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// SKELLAM DISTRIBUTION

class SkellamDistro : public SplitDistro_ {
public:
  SkellamDistro() : SplitDistro_(std::vector<double>(2, 0)) {}
  SkellamDistro(std::vector<double> p) : SplitDistro_(p) {}
  SkellamDistro(const SkellamDistro &rhs);
  SkellamDistro& operator=(const SkellamDistro &rhs);
  ~SkellamDistro() {}

  double get_mean() const {return params[0] - params[1];}
  size_t required_params() const {return 2;}
  double log_likelihood(double val) const;
  double log_likelihood(const double val, const double scale) const;
  void estimate_params_ml(const std::vector<double> &vals_a,
			  const std::vector<double> &vals_b);
  void estimate_params_ml(const std::vector<double> &vals_a,
			  const std::vector<double> &vals_b,
			  const std::vector<double> &probs);
  void estimate_params_ml(const std::vector<double> &vals_a,
			  const std::vector<double> &vals_b,
			  const std::vector<double> &scales,
			  const std::vector<double> &probs);
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// GAUSSIAN DISTRIBUTION

class GaussianDistro : public SplitDistro_ {
public:
  GaussianDistro() : SplitDistro_(std::vector<double>(2, 0)) {}
  GaussianDistro(std::vector<double> p) : SplitDistro_(p) {}
  GaussianDistro(const GaussianDistro &rhs);
  GaussianDistro& operator=(const GaussianDistro &rhs);
  ~GaussianDistro() {}

  double get_mean() const {return params[0];}
  size_t required_params() const {return 2;}
  double log_likelihood(double val) const;
  double log_likelihood(const double val, const double scale) const;
  void estimate_params_ml(const std::vector<double> &vals_a,
			  const std::vector<double> &vals_b);
  void estimate_params_ml(const std::vector<double> &vals_a,
			  const std::vector<double> &vals_b,
			  const std::vector<double> &probs);
  void estimate_params_ml(const std::vector<double> &vals_a,
			  const std::vector<double> &vals_b,
			  const std::vector<double> &scales,
			  const std::vector<double> &probs);
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// NEGATIVE BINOMIAL DIFFERENCE DISTRIBUTION

class NegBinomDiffDistro : public SplitDistro_ {
public:
  NegBinomDiffDistro() : SplitDistro_(std::vector<double>(4, 0)) {}
  NegBinomDiffDistro(std::vector<double> p) : SplitDistro_(p) {}
  NegBinomDiffDistro(const NegBinomDiffDistro &rhs);
  NegBinomDiffDistro& operator=(const NegBinomDiffDistro &rhs);
  ~NegBinomDiffDistro() {}

  double get_mean() const {return params[0] - params[2];}

  size_t required_params() const {return 4;}
  void hq_estimate_params_ml(const std::vector<double> &vals_a,
			     const std::vector<double> &vals_b);
  void andrew_estimate_params_ml(const std::vector<double> &vals_a,
				 const std::vector<double> &vals_b);
  void estimate_params_ml(const std::vector<double> &vals_a,
                          const std::vector<double> &vals_b);
  void hq_estimate_params_ml(const std::vector<double> &vals_a,
			     const std::vector<double> &vals_b,
			     const std::vector<double> &probs);
  void andrew_estimate_params_ml(const std::vector<double> &vals_a,
				 const std::vector<double> &vals_b,
				 const std::vector<double> &probs);
  void estimate_params_ml(const std::vector<double> &vals_a,
			  const std::vector<double> &vals_b,
			  const std::vector<double> &probs);
  void estimate_params_ml(const std::vector<double> &vals_a,
			  const std::vector<double> &vals_b,
			  const std::vector<double> &scales,
			  const std::vector<double> &probs);

  double log_likelihood(double val) const;
  double log_likelihood(const double val, const double scale) const;

  void set_params(const std::vector<double> &p) {
    SplitDistro_::set_params(p);
    ll_hash.clear();
    set_helpers();
  }
private:
  void set_helpers();
  static const double tolerance;
  static const double max_allowed_alpha;
  static const double min_allowed_alpha;
  static const double alpha_allowed_error;

  mutable std::tr1::unordered_map<int, double> ll_hash;

  double r1;
  double p1;

  double r2;
  double p2;

  double q1;
  double q2;
  double c;

  double r1_log_p1_M_lnG_r1_P_r2_log_p2_M_lnG_r2;
};

#endif
