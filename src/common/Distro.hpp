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

#ifndef DISTRO_HPP
#define DISTRO_HPP

#include <algorithm>
#include <vector>
#include <iostream>
#include <string>

#include <gsl/gsl_rng.h>

using std::vector;

class Distro_ {
public:
  
  Distro_();
  Distro_(const std::vector<double> p);
  Distro_(const Distro_ &);
  Distro_& operator=(const Distro_ &);
  virtual ~Distro_();
  
  void seed(int s);

  virtual size_t required_params() const = 0;
  virtual double sample() const = 0;
  virtual void estimate_params_ml(const std::vector<double> &) = 0;
  virtual void estimate_params_ml(const std::vector<double> &vals,
                                  const std::vector<double> &scales,
                                  const std::vector<double> &probs) = 0;
  virtual void estimate_params_ml(const std::vector<double> &,
				  const std::vector<double> &) = 0;
  virtual double log_likelihood(double val) const = 0;
  virtual double log_likelihood(const double &val, const double &scale) const = 0;
  double operator()(double val) const;
  double operator()(const std::vector<double> &) const;
  double log_likelihood(const std::vector<double> &vals) const;
    double log_likelihood(const std::vector<double> &vals,
                          const std::vector<double> &scales) const;
  double log_likelihood(std::vector<double>::const_iterator a,
			std::vector<double>::const_iterator b) const;
  
  std::string tostring() const;

  std::vector<double> get_params() const {return params;}

   virtual void set_params(const std::vector<double> &p) = 0;

  static double
  log_sum_log_vec(const std::vector<double> &vals, size_t limit);
  
protected:
  
  std::vector<double> params;
  std::vector<double> workspace_vals;
  std::vector<double> workspace_probs;

  gsl_rng *rng;
};

Distro_ *
distro_factory(std::string name, std::string params);

Distro_ *
distro_factory(std::string name);

class Distro {
public:
  
  Distro() : d(0) {}
  Distro(const std::string &name, const std::string &params);
  Distro(const std::string &name, const std::vector<double> &params);
  Distro(const std::string &s);
  Distro(const Distro &);
  Distro& operator=(const Distro &);
  ~Distro();
  
  void seed(int s) {d->seed(s);}
  
  double operator()() const {return d->sample();}
  double operator()(double val) const;
  double operator()(const std::vector<double> &) const;
  void estimate_params_ml(const std::vector<double> &);
  void estimate_params_ml(const std::vector<double> &vals,
                          const std::vector<double> &scales,
                          const std::vector<double> &probs);
  void estimate_params_ml(const std::vector<double> &,
			  const std::vector<double> &);
  double log_likelihood(const std::vector<double> &) const;
    double log_likelihood(const std::vector<double> &vals,
                          const std::vector<double> &scales) const;
    double log_likelihood(std::vector<double>::const_iterator a,
			std::vector<double>::const_iterator b) const;
  double log_likelihood(double val) const;
    double log_likelihood(const double &val, const double &scale) const;
  std::string tostring() const;
  std::vector<double> get_params() const {return d->get_params();}

  void set_params(const std::vector<double> &p) {d->set_params(p);}
  double sample() const {return d->sample();}
  
  static double 
  log_sum_log_vec(const std::vector<double> &vals, size_t limit);
  static bool has_params(const std::string &name);

private:
  std::string name;
  Distro_ *d;
};

std::ostream&
operator<<(std::ostream& s, const Distro& distro);

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// EXPONENTIAL

class ExpDistro : public Distro_ {
public:
  ExpDistro() : Distro_(std::vector<double>(1, 0)) {}
  ExpDistro(std::vector<double> p) : Distro_(p) {}
  ExpDistro(const ExpDistro &rhs);
  ExpDistro& operator=(const ExpDistro &rhs);
  ~ExpDistro() {}
  double sample() const;
  void set_params(const std::vector<double> &p);
  size_t required_params() const {return 1;}
  void estimate_params_ml(const std::vector<double> &vals);
    void estimate_params_ml(const std::vector<double> &vals,
                          const std::vector<double> &scales,
                          const std::vector<double> &probs);
  void estimate_params_ml(const std::vector<double> &vals,
			  const std::vector<double> &probs);
  double log_likelihood(double val) const;
    double log_likelihood(const double &val, const double &scale) const;
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// GAMMA

class Gamma : public Distro_ {
public:
  Gamma() : Distro_(std::vector<double>(2, 1.0)) {}
  Gamma(std::vector<double> p) : Distro_(p) {}
  Gamma(const Gamma &rhs);
  Gamma& operator=(const Gamma &rhs);
  ~Gamma() {}
  double sample() const;
  void set_params(const std::vector<double> &p);
  size_t required_params() const {return 2;}
  void estimate_params_ml(const std::vector<double> &vals);
    void estimate_params_ml(const std::vector<double> &vals,
                          const std::vector<double> &scales,
                          const std::vector<double> &probs);
  void estimate_params_ml(const std::vector<double> &vals,
			  const std::vector<double> &probs);
  double log_likelihood(double val) const;
    double log_likelihood(const double &val, const double &scale) const;

private:
  void check_params_and_set_helpers();
  double gsl_sf_lngamma_k_plus_k_log_theta;
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// POISSON


class PoisDistro : public Distro_ {
public:
  PoisDistro() : Distro_(std::vector<double>(1, 0)) {}
  PoisDistro(std::vector<double> p) : Distro_(p) {}
  PoisDistro(const PoisDistro &rhs);
  PoisDistro& operator=(const PoisDistro &rhs);
  ~PoisDistro() {}
  double sample() const;
  size_t required_params() const {return 1;}
    void set_params(const std::vector<double> &p);
  void estimate_params_ml(const std::vector<double> &vals);
  void estimate_params_ml(const std::vector<double> &vals,
			  const std::vector<double> &probs);
    void estimate_params_ml(const std::vector<double> &vals,
                          const std::vector<double> &scales,
                          const std::vector<double> &probs);
  double log_likelihood(double val) const;
    double log_likelihood(const double &val, const double &scale) const;
};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// GEOMETRIC

class GeoDistro : public Distro_ {
public:
  GeoDistro() : Distro_(std::vector<double>(1, 0)) {}
  GeoDistro(std::vector<double> p) : Distro_(p) {}
  GeoDistro(const GeoDistro &rhs);
  GeoDistro& operator=(const GeoDistro &rhs);
  ~GeoDistro() {}
  double sample() const;
  size_t required_params() const {return 1;}
    void set_params(const std::vector<double> &p);
  void estimate_params_ml(const std::vector<double> &vals);
  void estimate_params_ml(const std::vector<double> &vals,
			  const std::vector<double> &probs);
    void estimate_params_ml(const std::vector<double> &vals,
                          const std::vector<double> &scales,
                          const std::vector<double> &probs);
  double log_likelihood(double val) const;
    double log_likelihood(const double &val, const double &scale) const;
};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// NEGATIVE BINOMIAL

class NegBinomDistro : public Distro_ {
public:

  NegBinomDistro() : Distro_(std::vector<double>(2, 0)) {}
  NegBinomDistro(std::vector<double> p) : Distro_(p) {set_helpers();}
  NegBinomDistro(const NegBinomDistro &rhs);
  NegBinomDistro& operator=(const NegBinomDistro &rhs);
  ~NegBinomDistro() {}
  void set_helpers();
  double sample() const;
  size_t required_params() const {return 2;}
  void set_params(const std::vector<double> &p);
  void estimate_params_ml(const std::vector<double> &vals);
  void estimate_params_ml(const std::vector<double> &vals,
			  const std::vector<double> &probs);
  void estimate_params_ml(const std::vector<double> &vals,
                          const std::vector<double> &scales,
                          const std::vector<double> &probs);

  double log_likelihood(double val) const;
    double log_likelihood(const double &val, const double &scale) const;
  
private:

  static const double max_allowed_alpha;
  static const double min_allowed_alpha;
  static const double alpha_allowed_error;

  double n_helper;
  double p_helper;
  double n_log_p_minus_lngamma_n_helper;
  double log_q_helper;
};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Beta

class Beta : public Distro_ {
public:
  Beta() : Distro_(std::vector<double>(1, 1)) {}
  Beta(std::vector<double> p) : Distro_(p) {}
  Beta(const Beta &rhs);
  Beta& operator=(const Beta &rhs);
  ~Beta() {}
  double sample() const;
  size_t required_params() const {return 2;}
    void set_params(const std::vector<double> &p);
  void estimate_params_ml(const std::vector<double> &vals);
  void estimate_params_ml(const std::vector<double> &vals,
			  const std::vector<double> &weights);
    void estimate_params_ml(const std::vector<double> &lp_vals,
                          const std::vector<double> &lq_vals,
                          const std::vector<double> &weights);
  double log_likelihood(double val) const;
    double log_likelihood(const double &val, const double &scale) const;
private:
  double alpha, beta, lnbeta_helper;
};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Binomial

class Binom : public Distro_ {
public:
  Binom() : Distro_(std::vector<double>(1, 0.5)) {}
  Binom(std::vector<double> p) : Distro_(p) {}
  Binom(const Binom &rhs);
  Binom& operator=(const Binom &rhs);
  ~Binom() {}
  double sample() const;
  size_t required_params() const {return 2;}
    void set_params(const std::vector<double> &p);
  void estimate_params_ml(const std::vector<double> &vals);
  void estimate_params_ml(const std::vector<double> &vals,
			  const std::vector<double> &weights);
    void estimate_params_ml(const std::vector<double> &lp_vals,
                          const std::vector<double> &lq_vals,
                          const std::vector<double> &weights);
  double log_likelihood(double val) const;
    double log_likelihood(const double &val, const double &scale) const;
private:
  double p;
  int n;
};




////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// EMPIRICAL

class EmpiricalDistro : public Distro_ {
public:

  EmpiricalDistro() : Distro_(std::vector<double>(2, 0)) {}
  EmpiricalDistro(std::vector<double> p) : Distro_(p) {}
  EmpiricalDistro(const EmpiricalDistro &rhs);
  EmpiricalDistro& operator=(const EmpiricalDistro &rhs);
  ~EmpiricalDistro() {}
  
  void set_params(const std::vector<double> &p) {params = p;}
  double sample() const;
  size_t required_params() const {return 2;}
  void estimate_params_ml(const std::vector<double> &vals);
  void estimate_params_ml(const std::vector<double> &vals,
			  const std::vector<double> &probs);
  void estimate_params_ml(const std::vector<double> &vals,
                          const std::vector<double> &scales,
                          const std::vector<double> &probs);
  double log_likelihood(double val) const;
    double log_likelihood(const double &val, const double &scale) const;


private:

  static void get_breaks(std::vector<double> data,
			 size_t n_vals,
			 size_t n_class,
			 double max_val,
			 std::vector<double> &breaks);
  static void get_breaks(const std::vector<double> &data,
			 const std::vector<double> &weights,
			 size_t n_vals,
			 size_t n_class,
			 double max_val,
			 std::vector<double> &breaks);

  static void make_hist(const std::vector<double> &data,
			size_t n_vals,
			size_t n_class,
			double max_val,
			std::vector<double> &breaks,
			std::vector<double> &hist);
  
  static void make_weighted_hist(const std::vector<double> &data,
				 const std::vector<double> &weights,
				 size_t n_vals,
				 size_t n_class,
				 double max_val,
				 std::vector<double> &breaks,
				 std::vector<double> &hist);

  static double
  estimate_bandwidth(const std::vector<double> &vals);
  static double
  estimate_number_of_classes(const std::vector<double> &vals);
  
  static double
  estimate_bandwidth(const std::vector<double> &vals,
		     const std::vector<double> &probs);
  static double
  estimate_number_of_classes(const std::vector<double> &vals,
			     const std::vector<double> &probs);
  
  static void make_cumulative(std::vector<double> &vals);
  
  static size_t find_bin(const std::vector<double> &bins,
			 const double val);

  static const size_t max_classes = 10000;
  static const double MIN_PROB;
  
  std::vector<double> breaks;
  std::vector<double> log_hist;
  std::vector<double> hist;
  std::vector<double> cumulative;
};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// DISCRETE EMPIRICAL

class DiscEmpDistro : public Distro_ {
public:
  
  DiscEmpDistro() : Distro_(std::vector<double>(2, 0)) {}
  DiscEmpDistro(std::vector<double> p) : Distro_(p) {}
  DiscEmpDistro(const DiscEmpDistro &rhs);
  DiscEmpDistro& operator=(const DiscEmpDistro &rhs);
  void set_params(const std::vector<double> &p) {params = p;}
  ~DiscEmpDistro() {}
  double sample() const;
  size_t required_params() const {return 2;}
  void estimate_params_ml(const std::vector<double> &vals);
  void estimate_params_ml(const std::vector<double> &vals,
			  const std::vector<double> &probs);
  void estimate_params_ml(const std::vector<double> &vals,
                          const std::vector<double> &scales,
                          const std::vector<double> &probs);
  double log_likelihood(double val) const;
    double log_likelihood(const double &val, const double &scale) const;

private:

  static size_t find_bin(const std::vector<double> &bins, const double val);
  static void make_hist(const std::vector<double> &data, size_t n_vals, 
			size_t n_classes, double max_val, std::vector<double> &hist);
  
  static void make_weighted_hist(const std::vector<double> &data,
				 const std::vector<double> &weights, size_t n_vals,
				 size_t n_classes, double max_val, 
				 std::vector<double> &hist);
  
  static void make_cumulative(std::vector<double> &vals);
  
  static const double MIN_PROB;
  
  size_t n_classes;
  double max_val;
  std::vector<double> log_hist;
  std::vector<double> hist;
  std::vector<double> cumulative;
};

#endif
