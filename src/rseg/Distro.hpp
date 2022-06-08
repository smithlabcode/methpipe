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

#ifndef DISTRO_HPP
#define DISTRO_HPP

#include <algorithm>
#include <vector>
#include <iostream>
#include <string>

class Distro_ {
public:

  Distro_();
  Distro_(const std::vector<double> p);
  Distro_(const Distro_ &);
  Distro_& operator=(const Distro_ &);
  virtual ~Distro_();

  virtual size_t required_params() const = 0;
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

  void set_params(const std::vector<double> &p) {params = p;}


  static double
  log_sum_log_vec(const std::vector<double> &vals, size_t limit);

protected:
  std::vector<double> params;
  std::vector<double> workspace_vals;
  std::vector<double> workspace_probs;

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
// POISSON


class PoisDistro : public Distro_ {
public:
  PoisDistro() : Distro_(std::vector<double>(1, 0)) {}
  PoisDistro(std::vector<double> p) : Distro_(p) {}
  PoisDistro(const PoisDistro &rhs);
  PoisDistro& operator=(const PoisDistro &rhs);
  ~PoisDistro() {}
  size_t required_params() const {return 1;}
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

#endif
