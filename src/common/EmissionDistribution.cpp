/*
  Copyright (C) 2017 University of Southern California
  Authors: Andrew D. Smith and Benjamin E. Decato

  This file is part of methpipe.

  methpipe is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  methpipe is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with rmap; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include "EmissionDistribution.hpp"

using std::vector;
using std::pair;
using std::setw;
using std::max;
using std::min;
using std::cerr;
using std::endl;
using std::string;
using std::setprecision;

EmissionDistribution::EmissionDistribution() :
    alpha(1), beta(1), lnbeta_helper(gsl_sf_lnbeta(1, 1)) {}

EmissionDistribution::EmissionDistribution(const double a, const double b) :
    alpha(a), beta(b), lnbeta_helper(gsl_sf_lnbeta(a, b)) {}

EmissionDistribution::EmissionDistribution(const string &str) {
    std::istringstream iss(str, std::istringstream::in);
    string name;
    iss >> name >> alpha >> beta;
    if (name != "edtn" || alpha < 0 || beta < 0)
    {
        cerr << "EmissionDistribution::EmissionDistribution: "
             << "bad string representation of emission distribution: "
             << str << endl;
        throw "bad string representation of emission distribution";
    }
    lnbeta_helper = gsl_sf_lnbeta(alpha, beta);
}

EmissionDistribution::~EmissionDistribution() {}

string
EmissionDistribution::tostring() const {
    std::ostringstream os;
    os << "Emission dtn params: " << setprecision(4) << alpha << " "
       << setprecision(4) << beta;
    return os.str();
}


double
EmissionDistribution::sign(const double x) {
    return (x >= 0) ? 1.0 : -1.0;
}


double
EmissionDistribution::invpsi(const double tolerance, const double x) {
    double L = 1.0, Y = std::exp(x);
    while (L > tolerance)
    {
        Y += L*sign(x - gsl_sf_psi(Y));
        L /= 2.0;
    }
    return Y;
}


double
EmissionDistribution::movement(const double curr, const double prev) {
    return std::abs(curr - prev)/std::max(std::fabs(curr), std::fabs(prev));
}


void
EmissionDistribution::fit(const vector<double> &vals_a,
                      const vector<double> &vals_b, const vector<double> &p) {
    const double p_total = std::accumulate(p.begin(), p.end(), 0.0);
    const double alpha_rhs = inner_product(vals_a.begin(), vals_a.end(),
                                           p.begin(), 0.0)/p_total;
    const double beta_rhs = inner_product(vals_b.begin(), vals_b.end(),
                                          p.begin(), 0.0)/p_total;

    double prev_alpha = 0.0, prev_beta = 0.0;
    alpha = beta = 0.01;
    while (movement(alpha, prev_alpha) > tolerance &&
           movement(beta, prev_beta) > tolerance)
    {
        prev_alpha = alpha;
        prev_beta = beta;
        alpha = invpsi(tolerance, gsl_sf_psi(prev_alpha + prev_beta) + alpha_rhs);
        beta = invpsi(tolerance, gsl_sf_psi(prev_alpha + prev_beta) + beta_rhs);
    }
    lnbeta_helper = gsl_sf_lnbeta(alpha, beta);
}

Beta::Beta() : EmissionDistribution() {}
Beta::Beta(const double a, const double b) : EmissionDistribution(a,b) {}
Beta::Beta(const std::string &str) : EmissionDistribution(str) {}

double
Beta::operator()(const pair<double, double> &val) const {
    const double p = val.first/val.second;
    return (alpha-1.0)*log(p) + (beta-1.0)*log(1.0-p) - gsl_sf_lnbeta(alpha, beta);
}

double
Beta::log_likelihood(const pair<double, double> &val) const {
    const double p = val.first/val.second;
    return (alpha-1.0)*log(p) + (beta-1.0)*log(1.0-p) - gsl_sf_lnbeta(alpha, beta);
}

BetaBinomial::BetaBinomial() : EmissionDistribution() {}
BetaBinomial::BetaBinomial(const double a, const double b)
                             : EmissionDistribution(a,b) {}
BetaBinomial::BetaBinomial(const std::string &str)
                             : EmissionDistribution(str) {}

double
BetaBinomial::operator()(const pair<double, double> &val) const {
    const size_t x = static_cast<size_t>(val.first);
    const size_t n = static_cast<size_t>(x + val.second);
    return gsl_sf_lnchoose(n, x) +
        gsl_sf_lnbeta(alpha + x, beta + val.second) - lnbeta_helper;
}

double
BetaBinomial::log_likelihood(const pair<double, double> &val) const {
    const size_t x = static_cast<size_t>(val.first);
    const size_t n = static_cast<size_t>(x + val.second);
    return gsl_sf_lnchoose(n, x) +
        gsl_sf_lnbeta(alpha + x, beta + val.second) - lnbeta_helper;
}
