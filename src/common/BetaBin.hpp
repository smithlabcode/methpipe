/*
  Copyright (C) 2011 University of Southern California
  Authors: Andrew D. Smith, Song Qiang

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

#ifndef BETABIN_HPP
#define BETABIN_HPP

#include <utility>
#include <string>
#include <vector>

// struct betabin;
struct betabin 
{
    betabin();
    betabin(const double a, const double b);
    betabin(const std::string &str);
    double operator()(const std::pair<double, double> &val) const;
    double log_likelihood(const std::pair<double, double> &val) const;
    double sign(const double x);
    double invpsi(const double tolerance, const double x);
    double movement(const double curr, const double prev);
    void fit(const std::vector<double> &vals_a,
             const std::vector<double> &vals_b,
             const std::vector<double> &p);
    std::string tostring() const;
    double alpha;
    double beta;
    double lnbeta_helper;
    
    static const double tolerance;
};

#endif

