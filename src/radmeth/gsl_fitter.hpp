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

// Declares a function that fits regression models using GSL library.

#ifndef GSL_FITTER_HPP_
#define GSL_FITTER_HPP_

#include <vector>

class Regression;

bool gsl_fitter(Regression &r, 
              std::vector<double> initial_parameters = std::vector<double> ());

#endif // GSL_FITTER_HPP_
