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

// Defines the wand function which is described below.

#ifndef PIPELINE_HPP_
#define PIPELINE_HPP_

#include <sstream>
#include <string>

// The wand function executes the beta-binomial regresson pipeline. Parameters 
// design_encoding and count_table_encoding are streams with the design matrix 
// and the count table respectively. The test_factor_name parameter is the name 
// of the factor with respect to which the differential methylation should be 
// tested. All of the output is sent to the out output stream.

void wand(std::istream &design_encoding, std::istream &count_table_encoding, 
          std::string test_factor_name, std::ostream &out);

#endif //PIPELINE_HPP_
