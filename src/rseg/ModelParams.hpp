/*
 * Copyright (C) 2012 University of Southern California
 *                    Andrew D Smith and Qiang Song
 * Author: Qiang Song
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

#ifndef MODEL_PARAMS_HPP
#define MODEL_PARAMS_HPP

#include <string>
#include <vector>

template<class Distro_Type> void
read_param_file(const std::string &infile, size_t &n,
                std::vector<double> &start_trans,
                std::vector<std::vector<double> > &trans,
                std::vector<std::vector<double> > &end_trans,
                std::vector<Distro_Type> &distros);

template <class Distro_Type> void
write_param_file(const std::string &outfile,
                 const size_t &n,
                 const std::vector<double> &start_trans,
                 const std::vector<std::vector<double> > &trans,
                 const std::vector<std::vector<double> > &end_trans,
                 const std::vector<Distro_Type> &distros);

#include "ModelParams.cpp"

#endif

