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

#ifndef TWO_STATE_RESOLVE_MIXTURE_HPP
#define TWO_STATE_RESOLVE_MIXTURE_HPP

#include "Distro.hpp"
void
TwoStateResolveMixture(
    const std::vector<double> &values,
    const std::vector<double> &scales,
    const size_t max_iterations,
    const double tolerance, int VERBOSE,
    Distro &fg_distro, Distro &bg_distro, double &mixing);

#endif
