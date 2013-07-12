/*
 *    Copyright (C) 2012 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
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
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BAMTOOLS_INTERFACE_HPP
#define BAMTOOLS_INTERFACE_HPP

#include <string>
#include <vector>

#include <MappedRead.hpp>
#include <tr1/unordered_map>

#include <api/BamReader.h>
#include <api/BamAlignment.h>

void
ReadBAMFormatInput(const std::string &infile, std::vector<MappedRead> &reads);

void
BamAlignmentToMappedReadWithMapper(
    const std::tr1::unordered_map<size_t, std::string> &chrom_lookup,
    const BamTools::BamAlignment &ba, MappedRead &mr, std::string mapper);

#endif
