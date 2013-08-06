/*
 *    Part of RMAP software
 *
 *    Copyright (C) 2008 Cold Spring Harbor Laboratory, 
 *                       University of Southern California and
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

#ifndef RMAP_UTILS_HPP
#define RMAP_UTILS_HPP

typedef size_t MASK_t;

namespace rmap_bits {
  static const MASK_t low_bit = 1ul;
  static const MASK_t high_bit = static_cast<size_t>(0x8000000000000000ul);
  static const MASK_t all_ones = static_cast<size_t>(-1);
  static const MASK_t all_zeros = 0ul;
  static const size_t word_size = 64ul;
}

#endif
