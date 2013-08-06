/*
 *    Part of RMAP software
 *
 *    Copyright (C) 2009 University of Southern California and
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

#ifndef LOAD_PAIRED_END_READS_HPP
#define LOAD_PAIRED_END_READS_HPP

#include "smithlab_utils.hpp"
#include "FastRead.hpp"
#include "FastReadQuality.hpp"
#include "FastReadWC.hpp"

#include <vector>
#include <string>

void
load_reads_from_fasta_file(const std::string &filename, const size_t max_diffs,
			   size_t &read_width, 
			   std::vector<FastRead> &fast_reads_left,
			   std::vector<FastRead> &fast_reads_right,
			   std::vector<size_t> &read_words, std::vector<size_t> &read_index);

void
load_reads_from_fastq_file(const std::string &filename, const size_t max_diffs,
			   size_t &read_width, 
			   std::vector<FastRead> &fast_reads_left,
			   std::vector<FastRead> &fast_reads_right,
			   std::vector<size_t> &read_words, std::vector<size_t> &read_index);

void
load_reads_from_fastq_file(const std::string &filename, const size_t max_diffs,
			   size_t &read_width, 
			   std::vector<FastReadWC> &fast_reads_left,
			   std::vector<FastReadWC> &fast_reads_right,
			   std::vector<size_t> &read_words, std::vector<size_t> &read_index);

void
load_reads_from_fastq_file(const std::string &filename, const size_t max_diffs,
			   size_t &read_width, 
			   std::vector<FastReadQuality> &fast_reads_left,
			   std::vector<FastReadQuality> &fast_reads_right,
			   std::vector<size_t> &read_words, std::vector<size_t> &read_index);

void
load_reads_from_prb_file(const std::string &filename, const size_t max_diffs,
			 size_t &read_width, 
			 std::vector<FastReadWC> &fast_reads_left,
			 std::vector<FastReadWC> &fast_reads_right,
			 std::vector<size_t> &read_words, std::vector<size_t> &read_index);

void
load_reads_from_prb_file(const std::string &filename, const size_t max_diffs,
			 size_t &read_width, 
			 std::vector<FastReadQuality> &fast_reads_left,
			 std::vector<FastReadQuality> &fast_reads_right,
			 std::vector<size_t> &read_words, std::vector<size_t> &read_index);


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void
load_reads_from_fasta_file_cheat(const std::string &filename, const size_t max_diffs,
			   size_t &read_width, 
			   std::vector<FastRead> &fast_reads_left,
			   std::vector<FastRead> &fast_reads_right,
			   std::vector<size_t> &read_words_left, std::vector<size_t> &read_words_right, 
			   std::vector<size_t> &read_index_left, std::vector<size_t> &read_index_right);

void
load_reads_from_fastq_file_cheat(const std::string &filename, const size_t max_diffs,
			   size_t &read_width, 
			   std::vector<FastRead> &fast_reads_left,
			   std::vector<FastRead> &fast_reads_right,
			   std::vector<size_t> &read_words_left, std::vector<size_t> &read_words_right, 
			   std::vector<size_t> &read_index_left, std::vector<size_t> &read_index_right);

void
load_reads_from_fastq_file_cheat(const std::string &filename, const size_t max_diffs,
			   size_t &read_width, 
			   std::vector<FastReadWC> &fast_reads_left,
			   std::vector<FastReadWC> &fast_reads_right,
			   std::vector<size_t> &read_words_left, std::vector<size_t> &read_words_right, 
			   std::vector<size_t> &read_index_left, std::vector<size_t> &read_index_right);

void
load_reads_from_fastq_file_cheat(const std::string &filename, const size_t max_diffs,
			   size_t &read_width, 
			   std::vector<FastReadQuality> &fast_reads_left,
			   std::vector<FastReadQuality> &fast_reads_right,
			   std::vector<size_t> &read_words_left, std::vector<size_t> &read_words_right, 
			   std::vector<size_t> &read_index_left, std::vector<size_t> &read_index_right);

void
load_reads_from_prb_file_cheat(const std::string &filename, const size_t max_diffs,
			 size_t &read_width, 
			 std::vector<FastReadWC> &fast_reads_left,
			 std::vector<FastReadWC> &fast_reads_right,
			 std::vector<size_t> &read_words_left, std::vector<size_t> &read_words_right, 
			 std::vector<size_t> &read_index_left, std::vector<size_t> &read_index_right);

void
load_reads_from_prb_file_cheat(const std::string &filename, const size_t max_diffs,
			 size_t &read_width, 
			 std::vector<FastReadQuality> &fast_reads_left,
			 std::vector<FastReadQuality> &fast_reads_right,
			 std::vector<size_t> &read_words_left, std::vector<size_t> &read_words_right, 
			 std::vector<size_t> &read_index_left, std::vector<size_t> &read_index_right);

#endif
