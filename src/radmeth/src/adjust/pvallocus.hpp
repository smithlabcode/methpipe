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
 
#ifndef PVALLOCUS_HPP_
#define PVALLOCUS_HPP_

#include <vector>
#include <iostream>

struct PvalLocus {
	PvalLocus(std::size_t new_chrom_ind, std::size_t new_pos, double new_raw_pval)
		: chrom_ind(new_chrom_ind), pos(new_pos), raw_pval(new_raw_pval),
		  combined_pval(0), corrected_pval(0) {}
	std::size_t chrom_ind;
	std::size_t pos;
	double raw_pval;
	double combined_pval;
	double corrected_pval;
};

void initialize_pval_loci(std::istream &encoding, 
						 							std::vector<PvalLocus> &pval_loci);

void update_pval_loci(std::istream &input_encoding, 
								 			const std::vector<PvalLocus> &pval_loci, 
								 			std::ostream &output_loci_encoding);

#endif //PVALLOCUS_HPP_