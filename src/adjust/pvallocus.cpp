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

#include <iostream>
#include <sstream>

#include "pvallocus.hpp"

using std::vector;
using std::string;
using std::istream;
using std::ostream;

void 
initialize_pval_loci(istream &encoding, vector<PvalLocus> &pval_loci) {
	string record;
	string chrom;
	size_t begin;
	size_t end;
	string name;
	double pval;

	string prev_chrom;
	size_t chrom_index = 0;

	while(getline(encoding, record)) {
		//std::cout << record << std::endl;
		
		try {
  			std::istringstream iss(record);
  			iss.exceptions(std::ios::failbit);
  			iss >> chrom >> begin >> end >> name >> pval;
  		} catch (std::exception const & err) {
    		std::cerr << err.what() << std::endl << "Couldn't parse the line \""
    				  << record << "\"." << std::endl;
    		std::terminate();
  		}

  		// All loci not associated with valide p-values are skipped;
  		if (0 <= pval && pval <= 1) {
  			
  			if (prev_chrom.empty())
  				prev_chrom = chrom;

  			if (chrom != prev_chrom) {
  				prev_chrom = chrom;
  				++chrom_index;
  			}

			PvalLocus plocus(chrom_index, begin, pval);
  			pval_loci.push_back(plocus);
  		}
	}
}

void
update_pval_loci(istream &input_encoding, 
								 const vector<PvalLocus> &pval_loci, 
								 ostream &output_encoding) {
	string record;
	string chrom;
	size_t begin;
	size_t end;
	string name;
	double pval;

	string prev_chrom;
	size_t chrom_index = 0;

	vector<PvalLocus>::const_iterator cur_locus_iter = pval_loci.begin();

	while(getline(input_encoding, record)) {
		
		try {
  		std::istringstream iss(record);
  		iss.exceptions(std::ios::failbit);
  		iss >> chrom >> begin >> end >> name >> pval;
  	} catch (std::exception const & err) {
    	std::cerr << err.what() << std::endl << "Couldn't parse the line \"" 
    				  	<< record << "\"." << std::endl;
    	std::terminate();
  	}

    output_encoding << chrom << "\t" << begin << "\t" << end << "\t" << name;

    if (0 <= pval && pval <= 1) {
      output_encoding << ":" << cur_locus_iter->raw_pval << ":" 
                      << cur_locus_iter->combined_pval << "\t"
                      << cur_locus_iter->corrected_pval;
      cur_locus_iter++;
    } else {
      output_encoding << "\t" << pval;
    }

    if (cur_locus_iter != pval_loci.end())
      output_encoding << "\n";

  	}
}
