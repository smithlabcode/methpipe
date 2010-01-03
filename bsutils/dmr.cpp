/* dmrfinder: A program for identifying DMRs (differentially
 * methylated regions) based on a file showing probability of
 * differential methylation at each CpG or base.
 *
 * Copyright (C) 2009 University of Southern California
 *                    Andrew D Smith
 * Author: Andrew D. Smith
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

#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"

#include "FileIterator.hpp"

#include <gsl/gsl_sf_gamma.h>

#include <cmath>
#include <fstream>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;

using std::numeric_limits;
using std::ostream_iterator;
using std::ofstream;

class pred_significant_cpg
{
		double crit;
		bool LOWER_TAIL;
		bool UPPER_TAIL;
		bool get_significant;
public:
		pred_significant_cpg(double c, int TAIL_MARKER, bool get_sig_cpg = true) :
				crit(c),
				LOWER_TAIL(TAIL_MARKER & 1),
				UPPER_TAIL(TAIL_MARKER & 2),
				get_significant(get_sig_cpg) {}
		inline bool operator() (const GenomicRegion &cpg) const
		{
				bool is_significant = (LOWER_TAIL && cpg.get_score() < crit) ||
						(UPPER_TAIL && cpg.get_score() > 1 - crit);
				return get_significant ? is_significant : !is_significant;
		}
};

// seperate regions by chromosome and deserts
size_t get_reads_num(const string &name)
{
		return atoi(name.substr(name.find(':') + 1).c_str());
}

void
seperate_regions(vector<GenomicRegion> &cpgs,
				 vector<size_t> &reset_points,
				 const size_t desert_size)
{
// 		size_t j = 0;
// 		for (size_t i = 0; i < cpgs.size(); ++i)
// 				if (get_reads_num(cpgs[i].get_name()) > 0)
// 				{
// 						cpgs[j] = cpgs[i];
// 						++j;
// 				}
// 		cpgs.erase(cpgs.begin() + j, cpgs.end());
		
		reset_points.push_back(0);
		for (size_t i = 1; i < cpgs.size(); ++i)
		{
				const size_t dist = cpgs[i].distance(cpgs[i - 1]);
				if (dist > desert_size)
						reset_points.push_back(i);
		}
}


// funtions used to seek dmrs
// function decleration
void
seek_dmrs(const vector<GenomicRegion> &cpgs,
		  const vector<size_t> &reset_points,
		  vector<GenomicRegion> &dmrs,
		  const size_t width,
		  const double max_mismatch_ratio,
		  const double crit,
		  const int TAIL_MARKER);

void
seek_dmrs(const vector<GenomicRegion>::const_iterator &first,
		  const vector<GenomicRegion>::const_iterator &last,
		  vector<GenomicRegion> &dmrs,
		  const size_t min_width,
		  const double max_mismatch_ratio,
		  const double crit,
		  const int TAIL_MARKER);

// function definifition
// overal flow control
void
seek_dmrs(const vector<GenomicRegion> &cpgs,
		  const vector<size_t> &reset_points,
		  vector<GenomicRegion> &dmrs,
		  const size_t width,
		  const double max_mismatch_ratio,
		  const double crit,
		  const int TAIL_MARKER)
{
		for (size_t i = 0; i < reset_points.size() - 1; ++i)
				seek_dmrs(cpgs.begin() + reset_points[i],
						  cpgs.begin() + reset_points[i+1],
						  dmrs,
						  width,
						  max_mismatch_ratio,
						  crit,
						  TAIL_MARKER);
		return;
}

// greedy method to find dmrs in a consective series of CpGs 
void
seek_dmrs(const vector<GenomicRegion>::const_iterator &first,
		  const vector<GenomicRegion>::const_iterator &last,
		  vector<GenomicRegion> &dmrs,
		  const size_t min_width,
		  const double max_mismatch_ratio,
		  const double crit,
		  const int TAIL_MARKER)
{
		if (last - first < min_width) return; // range length is smaller than minimum width

		const size_t seed_size = 5;
		
		const pred_significant_cpg is_significant_cpg(crit, TAIL_MARKER);
		const pred_significant_cpg non_significant_cpg(crit, TAIL_MARKER, false);

		// find seed
		vector<GenomicRegion>::const_iterator seed_left
				= std::find_if(first, last, is_significant_cpg);
		vector<GenomicRegion>::const_iterator seed_right
				= std::find_if(seed_left, last, non_significant_cpg);
		while ((seed_left != last) && (seed_right - seed_left < seed_size))
		{
				seed_left = std::find_if(seed_right, last, is_significant_cpg);
				seed_right = std::find_if(seed_left, last, non_significant_cpg);
		}

		if (seed_left == last) return; // can not find a DMR seed in this range

		// extend seed to righthand side
		size_t mismatch = 0, count = seed_right - seed_left;
		while ((seed_right != last) && (static_cast<double>(mismatch) / count < max_mismatch_ratio))
		{

				++count;
				if (non_significant_cpg(*seed_right)) ++mismatch;
				seed_right++;
		}
		if (seed_right == last) seed_right--;
		while (non_significant_cpg(*seed_right))
		{
				--mismatch;
				--count;
				seed_right--; // ignore non sigsificant CpG's on right end
		}
		
		// extend seed to lefthand side
		while ((seed_left != first) && (static_cast<double>(mismatch) / count < max_mismatch_ratio))
		{

				++count;
				seed_left--;
				if (non_significant_cpg(*seed_left)) ++mismatch;
		}
		while (non_significant_cpg(*seed_left))
		{
				--mismatch;
				--count;
				++seed_left; // ignore non sigsificant CpG's on left end
		}
		
		// add new DMR
		if (count > min_width)
		{
				GenomicRegion region(*seed_left);
				region.set_end( seed_right->get_end() );
				double score = 0;
				for (vector<GenomicRegion>::const_iterator i = seed_left;
					 i < seed_right + 1;
					 ++i)
						score += i->get_score();
				region.set_score(score / ((seed_right - seed_left) + 1));
				dmrs.push_back(region);
		}
		
		seek_dmrs(seed_right + 1, last,
				  dmrs,
				  min_width, max_mismatch_ratio,
				  crit, TAIL_MARKER);
		return;
}

// // brutal force N^2
// void
// seek_dmrs(const vector<GenomicRegion>::const_iterator &first,
// 		  const vector<GenomicRegion>::const_iterator &last,
// 		  vector<GenomicRegion> &dmrs,
// 		  const size_t min_width,
// 		  const double max_mismatch_ratio,
// 		  const double crit,
// 		  const int TAIL_MARKER)
// {
// 		if (last - first < min_width) return; // range length is smaller than minimum width

// 		const size_t seed_size = 6;
		
// 		const pred_significant_cpg is_significant_cpg(crit, TAIL_MARKER);
// 		const pred_significant_cpg non_significant_cpg(crit, TAIL_MARKER, false);

// 		// find seed
// 		vector<GenomicRegion>::const_iterator seed_left
// 				= std::find_if(first, last, is_significant_cpg);
// 		vector<GenomicRegion>::const_iterator seed_right
// 				= std::find_if(seed_left, last, non_significant_cpg);
// 		while ((seed_left != last) && (seed_right - seed_left < seed_size))
// 		{
// 				seed_left = std::find_if(seed_right, last, is_significant_cpg);
// 				seed_right = std::find_if(seed_left, last, non_significant_cpg);
// 		}

// 		if (seed_left == last) return; // can not find a DMR seed in this range

// 		// exhausitive search
		
// 		for ()
		
// 		// add new DMR
// 		if (count > min_width)
// 		{
// 				GenomicRegion region(*seed_left);
// 				region.set_end( seed_right->get_end() );
// 				score = 0;
// 				for (vector<GenomicRegion>::const_iterator i = seed_left;
// 					 i < seed_right + 1;
// 					 ++i)
// 						score += i->get_score();
// 				region.set_score(score / ((seed_right - seed_left) + 1));
// 				dmrs.push_back(region);
// 		}
		
// 		seek_dmrs(seed_right + 1, last,
// 				  dmrs,
// 				  min_width, max_mismatch_ratio,
// 				  crit, TAIL_MARKER);
// 		return;
// }



int
main(int argc, const char **argv) {
		try {

				string outfile;
				double crit = 0.05;
				size_t min_cpgs = 20;
				double max_mismatch_ratio = 0.1;
				int TAIL_MARKER = 1;

				size_t desert_size = 2000;
				bool VERBOSE = false;

    
				/****************** COMMAND LINE OPTIONS ********************/
				OptionParser opt_parse("dmr", 
									   "A program for identifying DMRs (differentially "
									   "methylated regions) based on a file showing "
									   "probability of differential methylation at "
									   "each CpG or base. ",
									   "<cpg_meth_diffs_file>");
				opt_parse.add_opt("crit", 'c', "critical value (default: 0.05)", 
								  false, crit);
				opt_parse.add_opt("width", 'w', "width in terms of CpGs of min DMR size (default: 20)", 
								  false, min_cpgs);
				opt_parse.add_opt("mismatch ratio", 'r', "The maximum proportion of non-significant CpGs allowed in DMR (default 0.1)",
								  false, max_mismatch_ratio);
				opt_parse.add_opt("tail", 't',
								  "use which tail to determine signif.\n 1: LOWER TAIL; 2: UPPER_TAIL; 3: BOTH",
								  false, TAIL_MARKER);
				opt_parse.add_opt("desert", 'd', "size of desert",
								  false, desert_size);
				opt_parse.add_opt("out", 'o', "output file (BED format)", 
								  false, outfile);
				opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
				vector<string> leftover_args;
				opt_parse.parse(argc, argv, leftover_args);
				if (argc == 1 || opt_parse.help_requested()) {
						cerr << opt_parse.help_message() << endl
							 << opt_parse.about_message() << endl;
						return EXIT_SUCCESS;
				}
				if (opt_parse.about_requested()) {
						cerr << opt_parse.about_message() << endl;
						return EXIT_SUCCESS;
				}
				if (opt_parse.option_missing()) {
						cerr << opt_parse.option_missing_message() << endl;
						return EXIT_SUCCESS;
				}
				if (leftover_args.size() != 1) {
						cerr << opt_parse.help_message() << endl;
						return EXIT_SUCCESS;
				}
				const string cpgs_file = leftover_args.front();
				/****************** END COMMAND LINE OPTIONS *****************/

				// read in  methdiff output
				if (VERBOSE)
						cerr << "Reading CpGs ... ";
				vector<GenomicRegion> cpgs;
				ReadBEDFile(cpgs_file, cpgs);
				if (!check_sorted(cpgs))
						throw RMAPException("file not sorted: \"" + cpgs_file + "\"");
				if (VERBOSE)
						cerr << "Done" << endl;
    
				// seperate regions by chromosome and remove deserts
				if (VERBOSE)
						cerr << "Seperating regions ... ";
				vector<size_t> reset_points;
				seperate_regions(cpgs, reset_points, desert_size);
				if (VERBOSE)
						cerr << "done" << endl;

				// find DMR's
				if (VERBOSE)
						cerr << "Seeking differentially methylated regions ... ";
				vector<GenomicRegion> regions;
				seek_dmrs(cpgs, reset_points, regions,
						  min_cpgs, max_mismatch_ratio,
						  crit, TAIL_MARKER);
				if (VERBOSE)
						cerr << "done" << endl;

				// output DMR
				if (VERBOSE)
						cerr << "Writing result ... ";
				std::ostream *out = (outfile.empty()) ? &cout : 
						new std::ofstream(outfile.c_str());
				copy(regions.begin(), regions.end(), 
					 ostream_iterator<GenomicRegion>(*out, "\n"));
				if (out != &cout) delete out;
				if (VERBOSE)
						cerr << "done" << endl;
		}
		catch (RMAPException &e) {
				cerr << "ERROR:\t" << e.what() << endl;
				return EXIT_FAILURE;
		}
		catch (std::bad_alloc &ba) {
				cerr << "ERROR: could not allocate memory" << endl;
				return EXIT_FAILURE;
		}
		return EXIT_SUCCESS;
} 
