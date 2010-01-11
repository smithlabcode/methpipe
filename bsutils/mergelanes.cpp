/*    mergelanes: a program for sorting read sequences relative to an order
 *    given in a BED file (i.e. mapping locations)
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

#include <string>
#include <vector>
#include <unistd.h>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <queue>

#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"
#include "RNG.hpp"

#include <sys/time.h>
#include <sys/resource.h>

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::numeric_limits;
using std::ofstream;


template <class T>
class FileIterator
{
public:
		FileIterator(const std::string &f, const size_t bs);
		void increment_first()
		{
				if (++first == buffer.end()) {
						refill_buffer();
				}
		}
		typename std::vector<T>::const_iterator get_first() const {return first;}
		bool first_is_good() const {return (!in.eof() || first < buffer.end());}
  
private:
		std::ifstream in;
		std::vector<T> buffer;
		typename std::vector<T>::iterator first;
		void refill_buffer();
};

std::istream& 
operator>>(std::istream &in, GenomicRegion &region)
{
		static const size_t buffer_size = 10000; // Magic
		char buffer[buffer_size];
		in.getline(buffer, buffer_size);
		region = GenomicRegion(buffer);
		return in;
}

/*  THIS FUNCTION FILLS A BUFFER FOR GenomicRegion OBJECTS
 */
void 
fill_buffer(std::ifstream &in, const size_t buffer_start, 
			vector<GenomicRegion> &buffer)
{
		GenomicRegion tmp;
		size_t i = buffer_start;
		assert(buffer_start <= buffer.size());
		for (; i != buffer.size() && !in.eof(); ++i)
		{
				in >> tmp;
				buffer[i].swap(tmp);
				in.peek();
		}
		if (i < buffer.size())
				buffer.erase(buffer.begin() + i, buffer.end());
}

/*  THIS FUNCTION FILLS A BUFFER FOR THE ACTUAL READS, REPRESENTED AS
 *  STRINGS, AND MUST BE IN A FASTA FORMAT FILE
 */
void 
fill_buffer(std::ifstream &in, const size_t buffer_start, 
			vector<string> &buffer) {
		string tmp;
		size_t i = buffer_start;
		for (; i != buffer.size() && !in.eof(); ++i) {
				// the read name...
				in >> tmp; // DANGER: assumes that the name of the read has no
				// spaces in it!!
				// the read itself:
				in >> buffer[i];
				in.peek();
		}
		if (i < buffer.size())
				buffer.erase(buffer.begin() + i, buffer.end());
}


/*  THIS FUNCTION FILLS A BUFFER FOR THE ACTUAL READS, REPRESENTED AS
 *  STRINGS, AND MUST BE IN A FASTA FORMAT FILE
 */
struct FASTQRecord {
		FASTQRecord() {}
		FASTQRecord(const FASTQRecord &rhs) : 
				name(rhs.name), seq(rhs.seq), qual(rhs.qual) {}
		FASTQRecord(const string &n_,
					const string &s_,
					const string &q_) : name(n_), seq(s_), qual(q_) {}
		string name;
		string seq;
		string qual;
};

std::ostream&
operator<<(std::ostream &s, const FASTQRecord &rhs) {
		return s << rhs.name << '\n'
				 << rhs.seq << "\n+\n"
				 << rhs.qual;
}

void 
fill_buffer(std::ifstream &in, const size_t buffer_start, 
			vector<FASTQRecord> &buffer) 
{
		string tmp, read_name, read_seq, scores_seq;
		size_t i = buffer_start;
		for (; i != buffer.size() && !in.eof(); ++i)
		{
				// the read name...
				in >> read_name; // DANGER: assumes that the name of the read has no
				// spaces in it!!
				// the read itself:
				in >> read_seq;
				in >> tmp;
				in >> scores_seq;
				buffer[i] = FASTQRecord(read_name, read_seq, scores_seq);
				in.peek();
		}
		if (i < buffer.size())
				buffer.erase(buffer.begin() + i, buffer.end());
}


/* THIS REFILL BUFFER IS USED WHEN INCREMENTS TO EITHER THE FIRST OR
   THE LAST CURRENTLY USED ELEMENTS IN THE BUFFER HIT THE END OF THE
   BUFFER. HOPEFULLY THE FIRST ONE WILL NOT HIT THE END BEFORE THE
   LAST: THE FIRST IS ALWAYS SUPPOSED TO BE LESS THAN OR EQUAL TO THE
   LAST.
*/
template <class T> void
FileIterator<T>::refill_buffer() 
{
		first = buffer.begin();
		fill_buffer(in, 0, buffer);
}

template <class T>
FileIterator<T>::FileIterator(const std::string &f, const size_t bs) :
		buffer(std::vector<T>(bs)) 
{
		in.open(f.c_str());
		if (!in)
				throw RMAPException("cannot open input file " + f);
		fill_buffer(in, 0, buffer);
		first = buffer.begin();
}

struct ComparePairs 
{
		bool operator()(const pair<GenomicRegion, size_t> &a,
						const pair<GenomicRegion, size_t> &b) const {
				return !(a.first < b.first);
		}
};


int 
main(int argc, const char **argv) {

		try {

				bool VERBOSE = false;
				string map_outfile;
				string read_outfile;
				size_t random_number_seed = numeric_limits<size_t>::max();
    
				size_t BUFFER_SIZE = 10000ul;
    
				/****************** COMMAND LINE OPTIONS ********************/
				OptionParser opt_parse("mergelanes", 
									   "A program to merge different lanes in a shotgun bisulfite "
									   "sequencing experiment producing one large sorted file of all read "
									   "mapping locations (without duplicate 5' ends) and one large file "
									   "sorted similarly containing the corresponding sequences.",
									   "<fasta-reads-file-1> [<fasta-reads-file-2> ...]");
				opt_parse.add_opt("output", 'o', "Name of maps output file", 
								  true, map_outfile);
				opt_parse.add_opt("readout", 'r', "Name of reads output file", 
								  true, read_outfile);
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
				if (leftover_args.empty()) {
						cerr << opt_parse.help_message() << endl;
						return EXIT_SUCCESS;
				}
				const string mapped_file_names_file(leftover_args.front());
				vector<string> mapped_files;
				if (!mapped_file_names_file.empty())
						read_filename_file(mapped_file_names_file.c_str(), mapped_files);
				else {
						cerr << opt_parse.help_message() << endl;
						return EXIT_SUCCESS;
				}
				const string read_file_names_file(leftover_args.back());
				vector<string> reads_files;
				if (!read_file_names_file.empty())
						read_filename_file(read_file_names_file.c_str(), reads_files);
				else {
						cerr << opt_parse.help_message() << endl;
						return EXIT_SUCCESS;
				}
				/****************** END COMMAND LINE OPTIONS *****************/
    
				vector<FileIterator<GenomicRegion> *> itrs;
				if (VERBOSE)
						cerr << "[MAPPED READS FILES:]" << endl;
				for (size_t i = 0; i < mapped_files.size(); ++i)
				{
						if (VERBOSE)
								cerr << mapped_files[i] << endl;
						itrs.push_back(new FileIterator<GenomicRegion>(mapped_files[i], BUFFER_SIZE));
				}
    
				vector<FileIterator<FASTQRecord> *> read_itrs;
				if (VERBOSE)
						cerr << "[FASTQ READS FILES:]" << endl;
				for (size_t i = 0; i < reads_files.size(); ++i)
				{
						if (VERBOSE)
								cerr << reads_files[i] << endl;
						read_itrs.push_back(new FileIterator<FASTQRecord>(reads_files[i], BUFFER_SIZE));
				}
    
				std::priority_queue<pair<GenomicRegion, size_t>, 	
					vector<pair<GenomicRegion, size_t> >, ComparePairs> a;
				for (size_t i = 0; i < itrs.size(); ++i)
						a.push(make_pair(*itrs[i]->get_first(), i));
    
				ofstream out(map_outfile.c_str());
				ofstream read_out(read_outfile.c_str());

				vector<GenomicRegion> mapped_ties;
				vector<FASTQRecord> reads_ties;
				double score = std::numeric_limits<double>::max();
				const Runif rng(random_number_seed);
    
				while (!a.empty())
				{
						const size_t file_id = a.top().second;
						if (mapped_ties.empty() || 
							!(mapped_ties.front() < a.top().first)) // the same mapped location
						{
								double new_score = a.top().first.get_score(); // reads of better quality
								if (new_score < score)
								{
										new_score = a.top().first.get_score();
										mapped_ties.clear();
										reads_ties.clear();
								}
						} else {
								const size_t rand_idx = rng.runif(0ul, mapped_ties.size());
								out << mapped_ties[rand_idx] << '\n';
								read_out << reads_ties[rand_idx] << '\n';
								mapped_ties.clear();
								reads_ties.clear();
								score = std::numeric_limits<double>::max();
						}
						mapped_ties.push_back(a.top().first);
						reads_ties.push_back(*read_itrs[file_id]->get_first());
						a.pop();
						itrs[file_id]->increment_first();
						read_itrs[file_id]->increment_first();
						if (itrs[file_id]->first_is_good())
								a.push(make_pair(*itrs[file_id]->get_first(), file_id));
				}
				const size_t rand_idx = rng.runif(0ul, mapped_ties.size());
				out << mapped_ties[rand_idx] << '\n';
				read_out << reads_ties[rand_idx] << '\n';
		}
		catch (const RMAPException &e) {
				cerr << e.what() << endl;
				return EXIT_FAILURE;
		}
		catch (std::bad_alloc &ba) {
				cerr << "ERROR: could not allocate memory" << endl;
				return EXIT_FAILURE;
		}
		return EXIT_SUCCESS;
}
