 
#include <sstream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <iostream>
#include <exception>

#include "locus.hpp"

using std::string; using std::vector;

const double max_abs_error = 0.001;

Locus::Locus(std::string encoding) {
  try {
  std::istringstream iss(encoding);
  iss.exceptions(std::ios::failbit);
  iss >> chrom_ >> begin_ >> end_ >> name_ >> score_;
  } catch (std::exception const & err) {
    std::cerr << err.what() << std::endl
              << "Couldn't parse the line \"" << encoding << "\"." << std::endl;
    std::terminate();
  }
  
  //if (!iss) {
  //  std::string message("Cannot parse \"" + encoding + "\";" 
  //                        + " at least five columns are required.");
  //  throw std::logic_error(message);
  //}
  
  if (end_ <= begin_)
    throw std::logic_error("in \"" + encoding + "\" first coordinate should " 
                            + "be less than second");
}

bool 
Locus::operator== (const Locus &other_locus) const {
  return (chrom_ == other_locus.chrom_) &&
         (begin_ == other_locus.begin_) &&
         (end_ == other_locus.end_) && 
         (name_ == other_locus.name_) &&
         (fabs(score_ - other_locus.score_) <= max_abs_error);
}

void
read_loci(std::istream &loci_encoding, std::vector<Locus> &loci) {
  loci.clear();
  
  string locus_encoding;
  
  while( getline(loci_encoding, locus_encoding) ) {
    loci.push_back(Locus(locus_encoding));
  }
}

std::ostream&
operator<<(std::ostream& os, const Locus &locus) {
  os << locus.chrom_ << "\t"
     << locus.begin_ << "\t"
     << locus.end_ << "\t"
     << locus.name_ << "\t"
     << locus.score_;
  return os;
}

void 
get_iterators_to_good_loci(vector<Locus> &input_loci, 
                           vector<LocusIterator> &good_loci_iterators) {
  for (LocusIterator it = input_loci.begin(); it != input_loci.end(); ++it) {
    const double score = it->score();
    if (0.0 <= score && score <= 1.0) {
      good_loci_iterators.push_back(it);
    } 
  }
}


