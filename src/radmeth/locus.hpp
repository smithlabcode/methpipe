
#ifndef LOCUS_HPP_
#define LOCUS_HPP_

#include <string>
#include <sstream>
#include <vector>

class Locus {
public:
  Locus() : chrom_(""), begin_(0), end_(0), name_(""), score_(0) {}
  Locus(std::string encoding);
  
  std::string chrom() const {return chrom_;}
  size_t begin() const {return begin_;}
  size_t end() const {return end_;}
  std::string name() const {return name_;}
  double score() const {return score_;}
  void set_chrom(std::string chrom) { chrom_ = chrom; }
  void set_begin(size_t begin) { begin_ = begin; }
  void set_end(size_t end) { end_ = end; }
  void set_name(std::string name) { name_ = name; }  
  void set_score(double score) { score_ = score; }
  bool operator== (const Locus &other_locus) const;
  friend std::ostream& operator<<(std::ostream& os, const Locus &locus);
private:
  std::string chrom_;
  size_t begin_;
  size_t end_;
  std::string name_;
  double score_;
};

void read_loci(std::istream &loci_encoding, std::vector<Locus> &loci);

typedef std::vector<Locus>::iterator LocusIterator;

void get_iterators_to_good_loci(std::vector<Locus> &input_loci, 
                              std::vector<LocusIterator> &good_loci_iterators);

#endif //LOCUS_HPP_
