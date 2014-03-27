
#ifndef PROXIMAL_LOCI_HPP_
#define PROXIMAL_LOCI_HPP_

#include "locus.hpp"

class ProximalLoci {
public:
  ProximalLoci(std::vector<LocusIterator> &loci, size_t max_distance)
    : loci_(loci), max_distance_(max_distance), next_pos_(loci.begin()) {};
  bool get(std::vector<LocusIterator> &neighbors);
  LocusIterator cur_region() {return *(next_pos_ - 1);}

private:
  const std::vector<LocusIterator> &loci_;
  size_t max_distance_;
  std::vector<LocusIterator>::const_iterator next_pos_;
};

#endif //PROXIMAL_LOCI_HPP_
