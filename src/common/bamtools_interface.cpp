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

#include "bamtools_interface.hpp"

#include <tr1/unordered_map>
#include <string>
#include <vector>

#include "MappedRead.hpp"
#include "GenomicRegion.hpp"
#include "SAM.hpp"

#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <api/BamAux.h>

using std::tr1::unordered_map;
using std::string;
using std::vector;

using std::cerr;
using std::endl;

using BamTools::BamAlignment;
using BamTools::SamHeader;
using BamTools::RefVector;
using BamTools::BamReader;
using BamTools::RefData;
using BamTools::CigarOp;

void
apply_CIGAR_for_BAM(const string &seq, const string &qual,
    vector<CigarOp> CigarData, string &new_seq, string &new_qual) {
  assert(seq.size() == qual.size());
  assert(new_seq.size() == 0 && new_qual.size() == 0);
  size_t n;
  char op;
  size_t i = 0;

  for(size_t j = 0; j < CigarData.size(); j++) {
    op = CigarData[j].Type;
    n = CigarData[j].Length;
    switch (op)
    {
    case 'M':
      new_seq += seq.substr(i, n);
      new_qual += qual.substr(i, n);
      i += n;
      break;
    case 'I':
      i += n;
      break;
    case 'D':
      new_seq += string(n, 'N');
      new_qual += string(n, 'B');
      break;
    case 'S':
      i += n;
      break;
    case 'H':
      ;
      break;
    case 'P':
      ;
      break;
    case '=':
      ;
      break;
    case 'X':
      ;
      break;
    }
  }
  
  assert(i == seq.length());
  assert(new_seq.size() == new_qual.size());
}

void
BamAlignmentToMappedReadWithMapper(
    const unordered_map<size_t, string> &chrom_lookup, const BamAlignment &ba,
    MappedRead &mr, string mapper) {

  const unordered_map<size_t, string>::const_iterator 
    the_chrom(chrom_lookup.find(ba.RefID));
  if (the_chrom == chrom_lookup.end())
    throw SMITHLABException("no chrom with id: " + toa(ba.RefID));
  
  string chrom = the_chrom->second;
  size_t start = ba.Position;
  size_t end = ba.GetEndPosition();
  string name(ba.Name);
  uint8_t mismatch = 0;
  char strand = (ba.IsReverseStrand() ? '-' : '+');
  string seq = ba.QueryBases;
  string qual = ba.Qualities;

  // if a read is mapped to - strand, bismark stores the + strand seq of
  // reference genome rather than the seq of that read, however bsmap
  // stores the original read sequence. I'm not sure about the orientation of
  // CIGAR string in bismark. But we do need the original sequence in .mr
  if (strand == '-' && mapper.compare("bismark") == 0) {
    revcomp_inplace(seq);
    std::reverse(qual.begin(), qual.end());
  }
  // bismark adds /1 and /2 in mapping
  if (ba.IsPaired() && mapper.compare("bsmap") == 0) {
    if (ba.IsFirstMate()) {
      name = name + "/1";
    }
    else {
      name = name + "/2";
    }
  }

  string new_seq, new_qual;
  apply_CIGAR_for_BAM(seq, qual, ba.CigarData, new_seq, new_qual);

  if (mapper.compare("bsmap") == 0) {
    if (!ba.GetTag("NM",mismatch))
      throw SMITHLABException(ba.GetErrorString());
  }
  else if (mapper.compare("bismark") == 0) {
    int convert_count = 0;
    string convert_str;
    ba.GetTag("NM",mismatch);
    ba.GetTag("XM",convert_str);

    const char *temp = convert_str.c_str();
    while(*temp != '\0') {
      if (*temp == 'x' || *temp == 'h' || *temp == 'z')
        ++convert_count;
      ++temp;
    }
    mismatch = mismatch - convert_count;
  }
  
  mr.r = GenomicRegion(chrom, start, end, name, mismatch, strand);
  mr.seq = new_seq;
  mr.scr = new_qual;
}

static MappedRead
BamAlignmentToMappedRead(const unordered_map<size_t, string> &chrom_lookup,
			 const BamAlignment &ba, MappedRead &mr) {

  const unordered_map<size_t, string>::const_iterator 
    the_chrom(chrom_lookup.find(ba.RefID));
  if (the_chrom == chrom_lookup.end())
    throw SMITHLABException("no chrom with id: " + toa(ba.RefID));
  
  const string chrom = the_chrom->second;
  const size_t start = ba.Position;
  const size_t end = start + ba.Length;
  const string name(ba.Name);
  const float score = ba.MapQuality;
  const char strand = (ba.IsReverseStrand() ? '-' : '+');
  const string seq = ba.QueryBases;
  const string scr = ba.Qualities;
  mr.r = GenomicRegion(chrom, start, end, name, score, strand);
  mr.seq = seq;
  mr.scr = scr;

  return mr;
}

void
ReadBAMFormatInput(const string &infile, vector<MappedRead> &reads) {
  
  BamReader reader;
  reader.Open(infile);
  
  // Get header and reference
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  
  unordered_map<size_t, string> chrom_lookup;
  for (size_t i = 0; i < refs.size(); ++i)
    chrom_lookup[i] = refs[i].RefName;
  
  BamAlignment bam;
  MappedRead read;
  while (reader.GetNextAlignment(bam))
    reads.push_back(BamAlignmentToMappedRead(chrom_lookup, bam, read));
  reader.Close();
}
