/* Copyright (C) 2009 University of Southern California
 *                    Andrew D Smith
 * Author: Song Qiang, Andrew D. Smith
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

#include <numeric>
#include <cmath>
#include <fstream>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "ThreeStateHMM.hpp"
#include "MethpipeFiles.hpp"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;
using std::max;
using std::min;
using std::pair;
using std::make_pair;

using std::ostream_iterator;
using std::ofstream;


static void
load_cpgs(const bool VERBOSE,
          string cpgs_file, vector<SimpleGenomicRegion> &cpgs,
          vector<pair<double, double> > &meth, vector<size_t> &reads)
{
  if (VERBOSE)
    cerr << "[READING CPGS AND METH PROPS]" << endl;
  vector<GenomicRegion> cpgs_in;
  ReadBEDFile(cpgs_file, cpgs_in);
  if (!check_sorted(cpgs_in))
    throw SMITHLABException("CpGs not sorted in file \"" + cpgs_file + "\"");

  for (size_t i = 0; i < cpgs_in.size(); ++i)
    {
      cpgs.push_back(SimpleGenomicRegion(cpgs_in[i]));
      meth.push_back(make_pair(cpgs_in[i].get_score(), 0.0));
      const string r(cpgs_in[i].get_name());
      reads.push_back(atoi(r.substr(r.find_first_of(":") + 1).c_str()));
      meth.back().first = int(meth.back().first * reads.back());
      meth.back().second = int(reads.back() - meth.back().first);
    }
  if (VERBOSE)
    cerr << "TOTAL CPGS: " << cpgs.size() << endl
         << "MEAN COVERAGE: "
         << accumulate(reads.begin(), reads.end(), 0.0)/reads.size() << endl
         << endl;
}

static void
bin_by_basepair(const size_t bin_size,
                vector<SimpleGenomicRegion> &cpgs,
                vector< pair<double, double> > &meth,
                vector<size_t> &reads)
{
  vector<SimpleGenomicRegion> bins;
  vector<pair<double, double> > new_meth;
  vector<size_t> new_reads;

  bins.push_back(SimpleGenomicRegion(cpgs.front().get_chrom(),
                                     cpgs.front().get_start(),
                                     cpgs.front().get_start() + bin_size));
  new_meth.push_back(std::make_pair(0, 0));
  new_reads.push_back(0);
  size_t i = 0;
  while (i < cpgs.size())
    {
      if (bins.back().contains(cpgs[i]))
        {
          new_meth.back().first += meth[i].first;
          new_meth.back().second += meth[i].second;
          ++i;
        }
      else
        {
          new_reads.back() = static_cast<size_t>(
                                                 new_meth.back().first + new_meth.back().second);

          if (bins.back().same_chrom(cpgs[i]))
            bins.push_back(
                           SimpleGenomicRegion(bins.back().get_chrom(),
                                               bins.back().get_end(),
                                               bins.back().get_end() + bin_size));
          else
            bins.push_back(
                           SimpleGenomicRegion(cpgs[i].get_chrom(),
                                               cpgs[i].get_start(),
                                               cpgs[i].get_start() + bin_size));

          new_meth.push_back(std::make_pair(0, 0));
          new_reads.push_back(0);
        }
    }
  new_reads.back() = static_cast<size_t>(
                                         new_meth.back().first + new_meth.back().second);

  std::swap(cpgs, bins);
  std::swap(meth, new_meth);
  std::swap(reads, new_reads);
}

static void
bin_by_loci(const size_t bin_size,
            vector<SimpleGenomicRegion> &cpgs,
            vector< pair<double, double> > &meth,
            vector<size_t> &reads)
{
  vector<SimpleGenomicRegion> bins;
  vector<pair<double, double> > new_meth;
  vector<size_t> new_reads;

  bins.push_back(cpgs.front());
  new_meth.push_back(meth.front());
  new_reads.push_back(reads.front());
  size_t n = 1;

  size_t i = 1;
  while (i < cpgs.size())
    {
      if (bins.back().same_chrom(cpgs[i]) && n < bin_size)
        {
          bins.back().set_end(cpgs[i].get_end());
          new_meth.back().first += meth[i].first;
          new_meth.back().second += meth[i].second;
          ++n;
          ++i;
        }
      else
        {
          new_reads.back() = static_cast<size_t>(
                                                 new_meth.back().first + new_meth.back().second);

          bins.push_back(cpgs[i]);
          new_meth.push_back(meth[i]);
          new_reads.push_back(reads[i]);

          n = 1;
          ++i;
        }
    }
  new_reads.back() = static_cast<size_t>(
                                         new_meth.back().first + new_meth.back().second);

  std::swap(cpgs, bins);
  std::swap(meth, new_meth);
  std::swap(reads, new_reads);
}

template <class T, class U> static void
separate_regions(const bool VERBOSE, const size_t desert_size,
                 vector<SimpleGenomicRegion> &cpgs,
                 vector<T> &meth, vector<U> &reads,
                 vector<size_t> &reset_points)
{
  if (VERBOSE)
    cerr << "[SEPARATING BY CPG DESERT]" << endl;
  // eliminate the zero-read cpgs
  size_t j = 0;
  for (size_t i = 0; i < cpgs.size(); ++i)
    if (reads[i] > 0)
      {
        cpgs[j] = cpgs[i];
        meth[j] = meth[i];
        reads[j] = reads[i];
        ++j;
      }
  cpgs.erase(cpgs.begin() + j, cpgs.end());
  meth.erase(meth.begin() + j, meth.end());
  reads.erase(reads.begin() + j, reads.end());

  // segregate cpgs
  size_t prev_cpg = 0;
  for (size_t i = 0; i < cpgs.size(); ++i)
    {
      const size_t dist = (i > 0 && cpgs[i].same_chrom(cpgs[i - 1])) ?
        cpgs[i].get_start() - prev_cpg : numeric_limits<size_t>::max();
      if (dist > desert_size)
        reset_points.push_back(i);
      prev_cpg = cpgs[i].get_start();
    }
  reset_points.push_back(cpgs.size());
  if (VERBOSE)
    cerr << "CPGS RETAINED: " << cpgs.size() << endl
         << "DESERTS REMOVED: " << reset_points.size() - 2 << endl << endl;
}


static void
read_params_file(const string &params_file,
                 betabin &hypo_emission,
                 betabin &HYPER_emission,
                 betabin &HYPO_emission,
                 vector<vector<double> > &trans)
{
  std::ifstream in(params_file.c_str());
  string hypo_emission_str, HYPER_emission_str, HYPO_emission_str;
  std::getline(in, hypo_emission_str);
  std::getline(in, HYPER_emission_str);
  std::getline(in, HYPO_emission_str);

  trans.resize(3, vector<double>(3, 0.0));
  for (size_t i = 0; i < trans.size(); ++i)
    for (size_t j = 0; j < trans[i].size(); ++j)
      in >> trans[i][j];
  in.close();
}

static void
write_params_file(const string &params_file,
                  const betabin &hypo_emission,
                  const betabin &HYPER_emission,
                  const betabin &HYPO_emission,
                  const vector<vector<double> > &trans)
{
  std::ofstream out(params_file.c_str());
  out << hypo_emission.tostring() << endl
      << HYPER_emission.tostring() << endl
      << HYPO_emission.tostring() << endl;

  for (size_t i = 0; i < trans.size(); ++i)
    {
      for (size_t j = 0; j < trans[i].size(); ++j)
        out << trans[i][j] << "\t";
      out << endl;
    }
  out.close();
}


static void
build_domains(const bool VERBOSE,
              const vector<SimpleGenomicRegion> &cpgs,
              const vector<pair<double, double> > &meth,
              const vector<size_t> &reset_points,
              const vector<STATE_LABELS> &classes,
              vector<GenomicRegion> &domains)
{
  domains.clear();

  for (size_t i = 0; i < reset_points.size() - 1; ++i)
    {
      const size_t start = reset_points[i];
      const size_t end = reset_points[i + 1];

      GenomicRegion domain(cpgs[start]);
      STATE_LABELS prev_state = classes[start];
      size_t n = 1;
      string hmrcpgs = cpgs[start].tostring();
      double meth_sum =
        meth[start].first / (meth[start].first + meth[start].second);

      for (size_t j = start + 1; j < end; ++j)
        {
          // if ((prev_state == hypo && classes[j] == HYPO)
          //     || (prev_state == HYPO && classes[j] == hypo))
          //     cerr << "WARNING: inconsist state sequences"
          //         " from posterior decoding" << endl;

          if ((prev_state == hypo && classes[j] == hypo)
              || (prev_state != hypo && classes[j] != hypo))
            {
              ++n;
              meth_sum += meth[j].first / (meth[j].first + meth[j].second);
              hmrcpgs += ":" + cpgs[j].tostring();
            }
          else
            {
              domain.set_end(cpgs[j - 1].get_end());
              switch (prev_state) {
              case hypo:
                domain.set_name("hypo:" + smithlab::toa(n));
                break;
              case HYPER:
                domain.set_name("hyper:" + smithlab::toa(n));
                break;
              case HYPO:
                domain.set_name("hyper:" + smithlab::toa(n));
                break;
              }
              domain.set_score(meth_sum);
              domain.set_strand('+');
              if (prev_state == HYPER || prev_state == HYPO)
                domains.push_back(domain);

              domain = GenomicRegion(cpgs[j]);
              n = 1;
              hmrcpgs = cpgs[j].tostring();
              prev_state = classes[j];
              meth_sum = meth[j].first / (meth[j].first + meth[j].second);
            }
        }
      domain.set_end(cpgs[end - 1].get_end());
      switch (prev_state)
        {
        case hypo: domain.set_name("hypo:" + smithlab::toa(n)); break;
        case HYPER: domain.set_name("hyper:" + smithlab::toa(n)); break;
        case HYPO: domain.set_name("hyper:" + smithlab::toa(n)); break;
        }
      domain.set_score(meth_sum);
      domain.set_strand('+');
      if (prev_state == HYPER || prev_state == HYPO)
        domains.push_back(domain);
    }
}

static void
filter_domains(const bool VERBOSE, const double MIN_ACCUMULATIVE_METH,
               vector<GenomicRegion> &domains)
{
  size_t j = 0;
  for (size_t i = 0; i < domains.size(); ++i)
    if (domains[i].get_score() >= MIN_ACCUMULATIVE_METH)
      {
        domains[j] = domains[i];
        ++j;
      }
  domains.erase(domains.begin() + j, domains.end());
}

int
main(int argc, const char **argv)
{
  try
    {
      string outfile("/dev/stdout");
      string scores_file;
      string trans_file;

      size_t desert_size = 1000;
      size_t max_iterations = 10;

      // run mode flags
      bool VERBOSE = false;

      // corrections for small values (not parameters):
      double tolerance = 1e-10;

      double MIN_ACCUMULATIVE_METH = 4.0;
      size_t bin_size = 0;
      double BIN_BY_LOCI = false;
      bool USE_VITERBI_DECODING = false;

      string params_in_file;
      string params_out_file;

      /****************** COMMAND LINE OPTIONS ********************/
      OptionParser opt_parse(argv[0], "A program for segmenting DNA "
                             "methylation data"
                             "<cpg-BED-file>");
      opt_parse.add_opt("out", 'o', "output file (BED format)",
                        false, outfile);
      opt_parse.add_opt("scores", 's', "output file for posterior scores",
                        OptionParser::OPTIONAL, scores_file);
      opt_parse.add_opt("tolerance", 't', "Tolerance",
                        false, tolerance);
      opt_parse.add_opt("bin_size", 'b', "bin size (default: single base)",
                        OptionParser::OPTIONAL, bin_size);
      opt_parse.add_opt("bin_by_loci", '\0',
                        "binning genome by number "
                        "of loci (default by basepair)",
                        OptionParser::OPTIONAL, BIN_BY_LOCI);
      opt_parse.add_opt("desert", 'd', "desert size", false, desert_size);
      opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations);
      opt_parse.add_opt("viterbi", 'V', "Use Viterbi decoding",
                        OptionParser::OPTIONAL, USE_VITERBI_DECODING);
      opt_parse.add_opt("min-meth", 'M',
                        "Min accumulative methylation levels within HypeMR",
                        OptionParser::OPTIONAL, MIN_ACCUMULATIVE_METH);
      opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
      opt_parse.add_opt("params-in", 'P', "HMM parameters file",
                        false, params_in_file);
      opt_parse.add_opt("params-out", 'p', "HMM parameters file",
                        false, params_out_file);

      vector<string> leftover_args;
      opt_parse.parse(argc, argv, leftover_args);
      if (argc == 1 || opt_parse.help_requested())
        {
          cerr << opt_parse.help_message() << endl;
          return EXIT_SUCCESS;
        }
      if (opt_parse.about_requested())
        {
          cerr << opt_parse.about_message() << endl;
          return EXIT_SUCCESS;
        }
      if (opt_parse.option_missing())
        {
          cerr << opt_parse.option_missing_message() << endl;
          return EXIT_SUCCESS;
        }
      if (leftover_args.empty())
        {
          cerr << opt_parse.help_message() << endl;
          return EXIT_SUCCESS;
        }
      const string cpgs_file = leftover_args.front();
      /****************** END COMMAND LINE OPTIONS *****************/

      // separate the regions by chrom and by desert
      vector<SimpleGenomicRegion> cpgs;
      // vector<double> meth;
      vector<pair<double, double> > meth;
      vector<size_t> reads;
      if (methpipe::is_methpipe_file_single(cpgs_file))
        {
          if (VERBOSE)
            cerr << "[READING CPGS AND METH PROPS]" << endl;
          methpipe::load_cpgs(cpgs_file, cpgs, meth, reads);
          if (VERBOSE)
            cerr << "TOTAL CPGS: " << cpgs.size() << endl
                 << "MEAN COVERAGE: "
                 << accumulate(reads.begin(), reads.end(), 0.0)/reads.size()
                 << endl << endl;
        }
      else
        load_cpgs(VERBOSE, cpgs_file, cpgs, meth, reads);

      // if bin_size is greater than single base, do bining first
      // applicable to low coverage sample
      if (bin_size > 0)
        {
          if (BIN_BY_LOCI)
            bin_by_loci(bin_size, cpgs, meth, reads);
          else
            bin_by_basepair(bin_size, cpgs, meth, reads);
        }

      // separate the regions by chrom and by desert, and eliminate
      // those isolated CpGs
      vector<size_t> reset_points;
      separate_regions(VERBOSE, desert_size, cpgs, meth, reads, reset_points);

      ThreeStateHMM hmm(meth, reset_points, tolerance, max_iterations, VERBOSE);

      betabin hypo_emission, HYPER_emission, HYPO_emission;
      vector<vector<double> > trans(3, vector<double>(3, 0.0));
      trans[hypo][hypo] = 0.99;
      trans[hypo][HYPER] = 1 - 0.99;
      trans[HYPER][hypo] = 1 - 0.95;
      trans[HYPER][HYPER] = 0.95 * 0.5;
      trans[HYPER][HYPO] = 0.95 * 0.5;
      trans[HYPO][HYPER] = 0.6;
      trans[HYPO][HYPO] = 0.3;

      // double fdr_cutoff = std::numeric_limits<double>::max();

      if (!params_in_file.empty())
        {
          read_params_file(params_in_file,
                           hypo_emission, HYPER_emission, HYPO_emission,
                           trans);
        }
      else
        {
          const double n_reads =
            accumulate(reads.begin(), reads.end(), 0.0)/reads.size();
          double fg_alpha = 0.33*n_reads;
          double fg_beta = 0.67*n_reads;
          double bg_alpha = 0.67*n_reads;
          double bg_beta = 0.33*n_reads;

          hypo_emission = betabin(fg_alpha, fg_beta);
          HYPER_emission = betabin(bg_alpha, bg_beta);
          HYPO_emission = hypo_emission;
        }

      hmm.set_parameters(hypo_emission, HYPER_emission, HYPO_emission, trans);
      if (max_iterations > 0) hmm.BaumWelchTraining();
      hmm.get_parameters(hypo_emission, HYPER_emission, HYPO_emission, trans);

      if (!params_out_file.empty())
        write_params_file(params_out_file, hypo_emission, HYPER_emission,
                          HYPO_emission, trans);

      /***********************************
       * STEP 5: DECODE THE DOMAINS
       */
      vector<STATE_LABELS> classes;
      if (USE_VITERBI_DECODING)
        hmm.ViterbiDecoding();
      else
        hmm.PosteriorDecoding();
      hmm.get_classes(classes);

      /***********************************
       * STEP 6: WRITE THE RESULTS
       */
      vector<GenomicRegion> domains;
      build_domains(VERBOSE, cpgs, meth, reset_points, classes, domains);
      filter_domains(VERBOSE, MIN_ACCUMULATIVE_METH, domains);

      std::ofstream out(outfile.c_str());

      for (size_t i = 0; i < domains.size(); ++i)
        {
          out << domains[i] << endl;
        }
      out.close();

      if (!scores_file.empty())
        {
          if (USE_VITERBI_DECODING) hmm.PosteriorDecoding();
          vector<Triplet> scores;
          hmm.get_state_posteriors(scores);
          std::ofstream score_out(scores_file.c_str());
          for (size_t i = 0; i < cpgs.size(); ++i)
            {
              score_out << cpgs[i] << "\t";
              switch (classes[i])
                {
                case hypo:
                  score_out << "hypo" << "\t" << scores[i].hypo << endl;
                  break;
                case HYPER:
                  score_out << "HYPER" << "\t" << scores[i].HYPER << endl;
                  break;
                case HYPO:
                  score_out << "HYPO" << "\t" << scores[i].HYPO << endl;
                  break;
                }
            }
        }
    }
  catch (SMITHLABException &e)
    {
      cerr << "ERROR:\t" << e.what() << endl;
      return EXIT_FAILURE;
    }
  catch (std::bad_alloc &ba)
    {
      cerr << "ERROR: could not allocate memory" << endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
