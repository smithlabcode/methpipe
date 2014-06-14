/* Copyright (C) 2012 University of Southern California
 *                    Andrew D Smith
 * Author: Andrew D. Smith, Song Qiang
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

#define DEBUG

#include <sys/types.h>
#include <unistd.h>

#include <numeric>
#include <cmath>
#include <ctime>
#include <fstream>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "ThreeStateHDHMM.hpp"
#include "Distro.hpp"
#include "false_discovery_rate.hpp"
#include "contingency-table.hpp"
#include "RNG.hpp"
#include "nonparametric-test.hpp"
#include "ModelParams.hpp"

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

static const std::string STATE_LABEL_STRS_ARRAY[] = {"GAIN", "SAME", "LOSS"};
static const std::vector<std::string> STATE_LABEL_STRS(
    STATE_LABEL_STRS_ARRAY, STATE_LABEL_STRS_ARRAY + 3);

static void
get_meth_unmeth(const GenomicRegion &cpg, size_t &meth, size_t &unmeth) 
{
    const double prob = cpg.get_score();
    const string name(cpg.get_name());
    const size_t n_reads = atoi(name.substr(name.find_first_of(":") + 1).c_str());
    meth = static_cast<size_t>(prob * n_reads);
    unmeth = n_reads - meth;
}

static void
load_cpgs(const string &cpgs_file_a, const string &cpgs_file_b,
          vector<SimpleGenomicRegion> &cpgs,
          vector<size_t> &meth_a, vector<size_t> &unmeth_a,
          vector<size_t> &meth_b, vector<size_t> &unmeth_b,
          vector<double> &diffscores, const bool VERBOSE)
{
    static const size_t pseudocount = 1;
    
    vector<GenomicRegion> cpgs_a, cpgs_b;
    if (VERBOSE)
        cerr << "Reading input file " << cpgs_file_a << " ... ";
    ReadBEDFile(cpgs_file_a, cpgs_a);
    assert(check_sorted(cpgs_a));
    if (VERBOSE)
        cerr << " Done" << endl;

    if (VERBOSE)
        cerr << "Reading input file " << cpgs_file_b << " ... ";
    ReadBEDFile(cpgs_file_b, cpgs_b);
    assert(check_sorted(cpgs_b));
    if (VERBOSE)
        cerr << " Done" << endl;

    if (VERBOSE)
        cerr << "Calculating diffscores: ";
    size_t j = 0;
    for (size_t i = 0; i < cpgs_a.size(); ++i) 
    {
        while (j < cpgs_b.size() && cpgs_b[j] < cpgs_a[i]) ++j;
      
        if (cpgs_a[i].same_chrom(cpgs_b[j]) && 
            cpgs_a[i].get_start() == cpgs_b[j].get_start()) 
        {
            size_t m_a = 0, u_a = 0;
            get_meth_unmeth(cpgs_a[i], m_a, u_a);
            size_t m_b = 0, u_b = 0;
            get_meth_unmeth(cpgs_b[j], m_b, u_b);
	
            if (m_a + u_a > 0.0 && m_b + u_b > 0.0) 
            {
                const double diffscore = 
                    ContingencyTable::beta_population_greater(
                        m_b + pseudocount, u_b + pseudocount, 
                        m_a + pseudocount, u_a + pseudocount);

                cpgs.push_back(cpgs_a[i]);
                meth_a.push_back(m_a);
                unmeth_a.push_back(u_a);
                meth_b.push_back(m_b);
                unmeth_b.push_back(u_b);
                diffscores.push_back(diffscore);
            }
        }
    }
    if (VERBOSE)
        cerr << "Probes retained: " << diffscores.size() << endl;
}

template <class T> static void
separate_regions(const bool VERBOSE, const size_t desert_size, 
                 vector<SimpleGenomicRegion> &cpgs,
                 vector<T> &diffscores,
                 vector<size_t> &reset_points) 
{
    if (VERBOSE)
        cerr << "[SEPARATING BY CPG DESERT]" << endl;
  
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
pick_sample(const vector<double> &diffscores,
            const vector<size_t> &reset_points,
            const size_t training_size,
            vector<double> &diffscores_sample,
            vector<size_t> &reset_points_sample)
{
    // random training sample
    vector<size_t> idxs(reset_points.size() - 1);
    for (size_t i = 0; i < idxs.size(); ++i) idxs[i] = i;
    srand(time(0) + getpid());
    std::random_shuffle(idxs.begin(), idxs.end());

    size_t sample_size = 0;
    size_t i = 0;
    reset_points_sample.push_back(sample_size);
    while (i < idxs.size() && sample_size < training_size)
    {	
        const size_t idx = idxs[i];
        const size_t start = reset_points[idx];
        const size_t end = reset_points[idx + 1];

        std::copy(diffscores.begin() + start, diffscores.begin() + end,
                  std::back_inserter(diffscores_sample));
        sample_size += end - start;
        reset_points_sample.push_back(sample_size);
        ++i;
    }

    assert(diffscores_sample.size() <= diffscores.size());
    assert(diffscores_sample.size() == reset_points_sample.back());
}

static inline double
max(const double &a, const double &b, const double &c)
{
    return max(max(a, b), c);
}

static void
build_domains(const bool VERBOSE, 
              const vector<SimpleGenomicRegion> &cpgs,
              const vector<size_t> &meth_a,
              const vector<size_t> &unmeth_a,
              const vector<size_t> &meth_b,
              const vector<size_t> &unmeth_b,
              const vector<Triplet> &post_scores,
              const vector<size_t> &reset_points,
              const vector<STATE_LABELS> &classes,
              vector<GenomicRegion> &domains,
              vector<double> &domain_posterior_scores) 
{
    assert(cpgs.size() == post_scores.size());
    assert(cpgs.size() == classes.size());
    assert(cpgs.size() == reset_points.back());

    vector<double> diff_meths(cpgs.size());
    for (size_t i = 0; i < cpgs.size(); ++i)
        diff_meths[i] = static_cast<double>(meth_b[i]) / (meth_b[i]+unmeth_b[i])
            - static_cast<double>(meth_a[i]) / (meth_a[i]+unmeth_a[i]);

    for (size_t idx = 0; idx < reset_points.size() - 1; ++idx)
    {
        const size_t start  = reset_points[idx];
        const size_t end = reset_points[idx + 1];
        GenomicRegion d(cpgs[start]);
        size_t n_cpgs = 1;
        double diff_meth = diff_meths[start];
        double posterior_sum =
            max(post_scores[start].gain, post_scores[start].same,
                post_scores[start].loss);
        for (size_t i = start + 1; i < end; ++i)
        {
            if (classes[i] == classes[i - 1])
            {
                diff_meth += diff_meths[i];
                ++n_cpgs;
                posterior_sum +=
                    max(post_scores[i].gain, post_scores[i].same,
                        post_scores[i].loss);
            }
            else
            {
                d.set_end(cpgs[i - 1].get_end());
                d.set_score(diff_meth);
                d.set_name(STATE_LABEL_STRS[classes[i - 1]]
                           + ":" + smithlab::toa(n_cpgs));
                if (classes[i - 1] != SAME)
                {
                    domains.push_back(d);
                    domain_posterior_scores.push_back(posterior_sum);
                }
                
                d = GenomicRegion(cpgs[i]);
                n_cpgs = 1;
                diff_meth = diff_meths[i];
                posterior_sum =
                    max(post_scores[i].gain, post_scores[i].same,
                        post_scores[i].loss);
            }
        }
        
        d.set_end(cpgs[end - 1].get_end());
        d.set_score(diff_meth);
        d.set_name(STATE_LABEL_STRS[classes[end - 1]]
                   + ":" + smithlab::toa(n_cpgs));
        if (classes[end - 1] != SAME) 
        {
            domains.push_back(d);
            domain_posterior_scores.push_back(posterior_sum);
        }
    }
}

static void
calcualte_domain_p_values_by_wilcoxon_test(
    const vector<SimpleGenomicRegion> &cpgs,
    const vector<size_t> &meth_a, const vector<size_t> &unmeth_a,
    const vector<size_t> &meth_b, const vector<size_t> &unmeth_b,
    const vector<GenomicRegion> &domains, vector<double> &p_values)
{
    assert(check_sorted(cpgs));
    assert(check_sorted(domains));

    size_t j = 0;
    for (size_t i = 0; i < domains.size(); ++i)
    {
        const SimpleGenomicRegion simdom(domains[i]);
        vector<double> methylation_a, methylation_b;
        size_t domain_size = 0;
        
        while (j < cpgs.size() && !simdom.contains(cpgs[j])) ++j;

        while (j < cpgs.size() && simdom.contains(cpgs[j]))
        {	
            methylation_a.push_back(
                static_cast<double>(meth_a[j]) / (meth_a[j] + unmeth_a[j]));
            methylation_b.push_back(
                static_cast<double>(meth_b[j]) / (meth_b[j] + unmeth_b[j]));
            ++domain_size;
            ++j;
        }

        const double p_value =
            (domains[i].get_name().substr(0, 4)
             == STATE_LABEL_STRS[GAIN])
            ? NonParametricTest::wilcoxon_test(methylation_b, methylation_a)
            : NonParametricTest::wilcoxon_test(methylation_a, methylation_b);
        assert(p_value >= 0 && p_value <= 1);
        p_values.push_back(p_value);
    }

    assert(p_values.size() == domains.size());
}

static void
score_domain_by_wilcoxon_test(
    const vector<SimpleGenomicRegion> &cpgs,
    const vector<size_t> &meth_a,  const vector<size_t> &unmeth_a, 
    const vector<size_t> &meth_b,  const vector<size_t> &unmeth_b, 
    const double fdr, double fdr_cutoff, vector<GenomicRegion> &domains,
    const bool VERBOSE)
{
    if (VERBOSE)
        cerr << "Computing FDR cutoff with Wilcoxon signed-rank test" << endl;

    vector<double> p_values;
    calcualte_domain_p_values_by_wilcoxon_test(
        cpgs, meth_a, unmeth_a, meth_b, unmeth_b, domains, p_values);
    if (fdr_cutoff == numeric_limits<double>::max())
        fdr_cutoff = FDR::get_fdr_cutoff(p_values, fdr);
    if (VERBOSE)
        cerr << "cutoff = " << fdr_cutoff << endl;
    
    // filtering domains
    size_t j = 0;
    for (size_t i = 0; i < domains.size(); ++i)
        if (p_values[i] <= fdr_cutoff)
        {
            domains[i].set_name(domains[i].get_name() + ":"
                                + smithlab::toa(p_values[i]));
            domains[j] = domains[i];
            ++j;
        }
    domains.erase(domains.begin() + j, domains.end());
}


static void
calculate_random_scores_from_background(
    const vector<size_t> &meth_a,
    const vector<size_t> &unmeth_a,
    const vector<size_t> &meth_b,
    const vector<size_t> &unmeth_b,
    const size_t & cpg_num,
    const size_t & times, 
    vector<double> &random_scores) 
{
    assert(meth_a.size() == unmeth_a.size()
           && meth_a.size() == meth_b.size()
           && meth_a.size() == unmeth_b.size());
    
    vector<double> diffmeth(meth_a.size());
    for (size_t i = 0; i < meth_a.size(); ++i)
    {
        diffmeth[i] =
            static_cast<double>(meth_b[i]) / (meth_b[i] + unmeth_b[i]) 
            - static_cast<double>(meth_a[i]) / (meth_a[i] + unmeth_a[i]);
    }

    Runif rng(std::time(NULL) + getpid());
    for (size_t j = 0; j < times; ++j)
    {
        const size_t start = rng.runif(static_cast<size_t>(0),
                                       meth_a.size() - cpg_num);
        const double sum_diffmeth =
            std::accumulate(diffmeth.begin() + start,
                            diffmeth.begin() + start + cpg_num, 0.0);
        random_scores.push_back(fabs(sum_diffmeth));
    }
    
    std::sort(random_scores.begin(), random_scores.end());
}

class DomainSizeCmp:
    public std::binary_function<GenomicRegion, GenomicRegion, bool>
{
public:
    bool operator()(const GenomicRegion &lhs, const GenomicRegion &rhs) const
    {
        return atoi(lhs.get_name().substr(5).c_str())
            < atoi(rhs.get_name().substr(5).c_str());
    }
};

static void
score_domain_by_diff_meth_emp_p_value(
    const vector<SimpleGenomicRegion> &cpgs,
    const vector<size_t> &meth_a,  const vector<size_t> &unmeth_a, 
    const vector<size_t> &meth_b,  const vector<size_t> &unmeth_b, 
    const double fdr, double fdr_cutoff, vector<GenomicRegion> &domains,
    const bool VERBOSE)
{
    if (VERBOSE)
        cerr << "Computing FDR cutoff ... ";

    std::sort(domains.begin(), domains.end(), DomainSizeCmp());

    vector<double> p_values(domains.size());

    size_t domain_cpg_size = 0;
    vector<double> random_scores;
    for (size_t i = 0; i < domains.size(); ++i)
    {
        const size_t sz = atoi(domains[i].get_name().substr(5).c_str());
        const double domain_score = fabs(domains[i].get_score());

        assert(domain_cpg_size <= sz);
        if (domain_cpg_size < sz)
        {
            static const size_t random_scores_size = 50000;
            domain_cpg_size = sz;
            random_scores.clear();
            calculate_random_scores_from_background(
                meth_a, unmeth_a, meth_b, unmeth_b, domain_cpg_size,
                random_scores_size, random_scores);
        }
        
        p_values[i] = FDR::get_empirical_p_value(random_scores, domain_score);
    }

    if (fdr_cutoff == numeric_limits<double>::max())
        fdr_cutoff = FDR::get_fdr_cutoff(p_values, fdr);
    if (VERBOSE)
        cerr << "cutoff = " << fdr_cutoff << endl;

    // filtering domains
    size_t j = 0;
    for (size_t i = 0; i < domains.size(); ++i)
        if (p_values[i] <= fdr_cutoff)
        {
            domains[i].set_name(domains[i].get_name() + ":"
                                + smithlab::toa(p_values[i]));
            domains[j] = domains[i];
            ++j;
        }
    domains.erase(domains.begin() + j, domains.end());
    std::sort(domains.begin(), domains.end());
}

static void
write_scores_file(const string &scores_file,
                  const vector<SimpleGenomicRegion> &cpgs,
                  const vector<double> &diffscores,
                  const vector<Triplet> &posterior_scores)
{
    assert(cpgs.size() == diffscores.size());
    assert(cpgs.size() == posterior_scores.size());

    std::ofstream outf(scores_file.c_str());
    for (size_t i = 0; i < cpgs.size(); ++i) 
        outf << cpgs[i] << "\t" << diffscores[i] << "\t"
             << posterior_scores[i].gain << "\t"
             << posterior_scores[i].same << "\t"
             << posterior_scores[i].loss << endl;
    outf.close();
}

int
main(int argc, const char **argv) 
{
    try 
    {
        string outfile("/dev/stdout");
        string infile_a, infile_b;
        string scores_file;
        string params_in_file;
        string params_out_file;

        size_t desert_size = 1000;
        size_t max_iterations = 10;
        size_t training_size = 0;
        
        bool diff_meth_emp_p_value = false;
        bool domain_wilcoxon_test = false;
        double fdr = 0.05;
        double fdr_cutoff = std::numeric_limits<double>::max();
        
        // corrections for small values (not parameters):
        double tolerance = 1e-10;
        size_t MAX_LEN = 200;

        // run mode flags
        bool VERBOSE = false;
    
        /****************** COMMAND LINE OPTIONS ********************/
        OptionParser opt_parse(argv[0], "A program for segmenting DNA "
                               "methylation data");
        opt_parse.add_opt("input-a", 'A', "Methcount input file A", 
                          OptionParser::REQUIRED, infile_a);
        opt_parse.add_opt("input-b", 'B', "Methcount input file B", 
                          OptionParser::REQUIRED, infile_b);
        opt_parse.add_opt("out", 'o', "output file (default stdout)", 
                          OptionParser::OPTIONAL, outfile);
        opt_parse.add_opt("scores", 's', "scores file (WIG format)", 
                          OptionParser::OPTIONAL, scores_file);
        opt_parse.add_opt("params-in", '\0', "HMM parameters file",
                          OptionParser::OPTIONAL, params_in_file);
        opt_parse.add_opt("params-out", '\0', "HMM parameters file",
                          OptionParser::OPTIONAL, params_out_file);
        opt_parse.add_opt("desert", 'd', "desert size",
                          OptionParser::OPTIONAL, desert_size);
        opt_parse.add_opt("itr", 'i', "max iterations",
                          OptionParser::OPTIONAL, max_iterations); 
        opt_parse.add_opt("training-size", '\0', "The size of training sample",
                          OptionParser::OPTIONAL, training_size); 
        opt_parse.add_opt("max-len", 'L', "max foreground length",
                          OptionParser::OPTIONAL, MAX_LEN); 
        opt_parse.add_opt("domain-wilcoxon-test", '\0',
                          "Use Wilcoxon test on domains to compute p-values",
                          OptionParser::OPTIONAL, domain_wilcoxon_test); 
        opt_parse.add_opt("diff-meth", '\0',
                          "Use the differential methylation to compute p-values",
                          OptionParser::OPTIONAL, diff_meth_emp_p_value); 
        opt_parse.add_opt("fdr", 'F', "False discovery rate (default 0.05)",
                          OptionParser::OPTIONAL, fdr); 
        opt_parse.add_opt("fdr-cutoff", '\0',
                          "P-value cutoff based on false discovery rate",
                          OptionParser::OPTIONAL, fdr_cutoff); 
        opt_parse.add_opt("verbose", 'v', "print more run info", 
                          OptionParser::OPTIONAL, VERBOSE);

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

        /****************** END COMMAND LINE OPTIONS *****************/
    
        /***********************************
         * STEP 1: READ IN INPUT
         */
        // separate the regions by chrom and by desert
        vector<SimpleGenomicRegion> cpgs;
        vector<double> diffscores;
        vector<size_t> meth_a, unmeth_a, meth_b, unmeth_b;
        load_cpgs(infile_a, infile_b, cpgs,
                  meth_a, unmeth_a, meth_b, unmeth_b, diffscores, VERBOSE);
            
        // separate the regions by chrom and by desert, and eliminate
        // those isolated CpGs
        vector<size_t> reset_points;
        separate_regions(VERBOSE, desert_size, cpgs, diffscores, reset_points);
        ThreeStateHDHMM hmm(diffscores, reset_points, tolerance,
                            max_iterations, VERBOSE, MAX_LEN);
    

        /***********************************
         * STEP 1: MODEL INITIALIZATION AND TRAINING
         */
        size_t state_num = 3;
        vector<Distro> emissions;
        vector<Distro> durations;
        vector<vector<double> > trans(3, vector<double>(3, 0));
        if (!params_in_file.empty()) 
        {
            read_param_file(params_in_file, state_num,
                            trans, emissions, durations);
        }
        else 
        {
            emissions.push_back(Distro("beta 0.5 2"));
            emissions.push_back(Distro("beta 1.5 1.5"));
            emissions.push_back(Distro("beta 2 0.5"));

            durations.push_back(Distro("nbd 10 0.5"));
            durations.push_back(Distro("geo 0.02"));
            durations.push_back(Distro("nbd 10 0.5"));

            trans[GAIN][SAME] = trans[LOSS][SAME] = 1;
            trans[SAME][SAME] = 1 - durations[1].get_params().front();
            trans[SAME][GAIN] = trans[SAME][LOSS] = (1-trans[SAME][SAME]) * 0.5;
        }

        // model training
        if (max_iterations > 0)
        {
            if (training_size != 0)
            {
                // train with part of the dataset
                vector<double>  diffscores_sample;
                vector<size_t> reset_points_sample;
                pick_sample(diffscores, reset_points, training_size,
                            diffscores_sample, reset_points_sample);
                ThreeStateHDHMM hmm_training(
                    diffscores_sample, reset_points_sample,
                    tolerance, max_iterations, VERBOSE, MAX_LEN);
                
                hmm_training.set_parameters(
                    emissions.front(), emissions[1], emissions.back(),
                    durations.front(), durations[1], durations.back(),
                    trans);
                hmm_training.BaumWelchTraining();
                hmm_training.get_parameters(
                    emissions.front(), emissions[1], emissions.back(),
                    durations.front(), durations[1], durations.back(),
                    trans);
            }
            else
            {
                // train with the whole dataset
                hmm.set_parameters(
                    emissions.front(), emissions[1], emissions.back(),
                    durations.front(), durations[1], durations.back(),
                    trans);
                hmm.BaumWelchTraining();
            }
        }

        hmm.set_parameters(
            emissions.front(), emissions[1], emissions.back(),
            durations.front(), durations[1], durations.back(),
            trans);

        if (!params_out_file.empty()) 
        {
            write_param_file(params_out_file, state_num,
                             trans, emissions, durations);
        }
        
        /***********************************
         * STEP 3: DECODE THE DOMAINS
         */
        vector<STATE_LABELS> classes;
        vector<Triplet> scores;
        hmm.PosteriorDecoding();
        hmm.get_posterior_scores(scores, classes);
        vector<GenomicRegion> domains;
        vector<double> domain_posterior_scores;
        build_domains(VERBOSE, cpgs, meth_a, unmeth_a, meth_b, unmeth_b, scores,
                      reset_points, classes, domains, domain_posterior_scores);
        
        if (domain_wilcoxon_test)
            score_domain_by_wilcoxon_test(
                cpgs, meth_a, unmeth_a, meth_b, unmeth_b,
                fdr, fdr_cutoff, domains, VERBOSE);
        
        if (diff_meth_emp_p_value)
            score_domain_by_diff_meth_emp_p_value(
                cpgs, meth_a, unmeth_a, meth_b, unmeth_b,
                fdr, 0.01, domains, VERBOSE);

        /***********************************
         * STEP 6: WRITE THE RESULTS
         */
        if (!scores_file.empty())
        {
            if (VERBOSE)
                cerr << "Writing diffscore and posterior score file: "
                     << scores_file << endl;
            write_scores_file(scores_file, cpgs, diffscores, scores);
        }

        if (VERBOSE)
            cerr << "Writing differentially-methylated regions file: "
                 << outfile << endl;
        std::ofstream out(outfile.c_str());
        std::copy(domains.begin(), domains.end(),
                  std::ostream_iterator<GenomicRegion>(out, "\n"));
        out.close();
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


