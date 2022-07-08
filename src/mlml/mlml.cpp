/*   Copyright (C) 2014 University of Southern California and
 *                      Meng Zhou, Jenny Qu and Andrew D. Smith
 *
 *   Authors: Meng Zhou, Jenny Qu and Andrew D. Smith
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 */

#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <cassert>
#include <sstream>
#include <cmath>

#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "MethpipeSite.hpp"

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::max;
using std::min;
using std::round;
using std::ofstream;
using std::ifstream;
using std::runtime_error;

static void
wilson_ci_for_binomial(const double alpha, const double n,
                       const double p_hat, double &lower, double &upper) {
  const double z = gsl_cdf_ugaussian_Pinv(1 - alpha/2);
  const double denom = 1 + z*z/n;
  const double first_term = p_hat + z*z/(2*n);
  const double discriminant = p_hat*(1 - p_hat)/n + z*z/(4*n*n);
  lower = std::max(0.0, (first_term - z*std::sqrt(discriminant))/denom);
  upper = std::min(1.0, (first_term + z*std::sqrt(discriminant))/denom);
}

static int
binom_null(const double alpha, const double n,
           const double p_hat, const double p) {
  double lower;
  double upper;
  wilson_ci_for_binomial(alpha, n, p_hat, lower, upper);
  if (p < upper && p >lower) return 0;
  else return 1;
}

/* NOTATION:
 * p_m: probability of mC
 * p_h: probability of hmC
 *
 * TAB-seq VARIABLES (for 5hmC)
 * h: counts of C-reads
 * g: counts of T-reads
 *
 * oxBS-seq VARIABLES (for 5mC)
 * m: counts of C-reads
 * l: counts of T-reads
 *
 * BS-seq VARIABLES (for both 5hmC and 5mC)
 * t: counts of C-reads
 * u: counts of T-reads
 *
 * LATENT VARIABLES:
 * t1 (k in expectation): Cs that are from 5mC in BS-seq
 * g1 (j in expectation): Ts that are from 5mC in TAB-seq
 */


//////////////////////////
/// All 3 input files ////
//////////////////////////
static double
log_L(const size_t h, const size_t g, const size_t m, const size_t l,
      const size_t u, const size_t t, const double p_h, const double p_m) {
  double log_lkhd = gsl_sf_lnchoose(h+g, h) + gsl_sf_lnchoose(m+l, m) +
    gsl_sf_lnchoose(u+t, u);

  if (p_h > 0) log_lkhd += h*log(p_h);
  if (p_h < 1) log_lkhd += g*log(1-p_h);

  if (p_m > 0) log_lkhd += m*log(p_m);
  if (p_m < 1) log_lkhd += l*log(1-p_m);

  if (p_h + p_m < 1) log_lkhd += u*log(1-p_h-p_m);
  if (p_h + p_m > 0)log_lkhd += t*log(p_h+p_m);

  return log_lkhd;
}

//get start point if all 3 inputs are available
static void
get_start_point(const size_t t, const size_t u,
                const size_t m, const size_t l,
                const size_t h, const size_t g,
                const double tolerance,
                double &p_m, double &p_h) {

  p_m = static_cast<double>(m)/(m + l + tolerance) + tolerance;
  p_h = static_cast<double>(h)/(h + g + tolerance) + tolerance;
  double p_u = static_cast<double>(u)/(t + u + tolerance) + tolerance;
  double sum = p_m + p_h + p_u;
  p_m = p_m/sum;
  p_h = p_h/sum;
}

static void
expectation(const size_t a, const size_t x,
            const double p, const double q,
            vector<vector<double> > &coeff) {
  assert(p > 0.0 && q > 0.0);
  assert(p + q <= 1.0);
  coeff = vector<vector<double> >(x + 1, vector<double>(a + 1));

  const double log_p = log(p);
  const double log_1mpq = log(1.0 - p - q);
  const double log_q = log(q);
  const double log_1mq = log(1.0 - q);
  const double log_p_q = log(p + q);

  vector<double> a_c_j;
  for (size_t j = 0; j <= a; ++j) {
    a_c_j.push_back(gsl_sf_lnchoose(a, j) + log_q*(a - j) + log_p*j - log_p_q*a);
  }

  for (size_t k = 0; k <= x; ++k) {
    const double x_c_k =
      gsl_sf_lnchoose(x, k) + log_p*k + log_1mpq*(x - k) - log_1mq*x;
    for (size_t j = 0; j <= a; ++j)
      coeff[k][j] = exp(a_c_j[j] + x_c_k);
  }
}


static double
maximization(const size_t x, const size_t y,
             const size_t a, const size_t b,
             const vector<vector<double> > &coeff) {
  double num = y, denom = y + b;
  for (size_t k = 0; k <= x; ++k) {
    vector<double>::const_iterator c(coeff[k].begin());
    for (size_t j = 0; j <= a; ++j) {
      num += (*c)*(a - j);
      denom += (*c)*(a + x - k - j);
      ++c;
    }
  }
  return num/denom;
}


static double
update_p_m(const size_t x, const size_t y,
           const size_t z, const size_t w,
           const size_t a, const size_t b,
           const vector<vector<double> > &coeff) {
  double num = z;
  for (size_t k = 0; k <= x; ++k)
    for (size_t j = 0; j <= a; ++j)
      num += coeff[k][j]*(k + j);
  return num/(a + b + x + y + z + w);
}


static void
expectation_maximization(const bool DEBUG,
                         const size_t x, const size_t y,
                         const size_t z, const size_t w,
                         const size_t a, const size_t b,
                         const double tolerance,
                         double &p, double &q) {
  size_t iter = 0;
  double delta = std::numeric_limits<double>::max();

  if (DEBUG) {
    cerr << "t:" << a << ", u:" << b
         << ", m:" << z << ", l:" << w
         << ", h:" << y << ", g:" << x << endl
         << "p:" << p << ", q:" << q << endl;
  }

  do {
    vector<vector<double> > coeff;
    assert (p > 0 && q > 0);
    expectation(a, x, p, q, coeff);
    const double M = maximization(x, y, a, b, coeff);
    const double p_old = p, q_old = q;
    p = update_p_m(x, y, z, w, a, b, coeff);
    q = M*(1 - p);
    p = max(tolerance, min(p, 1-2*tolerance));
    q = max(tolerance, min(q, 1-tolerance-p));
    delta = max(fabs(p_old - p), fabs(q_old - q));

    iter ++;

  } while (delta > tolerance && iter <=500);

  if (DEBUG) {
    cerr << iter << '\t'
         << "p_m=" << p << '\t'
         << "p_h=" << q << '\t'
         << "log-likelihood=" << log_L(y,x,z,w,b,a,q, p) << endl;
  }
}

//////////////////////////
/*Only 2 input files*/////
//////////////////////////
static double
log_L(const size_t x, const size_t y,
      const size_t z, const size_t w,
      const double p, const double q) {
  assert(p+q <= 1);
  double log_lkhd = gsl_sf_lnchoose(x+y, x) + gsl_sf_lnchoose(z+w, z);
  if (p > 0) log_lkhd += x*log(p);
  if (p < 1) log_lkhd += y*log(1-p);

  if (q > 0) log_lkhd += z*log(q);
  if (q < 1) log_lkhd += w*log(1-q);

  return(log_lkhd);
}

/* Get start point for p and q, if only 2 inputs are available */
static void
get_start_point(const size_t x, const size_t y,
                const size_t z, const size_t w,
                const double tolerance,
                double &p, double &q) {

  p = tolerance + 1.0*x/(x + y);
  q = tolerance + 1.0*z/(z + w);
  if (p + q >= 1.0) {
    p = max(tolerance, (1.0 - tolerance)*p/(p+q));
    q = max(tolerance, 1.0 - p - tolerance);
  }

}

static void
expectation(const size_t y,
            const double p, const double q,
            vector<double> &coeff) {
  assert(p > 0.0 && q > 0.0);
  assert(p + q < 1.0);
  coeff.clear();
  const double log_1mpq = log(1.0 - p - q);
  const double log_q = log(q);
  const double log_1mp = log(1.0 - p);
  for (size_t j = 0; j <= y; ++j)
    coeff.push_back(exp(gsl_sf_lnchoose(y, j) + log_q*j
                        + log_1mpq*(y-j) - log_1mp*y));
}

static double
maximization(const size_t x, const size_t y,
             const vector<double> &coeff) {
  double num = x, denom = x+y;
  for (size_t j = 0; j <= y; ++j) {
    denom -= coeff[j]*j;
  }
  return num/denom;
}

static double
update_q(const size_t x, const size_t y,
         const size_t z, const size_t w,
         const vector<double> &coeff) {
  double num = z;
  double denom = x+y+z+w;
  for (size_t j =0; j <=y; ++j) {
    num += coeff[j]*j;
  }
  return num/denom;
}

static void
expectation_maximization(const bool DEBUG,
                         const size_t x, const size_t y,
                         const size_t z, const size_t w,
                         const double tolerance,
                         double &p, double &q) {
  size_t iter = 0;
  double delta = std::numeric_limits<double>::max();

  do {
    vector<double> coeff;
    expectation(y, p, q, coeff);
    const double M = maximization(x, y, coeff);
    const double p_old = p, q_old = q;
    q = update_q(x, y, z, w, coeff);
    p = M*(1 - q);
    p = max(tolerance, min(p, 1-2*tolerance));
    q = max(tolerance, min(q, 1-tolerance-p));
    delta = max(fabs(p_old - p), fabs(q_old - q));
    iter ++;
  }
  while (delta > tolerance && iter <= 500);

  if (DEBUG) {
    cerr << iter << '\t'
         << "p=" << p << '\t'
         << "q=" << q << '\t'
         <<"log-likelihood=" << log_L(x, y, z, w, p, q) << endl;
  }
}


static void
process_three_types(const double alpha,
                    const double tolerance,
                    const bool FLAG,
                    const string &hydroxy_file,
                    const string &bs_seq_file,
                    const string &oxbs_file,
                    const string &outfile,
                    const string &out_methcount_pseudo_h,
                    const string &out_methcount_pseudo_m,
                    size_t &total_sites,
                    size_t &overshoot_sites,
                    size_t &conflict_sites) {

  ifstream h_in(hydroxy_file);
  ifstream b_in(bs_seq_file);
  ifstream o_in(oxbs_file);

  ofstream out(outfile);
  ofstream out_m, out_h;
  if (!out_methcount_pseudo_m.empty())
    out_m.open(out_methcount_pseudo_m);
  if (!out_methcount_pseudo_h.empty())
    out_h.open(out_methcount_pseudo_h);

  MSite r, s, o;
  while (h_in >> r && b_in >> s && o_in >> o) {

    if (r.n_reads > 500) r.n_reads = 500;
    const size_t h = r.n_meth();
    const size_t g = r.n_unmeth();

    if (s.n_reads > 500) s.n_reads = 500;
    const size_t t = s.n_meth();
    const size_t u = s.n_unmeth();

    if (o.n_reads > 500) o.n_reads = 500;
    const size_t m = o.n_meth();
    const size_t l = o.n_unmeth();

    assert(r.chrom == s.chrom && r.chrom == o.chrom &&
           r.pos == o.pos && r.pos == s.pos);

    total_sites++;
    double p_m, p_h, p_u;
    int CONFLICT = 0, cflt_m, cflt_h, cflt_u;
    double p_m_hat, p_h_hat, p_u_hat;

    size_t x = 0, y = 0, z = 0, w = 0, a = 0, b= 0;

    if ((h + g > 0 && u + t > 0) ||
        (h + g > 0 && m + l > 0) ||
        (m + l > 0 && u + t > 0)) {

      x = g; y = h;
      z = m; w = l;
      a = t; b = u;

      p_h_hat = static_cast<double>(y)/(x + y);
      p_m_hat = static_cast<double>(z)/(z + w);
      p_u_hat = static_cast<double>(b)/(a + b);

      // use frequent method result if no overshoot
      if (p_h_hat + p_m_hat + p_u_hat == 1.0) {
        out << r.chrom << '\t' << r.pos << '\t'
            << r.pos + 1 << '\t' <<p_m_hat << '\t'
            << p_h_hat << '\t' << p_u_hat << "\t0" << endl;

        // write out pseudo methcount files for mC and hmC
        if (!out_methcount_pseudo_m.empty())
          out_m << r.chrom << '\t' << r.pos << "\t+\tCpG\t"
                << p_m_hat << "\t" << a+b+x+y+z+w << endl;
        if (!out_methcount_pseudo_h.empty())
          out_h << r.chrom << '\t' << r.pos << "\t+\tCpG\t"
                << p_h_hat << "\t" << a+b+x+y+z+w << endl;
      }
      else {
        overshoot_sites++;
        get_start_point(t, u, m, l, h, g, tolerance, p_m, p_h);
        expectation_maximization(false, x, y, z, w, a, b, tolerance, p_m, p_h);

        p_u = 1 - p_m - p_h;
        if (p_h <= 2.0*tolerance) p_h = 0.0;
        if (p_m <= 2.0*tolerance) p_m = 0.0;
        if (p_u <= 2.0*tolerance) p_u = 0.0;
        if (p_m >= 1.0-2.0*tolerance) p_m = 1.0;
        if (p_h >= 1.0-2.0*tolerance) p_h = 1.0;
        if (p_u >= 1.0-2.0*tolerance) p_u = 1.0;

        if (p_h_hat + p_m_hat + p_u_hat != 1 && FLAG) {
          cflt_h = binom_null(alpha, static_cast<double>(x+y), p_h_hat, p_h);
          cflt_m = binom_null(alpha, static_cast<double>(z+w), p_m_hat, p_m);
          cflt_u = binom_null(alpha, static_cast<double>(a+b), p_u_hat, p_u);
          CONFLICT = cflt_m + cflt_u + cflt_h;
        }


        out << r.chrom << '\t' << r.pos << '\t'
            << r.pos + 1 << '\t' << p_m << '\t'
            << p_h << "\t" << p_u << "\t" << CONFLICT << endl;

        if (CONFLICT > 1)
          conflict_sites++;

        // write out pseudo methcount files for mC and hmC
        if (!out_methcount_pseudo_m.empty())
          out_m << r.chrom << '\t' << r.pos << "\t+\tCpG\t"
                << p_m << '\t' << a+b+x+y+z+w << endl;
        if (!out_methcount_pseudo_h.empty())
          out_h << r.chrom << '\t' << r.pos << "\t+\tCpG\t"
                << p_h << '\t' << a+b+x+y+z+w << endl;
      }
    }
    else { //observation from only one experiment
      out << r.chrom << '\t' << r.pos << '\t'
          << r.pos +1 << '\t' << "NA\tNA\tNA\tNA" << endl;

      // write out pseudo methcount files for mC and hmC
      if (!out_methcount_pseudo_m.empty()) {
        int coverage = m + l;
        const double level = coverage > 0 ? static_cast<double>(m)/coverage : 0.0;
        out_m << r.chrom << '\t' << r.pos << "\t+\tCpG\t"
              << level << "\t" << coverage << endl;
      }
      if (!out_methcount_pseudo_h.empty()) {
        int coverage = h + g;
        double level = coverage > 0 ? static_cast<double>(h)/coverage : 0.0;
        out_h << r.chrom << '\t' << r.pos << "\t+\tCpG\t"
              << level << "\t" << coverage << endl;
      }
    }
  }
}


static void
process_two_types(const double alpha,
                  const double tolerance,
                  const bool FLAG,
                  const string &hydroxy_file,
                  const string &bs_seq_file,
                  const string &oxbs_file,
                  const string &outfile,
                  const string &out_methcount_pseudo_h,
                  const string &out_methcount_pseudo_m,
                  size_t &total_sites,
                  size_t &overshoot_sites,
                  size_t &conflict_sites) {

  ofstream out(outfile);
  ofstream out_m, out_h;
  if (!out_methcount_pseudo_m.empty())
    out_m.open(out_methcount_pseudo_m);
  if (!out_methcount_pseudo_h.empty())
    out_h.open(out_methcount_pseudo_h);

  std::ifstream f_in, s_in;
  bool f_rev = false, s_rev = false;
  if (oxbs_file.empty()) {
    s_rev = true;
    f_in.open(hydroxy_file);
    s_in.open(bs_seq_file);
  }
  else if (hydroxy_file.empty()) {
    f_rev = true;
    f_in.open(bs_seq_file);
    s_in.open(oxbs_file);
  }
  else {
    f_in.open(oxbs_file);
    s_in.open(hydroxy_file);
  }

  MSite f, s;
  while (f_in >> f && s_in >> s) {

    if (f.chrom != s.chrom || f.pos != s.pos)
        throw runtime_error("error: chrom or pos mismatch in both files");

    if (f.n_reads > 500) f.n_reads = 500;
    size_t x = f.n_meth();
    size_t y = f.n_unmeth();
    if (f_rev) std::swap(x, y);

    if (s.n_reads > 500) s.n_reads = 500;
    size_t z = s.n_meth();
    size_t w = s.n_unmeth();
    if (s_rev) std::swap(z, w);

    total_sites++;
    double p = 0.0, q = 0.0, r = 0.0;
    int CONFLICT = 0, cflt1, cflt2;
    if (x + y > 0 && z + w > 0) {
      if (static_cast<double>(x)/(x+y) +
          static_cast<double>(z)/(z+w) <= 1.0) {
        p = static_cast<double>(x)/(x+y);
        q = static_cast<double>(z)/(z+w);
        r = 1.0 - p - q;
      }
      else {
        overshoot_sites++;
        get_start_point(x, y, z, w, tolerance, p, q);
        expectation_maximization(false, x, y, z, w, tolerance, p, q);
        r = 1.0 - p - q;
        if (p <= 2.0*tolerance) p = 0.0;
        if (q <= 2.0*tolerance) q = 0.0;
        if (r <= 2.0*tolerance) r = 0.0;
        if (p >= 1.0 - 2.0*tolerance) p = 1.0;
        if (q >= 1.0 - 2.0*tolerance) q = 1.0;
        if (r >= 1.0 - 2.0*tolerance) r = 1.0;
        if (FLAG) {
          const double p_hat1 = static_cast<double>(x)/(x+y);
          cflt1 = binom_null(alpha, static_cast<double>(x+y), p_hat1, p);
          const double p_hat2 = static_cast<double>(z)/(z+w);
          cflt2 = binom_null(alpha, static_cast<double>(z+w), p_hat2, q);
          CONFLICT = cflt1 + cflt2;
        }
      }

      out << f.chrom << '\t' << f.pos << '\t' << f.pos + 1 << '\t';
      if (oxbs_file.empty()) out << r << '\t' << p << '\t' << q << '\t';
      else if (hydroxy_file.empty()) out << q << '\t' << r << '\t' << p << '\t';
      else out << p << '\t' << q << '\t' << r << '\t';
      out << CONFLICT << endl;

      if (CONFLICT > 1)
        conflict_sites++;

      // write out pseudo methcount files for mC and hmC
      if (!out_methcount_pseudo_h.empty()) {
        out_h << f.chrom << '\t' << f.pos << "\t+\tCpG\t";
        if (oxbs_file.empty()) out_h << p;
        else if (hydroxy_file.empty()) out_h << r;
        else out_h << q;
        out_h << '\t' << x + y + z+ w << endl;
      }

      if (!out_methcount_pseudo_m.empty()) {
        out_m << f.chrom << '\t' << f.pos << "\t+\tCpG\t";
        if (oxbs_file.empty()) out_m << r;
        else if (hydroxy_file.empty()) out_m << q;
        else out_m << p;
        out_m << '\t' << x + y + z+ w << endl;
      }
    }
    else { // at most one input file has non-zero coverage

      out << f.chrom << '\t' << f.pos << '\t'
          << f.pos + 1 << "\tNA\tNA\tNA\tNA" << endl;

      // write out pseudo methcount files for mC and hmC
      if (!out_methcount_pseudo_h.empty()) {
        int coverage = 0;
        double level = 0.0;
        if (oxbs_file.empty() && x+y > 0) {
          coverage = x + y;
          level = static_cast<double>(x)/(coverage);
        }
        else if (bs_seq_file.empty() && z+w > 0) {
          coverage = z + w;
          level = static_cast<double>(z)/(coverage);
        }
        out_h << s.chrom << '\t' << s.pos << "\t+\tCpG\t"
              << level << '\t' << coverage << endl;
      }

      if (!out_methcount_pseudo_m.empty()) {
        int coverage = 0;
        double level = 0.0;
        if (bs_seq_file.empty() && x + y > 0) {
          coverage = x + y;
          level = static_cast<double>(x)/(coverage);
        }
        else if (hydroxy_file.empty() && z + w > 0) {
          coverage = z + w;
          level = static_cast<double>(z)/(coverage);
        }
        out_m << s.chrom << '\t' << s.pos << "\t+\tCpG\t"
              << level << '\t' << coverage << endl;
      }
    }
  }
}


int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    bool FLAG = true;
    string oxbs_file;
    string hydroxy_file;
    string bs_seq_file;
    string outfile;
    string out_methcount_pseudo_h;
    string out_methcount_pseudo_m;
    double alpha = 0.05;
    static double tolerance = 1e-10;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "program to estimate "
                      "methylation levels", "at least two input files");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("bsseq", 'u', "input BS-seq methcounts file",
                      false , bs_seq_file);
    opt_parse.add_opt("tabseq", 'h', "input TAB-seq methcounts file",
                      false , hydroxy_file);
    opt_parse.add_opt("oxbsseq", 'm', "input oxBS-seq methcounts file",
                      false , oxbs_file);
    opt_parse.add_opt("tolerance", 't', "EM convergence threshold",
                      false , tolerance);
    opt_parse.add_opt("alpha", 'a',
                      "significance level of binomial test for each site",
                      false, alpha);
    opt_parse.add_opt("outh", 'H',
                      "hmC pseudo methcount output file (default: null)",
                      false, out_methcount_pseudo_h);
    opt_parse.add_opt("outm", 'M',
                      "mC pseudo methcount output file (default: null)",
                      false, out_methcount_pseudo_m);
    opt_parse.add_opt("verbose", 'v', "print run statistics", false, VERBOSE);
    opt_parse.set_show_defaults();

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
    if (leftover_args.size() >0) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (alpha <= 0.0 || alpha >= 1.0) {
      cerr << "Please specify a value in (0, 1) for -a option." << endl;
      return EXIT_SUCCESS;
    }
    if ((oxbs_file.empty() && hydroxy_file.empty()) ||
        (oxbs_file.empty() && bs_seq_file.empty()) ||
        (bs_seq_file.empty() && hydroxy_file.empty())) {
      cerr << "Please specify at least 2 bed files as input" << endl;
      return EXIT_SUCCESS;
    }
    tolerance = max(1e-15, min(tolerance, 0.1));
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "Output format:" << endl
           << "chrom    start   end     pm      "
           << "ph      pu      #_of_conflict" << endl;

    size_t total_sites = 0;
    size_t overshoot_sites = 0;
    size_t conflict_sites = 0;

    if (!hydroxy_file.empty() && !bs_seq_file.empty() && !oxbs_file.empty())
      process_three_types(alpha, tolerance, FLAG,
                          hydroxy_file, bs_seq_file, oxbs_file, outfile,
                          out_methcount_pseudo_h, out_methcount_pseudo_m,
                          total_sites, overshoot_sites, conflict_sites);
    else
      process_two_types(alpha, tolerance, FLAG,
                        hydroxy_file, bs_seq_file, oxbs_file, outfile,
                        out_methcount_pseudo_h, out_methcount_pseudo_m,
                        total_sites, overshoot_sites, conflict_sites);

    if (VERBOSE)
      cerr << "total sites: " << total_sites << endl
           << "sites with overshoot: " << overshoot_sites
           << " (" << 1.0*overshoot_sites/total_sites*100 << "%)" << endl
           << "sites conflicting to at least two input: " << conflict_sites
           << " (" << 1.0*conflict_sites/total_sites*100 << "%)" << endl;
  }
  catch (const std::runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
