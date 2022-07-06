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
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <algorithm>

// GSL headers
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>

// smithlab headers
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"

using std::string;
using std::vector;
using std::istringstream;
using std::cerr;
using std::cout;
using std::endl;
using std::istream;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::runtime_error;
using std::sort;
using std::min;
using std::copy;

/***************** REGRESSION *****************/
struct Design {
  vector<string> factor_names;
  vector<string> sample_names;
  vector<vector<double> > matrix;
};

struct SiteProportions {
  string chrom;
  size_t position;
  string strand;
  string context;
  vector<size_t> total;
  vector<size_t> meth;
};

struct Regression {
  Design design;
  SiteProportions props;
  double max_loglik;
};

istream&
operator>> (istream &is, Design &design) {
  string header_encoding;
  getline(is, header_encoding);

  istringstream header_is(header_encoding);
  string header_name;
  while (header_is >> header_name)
    design.factor_names.push_back(header_name);

  string row;
  while (getline(is, row)) {

    if (row.empty())
      continue;

    istringstream row_is(row);
    string token;
    row_is >> token;
    design.sample_names.push_back(token);

    vector<double> matrix_row;
    while (row_is >> token) {
      if (token.length() == 1 && (token == "0" || token == "1"))
        matrix_row.push_back(token == "1");
      else
        throw runtime_error("only binary factor levels are allowed:\n"
                            + row);
    }

    if (matrix_row.size() != design.factor_names.size())
      throw runtime_error("each row must have as many columns as "
                          "factors:\n" + row);

    design.matrix.push_back(vector<double>());
    swap(design.matrix.back(), matrix_row);
  }
  return is;
}

ostream&
operator<< (ostream &os, const Design &design) {
  for(size_t factor = 0; factor < design.factor_names.size(); ++factor) {
    os << design.factor_names[factor];
    if (factor + 1 != design.factor_names.size())
      os << "\t";
  }
  os << endl;

  for(size_t sample = 0; sample < design.sample_names.size(); ++sample) {
    os << design.sample_names[sample] << "\t";
    for(size_t factor = 0; factor < design.factor_names.size(); ++factor) {
      os << design.matrix[sample][factor];
      if (factor + 1 != design.factor_names.size())
        os << "\t";
    }
    os << "\n";
  }
  return os;
}

void
remove_factor(Design &design, size_t factor) {
  design.factor_names.erase(design.factor_names.begin() + factor);
  for (size_t sample = 0; sample < design.sample_names.size(); ++sample)
    design.matrix[sample].erase(design.matrix[sample].begin() + factor);
}

// Parses a natural number from its string representation. Throws exception if
// the string does not encode one.
static size_t
parse_natural_number(string encoding) {
  istringstream iss(encoding);
  size_t number;
  iss >> number;
  if (!iss)
    throw runtime_error("The token \"" +encoding + "\" "
                        "does not encode a natural number");
  return number;
}

istream&
operator>>(istream &table_encoding, SiteProportions &props) {
  props.chrom.clear();
  props.position = 0;
  props.strand.clear();
  props.context.clear();
  props.meth.clear();
  props.total.clear();

  string row_encoding;
  getline(table_encoding, row_encoding);

  // Skip lines contining only the newline character (e.g. the last line of the
  // proportion table).
  if(row_encoding.empty())
    return table_encoding;

  // Get the row name (which must be specified like this: "chr:position") and
  // parse it.
  istringstream row_stream(row_encoding);
  string row_name_encoding;
  row_stream >> row_name_encoding;

  // Every row must start an identifier consisiting of genomic loci of the
  // corresponding site. Here we check this identifier has the correct number
  // of colons.
  const size_t num_colon =
    count(row_name_encoding.begin(), row_name_encoding.end(), ':');

  if (num_colon != 3)
    throw runtime_error("Each row in the count table must start with "
                        "a line chromosome:position:strand:context."
                        "Got \"" + row_name_encoding + "\" instead." );

  // First parse the row identifier.
  istringstream name_stream(row_name_encoding);
  getline(name_stream, props.chrom, ':');

  if (props.chrom.empty())
    throw runtime_error("Error parsing " + row_name_encoding +
                        ": chromosome name is missing.");

  string position_encoding;

  getline(name_stream, position_encoding, ':');
  props.position = parse_natural_number(position_encoding);
  getline(name_stream, props.strand, ':');
  getline(name_stream, props.context, ':');

  // After parsing the row identifier, parse count proportions.
  size_t total_count, meth_count;

  while (row_stream >> total_count >> meth_count) {
    props.total.push_back(total_count);
    props.meth.push_back(meth_count);
  }

  if (!row_stream.eof())
    throw runtime_error("Some row entries are not natural numbers: " +
                        row_stream.str());

  if (props.total.size() != props.meth.size())
    throw runtime_error("This row does not encode proportions"
                        "correctly:\n" + row_encoding);
  return table_encoding;
}

static double
pi(Regression *reg, size_t sample, const gsl_vector *parameters) {
  double dot_prod = 0;

  for(size_t factor = 0; factor < reg->design.factor_names.size(); ++factor)
    dot_prod +=
      reg->design.matrix[sample][factor]*gsl_vector_get(parameters, factor);

  double p = exp(dot_prod)/(1 + exp(dot_prod));

  return p;
}

static double
neg_loglik(const gsl_vector *parameters, void *object) {
  Regression *reg = (Regression *)(object);
  const size_t num_parameters = reg->design.factor_names.size() + 1;

  double log_lik = 0;

  //dispersion parameter phi is the last element of parameter vector
  const double dispersion_param = gsl_vector_get(parameters,
                                                 num_parameters - 1);
  const double phi = exp(dispersion_param)/(1 + exp(dispersion_param));

  for(size_t s = 0; s < reg->design.sample_names.size(); ++s) {
    const double n_s = reg->props.total[s];
    const double y_s = reg->props.meth[s];
    const double p_s = pi(reg, s, parameters);

    for(int k = 0; k < y_s; ++k) {
      log_lik += log((1 - phi)*p_s + phi*k);
    }

    for(int k = 0; k < n_s - y_s; ++k) {
      log_lik += log((1 - phi)*(1 - p_s) + phi*k);
    }

    for(int k = 0; k < n_s; ++k) {
      log_lik -= log(1 + phi*(k - 1));
    }
  }

  return (-1)*log_lik;
}

static void
neg_gradient(const gsl_vector *parameters, void *object,
             gsl_vector *output) {

  Regression *reg = (Regression *)(object);
  const size_t num_parameters = reg->design.factor_names.size() + 1;

  const double dispersion_param = gsl_vector_get(parameters,
                                                 num_parameters - 1);

  const double phi = exp(dispersion_param)/(1 + exp(dispersion_param));

  for(size_t f = 0; f < num_parameters; ++f) {

    double deriv = 0;

    for(size_t s = 0; s < reg->design.sample_names.size(); ++s) {
      int n_s = reg->props.total[s];
      int y_s = reg->props.meth[s];
      double p_s = pi(reg, s, parameters);

      double term = 0;

      //a parameter linked to p
      if(f < reg->design.factor_names.size()) {
        double factor = (1 - phi)*p_s*(1 - p_s)*reg->design.matrix[s][f];
        if (factor == 0) continue;

        for(int k = 0; k < y_s; ++k)
          term += 1/((1 - phi)*p_s + phi*k);

        for(int k = 0; k < n_s - y_s; ++k)
          term -= 1/((1 - phi)*(1 - p_s) + phi*k);

        deriv += term*factor;
      } else { // the parameter linked to phi
        for(int k = 0; k < y_s; ++k)
          term += (k - p_s)/((1 - phi)*p_s + phi*k);

        for(int k = 0; k < n_s - y_s; ++k)
          term += (k - (1 - p_s))/((1 - phi)*(1 - p_s) + phi*k);

        for(int k = 0; k < n_s; ++k) {
          term -= (k - 1)/(1 + phi*(k - 1));
        }

        deriv += term * phi * (1 - phi);
      }
    }

    gsl_vector_set(output, f, deriv);
  }

  gsl_vector_scale(output, -1.0);
}

static void
neg_loglik_and_grad(const gsl_vector *parameters,
                    void *object,
                    double *loglik_val,
                    gsl_vector *d_loglik_val) {

  *loglik_val = neg_loglik(parameters, object);
  neg_gradient(parameters, object, d_loglik_val);
}

bool
fit(Regression &r, vector<double> initial_parameters) {
  const size_t num_parameters = r.design.factor_names.size() + 1;

  if (initial_parameters.empty()) {
    for(size_t ind = 0; ind < num_parameters - 1; ++ind)
      initial_parameters.push_back(0.0);
    initial_parameters.push_back(-2.5);
  }

  if (initial_parameters.size() != num_parameters)
    throw runtime_error("Wrong number of initial parameters.");

  int status = 0;

  size_t iter = 0;

  gsl_multimin_function_fdf loglik_bundle;

  loglik_bundle.f = &neg_loglik;
  loglik_bundle.df = &neg_gradient;
  loglik_bundle.fdf = &neg_loglik_and_grad;
  loglik_bundle.n = num_parameters;
  loglik_bundle.params = (void *)&r;

  gsl_vector *parameters = gsl_vector_alloc(num_parameters);

  for (size_t parameter = 0; parameter < initial_parameters.size();
       ++parameter) {
    gsl_vector_set(parameters, parameter, initial_parameters[parameter]);
  }

  const gsl_multimin_fdfminimizer_type *T;

  //can also try gsl_multimin_fdfminimizer_conjugate_pr;
  T = gsl_multimin_fdfminimizer_conjugate_fr;

  gsl_multimin_fdfminimizer *s;
  s = gsl_multimin_fdfminimizer_alloc (T, num_parameters);

  gsl_multimin_fdfminimizer_set (s, &loglik_bundle, parameters, 0.001, 1e-4);

  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);

    if (status)
      break;

    status = gsl_multimin_test_gradient (s->gradient, 1e-4);
  }
  while (status == GSL_CONTINUE && iter < 700);
  //It it reasonable to reduce the number of iterations to 500?

  r.max_loglik = (-1)*neg_loglik(s->x, &r);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(parameters);

  return status == GSL_SUCCESS;
}
bool
fit(Regression &r) {
  return fit(r, vector<double>());
}

/***************** COMBINE P-VALUES *****************/

struct PvalLocus {
  size_t pos;
  double raw_pval;
  double combined_pval;
  double corrected_pval;
};

class BinForDistance {
public:
  BinForDistance(string spec_string);

  size_t which_bin(size_t value) const;

  size_t num_bins() const {return num_bins_;}
  size_t invalid_bin() const {return invalid_bin_;}
  size_t max_dist() const {return max_dist_;}

private:
  size_t min_dist_;
  size_t max_dist_;
  size_t bin_size_;
  size_t num_bins_;
  size_t invalid_bin_;
};

class ProximalLoci {
public:
  ProximalLoci(vector<PvalLocus> &loci, size_t max_distance)
    : loci_(loci), max_distance_(max_distance), next_pos_(loci.begin()) {};
  bool get(vector<PvalLocus> &neighbors);
  PvalLocus cur_region() {return *(next_pos_ - 1);}

private:
  const vector<PvalLocus> &loci_;
  size_t max_distance_;
  vector<PvalLocus>::const_iterator next_pos_;
};

class DistanceCorrelation {
public:
  DistanceCorrelation(BinForDistance bin_for_dist)
    : bin_for_dist_(bin_for_dist) {};
  vector<double> correlation_table(const vector<PvalLocus> &loci);

private:
  double correlation(const vector<double> &x,
                      const vector<double> &y);
  void bin(const vector<PvalLocus> &loci);
  vector< vector<double> > x_pvals_for_bin_;
  vector< vector<double> > y_pvals_for_bin_;
  const BinForDistance bin_for_dist_;
};

static double
to_zscore(double pval) {
  static const double local_epsilon = 1e-6;

  if (pval > 1.0 - local_epsilon)
    pval = 1.0 - local_epsilon;
  else if (pval < local_epsilon)
    pval = local_epsilon;

  return gsl_cdf_ugaussian_Pinv(1.0 - pval);
}

static double
stouffer_liptak(const vector<vector<double> > &corr_mat, vector<double> &pvals) {

  double correction = 0;
  for (size_t row_ind = 0; row_ind < corr_mat.size(); ++row_ind)
    for (size_t col_ind = row_ind + 1; col_ind < corr_mat.size(); ++col_ind)
      correction += corr_mat[row_ind][col_ind];

  vector<double> zscores;
  transform(pvals.begin(), pvals.end(), std::back_inserter(zscores), to_zscore);

  const double sum = std::accumulate(zscores.begin(), zscores.end(), 0.0);
  const double test_stat =
    sum/sqrt(static_cast<double>(pvals.size()) + 2.0*correction);

  return 1.0 - gsl_cdf_gaussian_P(test_stat, 1.0);
}

void
update_pval_loci(std::istream &input_encoding,
                 const vector<PvalLocus> &pval_loci,
                 std::ostream &output_encoding) {

  string record, chrom, name, sign;
  size_t position, coverage_factor, meth_factor, coverage_rest, meth_rest;
  double pval;

  vector<PvalLocus>::const_iterator cur_locus_iter = pval_loci.begin();

  while (getline(input_encoding, record)) {
    // ADS: this seems not to be done well; the code should exit in a
    // "normal" state if bad parse
    try {
      std::istringstream iss(record);
      iss.exceptions(std::ios::failbit);
      iss >> chrom >> position >> sign >> name >> pval
          >> coverage_factor >> meth_factor >> coverage_rest >> meth_rest;
    }
    catch (std::exception const & err) {
      cerr << err.what() << endl << "could not parse line:\n"
           << record << endl;
      std::terminate();
    }

    output_encoding << chrom << "\t" << position << "\t" << sign << "\t"
                    << name << "\t" << pval << "\t";

    if (0.0 <= pval && pval <= 1.0) {
      output_encoding << cur_locus_iter->combined_pval << "\t"
                      << cur_locus_iter->corrected_pval << "\t";
      cur_locus_iter++;
    }
    else output_encoding << -1 << "\t" << -1 << pval << "\t"; // MAGIC??

    output_encoding << coverage_factor << "\t" << meth_factor << "\t"
                    << coverage_rest << "\t" << meth_rest << endl;
  }
}

BinForDistance::BinForDistance(std::string spec_string) {
  std::replace(spec_string.begin(), spec_string.end(), ':', ' ');

  std::istringstream iss(spec_string);
  iss >> min_dist_ >> max_dist_ >> bin_size_;

  num_bins_ = (max_dist_ - min_dist_) / bin_size_;
  invalid_bin_ = num_bins_ + 1;
}

size_t
BinForDistance::which_bin(size_t value) const {
  if (value < min_dist_)
    return invalid_bin_;

  const size_t bin = (value - min_dist_)/bin_size_;

  //Bin numbering is 0 based.
  if (bin >= num_bins_)
    return invalid_bin_;

  return bin;
}

bool
ProximalLoci::get(vector<PvalLocus> &neighbors) {

  if (next_pos_ == loci_.end())
    return false;

  vector<PvalLocus>::const_iterator cur_pos = next_pos_;
  neighbors.clear();
  neighbors.push_back(*cur_pos);

  if ( cur_pos != loci_.begin() ) {
    vector<PvalLocus>::const_iterator up_pos = cur_pos;
    bool too_far = false;

    do {
      --up_pos;
      size_t up_dist = cur_pos->pos - (up_pos->pos + 1);

      if(up_dist <= max_distance_) {
          neighbors.push_back(*up_pos);
      } else
        too_far = true;

    } while (!too_far && up_pos != loci_.begin());
  }

  std::reverse(neighbors.begin(), neighbors.end());

  if (cur_pos != loci_.end() - 1) {
    bool too_far = false;
    vector<PvalLocus>::const_iterator down_pos = cur_pos;

    do {
      ++down_pos;
      size_t down_dist = down_pos->pos - (cur_pos->pos + 1);

      if (down_dist <= max_distance_) {
          neighbors.push_back(*down_pos);
      }
      else too_far = true;

    } while (!too_far && down_pos != loci_.end() - 1);
  }

  ++next_pos_;
  return true;
}

void
DistanceCorrelation::bin(const vector<PvalLocus> &loci) {
  x_pvals_for_bin_.clear();
  y_pvals_for_bin_.clear();
  vector<PvalLocus>::const_iterator it = loci.begin();

  while (it != loci.end()) {
    vector<PvalLocus>::const_iterator forward_it = it + 1;
    bool too_far = false;

    while (forward_it != loci.end() && !too_far) {
      const size_t dist = forward_it->pos - (it->pos + 1);
      const size_t bin = bin_for_dist_.which_bin(dist);

      //check if the appropriate bin exists
      if (bin != bin_for_dist_.invalid_bin()) {
        x_pvals_for_bin_[bin].push_back(to_zscore(it->raw_pval));
        y_pvals_for_bin_[bin].push_back(to_zscore(forward_it->raw_pval));
      }

      if (dist > bin_for_dist_.max_dist())
        too_far = true;

      ++forward_it;
    }

    ++it;
  }
}

double
DistanceCorrelation::correlation(const vector<double> &x,
                                    const vector<double> &y) {
  //Correlation is 0 when all bins are empty.
  if (x.size() <= 1)
    return 0;

  gsl_vector_const_view gsl_x = gsl_vector_const_view_array(&x[0], x.size());
  gsl_vector_const_view gsl_y = gsl_vector_const_view_array(&y[0], y.size());
  const size_t stride = 1;
  double corr = gsl_stats_correlation( gsl_x.vector.data, stride,
                                         gsl_y.vector.data, stride,
                                         x.size());
  return corr;
}


vector<double>
DistanceCorrelation::correlation_table(const vector<PvalLocus> &loci) {
  const size_t num_bins = bin_for_dist_.num_bins();
  x_pvals_for_bin_.resize(num_bins);
  y_pvals_for_bin_.resize(num_bins);
  bin(loci);
  vector<double> correlation_table;

  for (size_t bin = 0; bin < num_bins; ++bin) {
    const double corr = correlation(x_pvals_for_bin_[bin],
                                      y_pvals_for_bin_[bin]);
    correlation_table.push_back(corr);
  }

  return correlation_table;
}

void
distance_corr_matrix(BinForDistance bin_for_dist,
                     const std::vector<double> &acor_for_bin,
                     const std::vector<PvalLocus> &neighbors,
                     std::vector< std::vector<double> > &corr_matrix) {
  corr_matrix.clear();
  const size_t num_neighbors = neighbors.size();

  corr_matrix.resize(num_neighbors);

  for (std::vector<std::vector<double> >::iterator row = corr_matrix.begin();
        row != corr_matrix.end(); ++row)
    row->resize(num_neighbors);

  for (size_t row = 0; row < num_neighbors; ++row) {
    corr_matrix[row][row] = 1.0;
    const size_t row_locus = neighbors[row].pos + 1;

    for (size_t col = row + 1; col < num_neighbors; ++col) {
      const size_t col_locus = neighbors[col].pos;
      const size_t dist = col_locus - row_locus;

      const size_t bin = bin_for_dist.which_bin(dist);

      if (bin == bin_for_dist.invalid_bin()) {
        corr_matrix[row][col] = 0;
      } else
        corr_matrix[row][col] = corr_matrix[col][row] = acor_for_bin[bin];
    }
  }
}

void
combine_pvals(vector<PvalLocus> &loci, const BinForDistance &bin_for_distance) {
  DistanceCorrelation distance_correlation(bin_for_distance);
  vector<double> correlation_for_bin =
                                  distance_correlation.correlation_table(loci);
  ProximalLoci proximal_loci(loci, bin_for_distance.max_dist());
  vector<double> combined_pvalues;
  vector<PvalLocus> neighbors;
  size_t i = 0;

  while (proximal_loci.get(neighbors)) {
   vector< vector<double> > correlation_matrix;
   vector<double> p_vals;

   for (vector<PvalLocus>::const_iterator it = neighbors.begin();
        it != neighbors.end(); ++it) {
      double pval = it->raw_pval;
      p_vals.push_back(pval);
    }

    distance_corr_matrix(bin_for_distance, correlation_for_bin,
                         neighbors, correlation_matrix);
    double combined_pval = stouffer_liptak(correlation_matrix, p_vals);
    loci[i].combined_pval = combined_pval;

    i++;
  }
}
/***************** RADMETH ALGORITHM *****************/
static bool
lt_locus_pval(const PvalLocus &r1, const PvalLocus &r2) {
  return r1.combined_pval < r2.combined_pval;
}

static bool
ls_locus_position(const PvalLocus &r1, const PvalLocus &r2) {
  return r1.pos < r2.pos;
}

void
fdr(vector<PvalLocus> &loci) {

  sort(begin(loci), end(loci), lt_locus_pval);

  for (size_t idx = 0; idx < loci.size(); ++idx) {
    const double current_score = loci[idx].combined_pval;
    //Assign a new one.
    const double corrected_pval = loci.size()*current_score/(idx + 1);
    loci[idx].corrected_pval = corrected_pval;
  }

  for (vector<PvalLocus>::reverse_iterator it = loci.rbegin() + 1;
       it != loci.rend(); ++it) {

    const PvalLocus &prev_locus = *(it - 1);
    PvalLocus &cur_locus = *(it);

    cur_locus.corrected_pval =
      min(prev_locus.corrected_pval, cur_locus.corrected_pval);
  }

  for (auto it(begin(loci)); it != end(loci); ++it)
    it->corrected_pval = min(it->corrected_pval, 1.0);

  // restore original order
  sort(begin(loci), end(loci), ls_locus_position);
}

// Splits a string using white-space characters as delimeters.
static vector<string>
split(string input) {
  istringstream iss(input);
  string token;
  vector<string> tokens;

  while (iss >> token)
    tokens.push_back(token);

  return tokens;
}

// Given the maximum likelihood estimates of the full and reduced models, the
// function outputs the p-value of the log-likelihood ratio. *Note* that it is
// assumed that the reduced model has one fewer factor than the reduced model.
double
loglikratio_test(double null_loglik, double full_loglik) {

  // The log-likelihood ratio statistic.
  const double log_lik_stat = -2*(null_loglik - full_loglik);

  // It is assumed that null model has one fewer factor than the full model.
  // Hence the number of degrees of freedom is 1.
  const size_t degrees_of_freedom = 1;

  // Log-likelihood ratio statistic has a chi-sqare distribution.
  double chisq_p = gsl_cdf_chisq_P(log_lik_stat, degrees_of_freedom);
  const double pval = 1.0 - chisq_p;

  return pval;
}

bool
has_low_coverage(const Regression &reg, const size_t test_factor) {

  bool is_covered_in_test_factor_samples = false;
  bool is_covered_in_other_samples = false;

  for (size_t sample = 0; sample < reg.design.sample_names.size(); ++sample) {
    if (reg.design.matrix[sample][test_factor] == 1) {
      if (reg.props.total[sample] != 0)
        is_covered_in_test_factor_samples = true;
    }
    else {
      if (reg.props.total[sample] != 0)
        is_covered_in_other_samples = true;
    }
  }

  return !is_covered_in_test_factor_samples || !is_covered_in_other_samples;
}

bool
has_extreme_counts(const Regression &reg) {

  bool is_maximally_methylated = true;
  bool is_unmethylated = true;

  for (size_t sample = 0; sample < reg.design.sample_names.size(); ++sample) {
    if (reg.props.total[sample] != reg.props.meth[sample])
      is_maximally_methylated = false;

    if (reg.props.meth[sample] != 0)
      is_unmethylated = false;
  }

  return is_maximally_methylated || is_unmethylated;
}

/***********************************************************************
 * Run beta-binoimial regression using the specified table with
 * proportions and design matrix
 */
static int
run_regression(int argc, const char **argv) {

  const string command_name = argv[0];

  string outfile;
  string test_factor_name;
  bool VERBOSE = false;

  OptionParser opt_parse(command_name,
                         "calculate differential methylation scores",
                         "<design-matrix> <data-matrix>");
  opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                    false, outfile);
  opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
  opt_parse.add_opt("factor", 'f', "a factor to test", true, test_factor_name);

  vector<string> leftover_args;
  opt_parse.parse(argc, argv, leftover_args);
  if (argc == 2 || opt_parse.help_requested()) {
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
  if (leftover_args.size() != 2) {
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }
  const string design_filename(leftover_args.front());
  const string table_filename(leftover_args.back());

  ifstream design_file(design_filename);
  if (!design_file)
    throw runtime_error("could not open file: " + design_filename);

  ifstream table_file(table_filename);
  if (!table_file)
    throw runtime_error("could not open file: " + table_filename);

  ofstream of;
  if (!outfile.empty()) of.open(outfile);
  ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

  // initialize full design matrix from file
  Regression full_regression;
  design_file >> full_regression.design;

  // Check that provided test factor name exists and find its index.
  // Identify test factors with their indexes to simplify naming
  vector<string>::const_iterator test_factor_it =
    find(begin(full_regression.design.factor_names),
              end(full_regression.design.factor_names), test_factor_name);

  if (test_factor_it == end(full_regression.design.factor_names))
    throw runtime_error("Error: " + test_factor_name +
                        " is not a part of the design specification.");

  const size_t test_factor = test_factor_it -
    begin(full_regression.design.factor_names);

  Regression null_regression;
  null_regression.design = full_regression.design;
  remove_factor(null_regression.design, test_factor);

  // Make sure that the first line of the proportion table file contains
  // names of the samples. Throw an exception if the names or their order
  // in the proportion table does not match those in the full design matrix.
  string sample_names_encoding;
  getline(table_file, sample_names_encoding);

  if (full_regression.design.sample_names != split(sample_names_encoding))
    throw runtime_error(sample_names_encoding + " does not match factor "
                        "names or their order in the design matrix. "
                        "Please verify that the design matrix and the "
                        "proportion table are correctly formatted.");

  // Performing the log-likelihood ratio test on proportions from each row
  // of the proportion table.
  while (table_file >> full_regression.props) {

    if (full_regression.design.sample_names.size() !=
        full_regression.props.total.size())
      throw runtime_error("found row with wrong number of columns");

    size_t coverage_factor = 0, coverage_rest = 0,
      meth_factor = 0, meth_rest = 0;

    for(size_t s = 0; s < full_regression.design.sample_names.size(); ++s) {
      if (full_regression.design.matrix[s][test_factor] != 0) {
        coverage_factor += full_regression.props.total[s];
        meth_factor += full_regression.props.meth[s];
      }
      else {
        coverage_rest += full_regression.props.total[s];
        meth_rest += full_regression.props.meth[s];
      }
    }

    out << full_regression.props.chrom << "\t"
        << full_regression.props.position << "\t"
        << full_regression.props.strand << "\t"
        << full_regression.props.context << "\t";

    // Do not perform the test if there's no coverage in either all case or
    // all control samples. Also do not test if the site is completely
    // methylated or completely unmethylated across all samples.
    if (has_low_coverage(full_regression, test_factor)) {
      out << -1;
    }
    else if (has_extreme_counts(full_regression)) {
      out << -1;
    }
    else {
      fit(full_regression);
      null_regression.props = full_regression.props;
      fit(null_regression);
      const double pval = loglikratio_test(null_regression.max_loglik,
                                           full_regression.max_loglik);

      // If error occured in fitting (p-val = nan or -nan).
      out << ((pval != pval) ? -1 : pval);
    }
    out << "\t" << coverage_factor << "\t" << meth_factor
        << "\t" << coverage_rest << "\t" << meth_rest << endl;
  }

  return EXIT_SUCCESS;
}

static int
run_adjust(int argc, const char **argv) {

  // first argument is name of command
  const string command_name = argv[0];

  string outfile;
  string bin_spec = "1:200:1";
  bool VERBOSE = false;

  /**************** GET COMMAND LINE ARGUMENTS *************************/
  OptionParser opt_parse(command_name,
                         "compute adjusted p-values using autocorrelation",
                         "<regression-output>");
  opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                    false, outfile);
  opt_parse.add_opt("bins", 'b', "corrlation bin specs", false , bin_spec);
  opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
  vector<string> leftover_args;
  opt_parse.parse(argc, argv, leftover_args);
  if (argc == 2 || opt_parse.help_requested()) {
    cerr << opt_parse.help_message() << endl;
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
  const string bed_filename = leftover_args.front();
  /*********************************************************************/

  BinForDistance bin_for_dist(bin_spec);

  ifstream bed_file(bed_filename);
  if (!bed_file)
    throw "could not open file: " + bed_filename;

  if (VERBOSE)
    cerr << "[reading input]" << endl;

  // Read in all p-value loci. The loci that are not correspond to valid
  // p-values (i.e. values in [0, 1]) are skipped.
  vector<PvalLocus> pvals;
  string input_line, prev_chrom;

  size_t chrom_offset = 0;

  while (getline(bed_file, input_line)) {

    istringstream iss(input_line);
    string chrom, sign, name;
    size_t position;
    double pval;
    if (!(iss >> chrom >> position >> sign >> name >> pval))
      throw runtime_error("failed to parse line: " + input_line);

    // Skip loci that do not correspond to valid p-values.
    if (0 <= pval && pval <= 1) {
      // locus is on new chrom.
      if (!prev_chrom.empty() && prev_chrom != chrom)
        chrom_offset += pvals.back().pos;

      PvalLocus plocus;
      plocus.raw_pval = pval;
      plocus.pos = chrom_offset + bin_for_dist.max_dist() + 1 + position;

      pvals.push_back(plocus);
      prev_chrom = chrom;
    }
  }

  if (VERBOSE)
    cerr << "[combining p-values]" << endl;
  combine_pvals(pvals, bin_for_dist);

  if (VERBOSE)
    cerr << "[running multiple test adjustment]" << endl;
  fdr(pvals);

  ofstream of;
  if (!outfile.empty()) of.open(outfile);
  ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

  ifstream original_bed_file(bed_filename);

  update_pval_loci(original_bed_file, pvals, out);

  //TODO: Check that the regions do not overlap & sorted
  return EXIT_SUCCESS;
}

// Attemps to find the next significant CpG site. Returns true if one was found
// and flase otherwise.
static bool
read_next_significant_cpg(istream &cpg_stream, GenomicRegion &cpg,
                          double cutoff, bool &skipped_any, bool &n_sig_sites,
                          size_t &test_cov, size_t &test_meth,
                          size_t &rest_cov, size_t &rest_meth) {
  GenomicRegion region;
  skipped_any = false;
  n_sig_sites = false;
  string cpg_encoding;

  while (getline(cpg_stream, cpg_encoding)) {
    string record, chrom, name, sign;
    size_t position;
    double raw_pval, adjusted_pval, corrected_pval;

    istringstream iss(cpg_encoding);
    iss.exceptions(std::ios::failbit);
    iss >> chrom >> position >> sign >> name >> raw_pval
        >> adjusted_pval >> corrected_pval
        >> test_cov >> test_meth >> rest_cov >> rest_meth;

    if (0 <= corrected_pval && corrected_pval < cutoff) {
      cpg.set_chrom(chrom);
      cpg.set_start(position);
      cpg.set_end(position + 1);
      n_sig_sites = (0 <= raw_pval && raw_pval < cutoff);
      return true;
    }
    skipped_any = true;
  }

  return false;
}

static void
merge(istream &cpg_stream, ostream &dmr_stream, double cutoff) {

  GenomicRegion dmr;
  dmr.set_name("dmr");

  size_t dmr_test_cov = 0;
  size_t dmr_test_meth = 0;
  size_t dmr_rest_cov = 0;
  size_t dmr_rest_meth = 0;

  size_t test_cov = 0;
  size_t test_meth = 0;
  size_t rest_cov = 0;
  size_t rest_meth = 0;

  // Find the first significant CpG, or terminate the function if none exist.
  bool skipped_last_cpg, n_sig_sites;
  if (!read_next_significant_cpg(cpg_stream, dmr, cutoff, skipped_last_cpg,
                                 n_sig_sites, test_cov, test_meth,
                                 rest_cov, rest_meth))
    return;

  dmr.set_score(n_sig_sites);
  dmr_test_cov += test_cov;
  dmr_test_meth += test_meth;
  dmr_rest_cov += rest_cov;
  dmr_rest_meth += rest_meth;

  GenomicRegion cpg;
  cpg.set_name("dmr");

  while (read_next_significant_cpg(cpg_stream, cpg, cutoff, skipped_last_cpg,
                                   n_sig_sites, test_cov, test_meth,
                                   rest_cov, rest_meth)) {

    if (skipped_last_cpg || cpg.get_chrom() != dmr.get_chrom()) {
      if (dmr.get_score() != 0)
        dmr_stream << dmr.get_chrom() << '\t'
                   << dmr.get_start() << '\t'
                   << dmr.get_end()   << '\t'
                   << dmr.get_name()  << '\t'
                   << dmr.get_score() << '\t'
                   << double(dmr_test_meth)/dmr_test_cov -
          double(dmr_rest_meth)/dmr_rest_cov << endl;
      dmr = cpg;
      dmr.set_score(n_sig_sites);
      dmr_test_cov = test_cov;
      dmr_test_meth = test_meth;
      dmr_rest_cov = rest_cov;
      dmr_rest_meth = rest_meth;
    }
    else {
      dmr.set_end(cpg.get_end());
      dmr.set_score(dmr.get_score() + n_sig_sites);
      dmr_test_cov += test_cov;
      dmr_test_meth += test_meth;
      dmr_rest_cov += rest_cov;
      dmr_rest_meth += rest_meth;
    }
  }
  if (dmr.get_score() != 0) {
    dmr_stream << dmr.get_chrom() << '\t'
               << dmr.get_start() << '\t'
               << dmr.get_end()   << '\t'
               << dmr.get_name()  << '\t'
               << dmr.get_score() << '\t'
               << double(dmr_test_meth)/dmr_test_cov -
      double(dmr_rest_meth)/dmr_rest_cov << endl;
  }
}

static int
run_merge(int argc, const char **argv) {

  // first argument is name of command
  const string command_name = argv[0];

  /* FILES */
  string outfile;
  string bin_spec = "1:200:25";
  double cutoff = 0.01;

  /**************** GET COMMAND LINE ARGUMENTS *************************/
  OptionParser opt_parse(command_name,
                         "merge significantly differentially"
                         " methylated CpGs into DMRs",
                         "<bed-file-in-radmeth-format>");
  opt_parse.add_opt("output", 'o',
                    "output file (default: stdout)", false, outfile);
  opt_parse.add_opt("cutoff", 'p', "P-value cutoff (default: 0.01)",
                    false , cutoff);
  vector<string> leftover_args;
  opt_parse.parse(argc, argv, leftover_args);
  if (argc == 1 || opt_parse.help_requested()) {
    cerr << opt_parse.help_message() << endl;
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
  const string bed_filename = leftover_args.front();
  /************************************************************************/

  ofstream of;
  if (!outfile.empty()) of.open(outfile);
  ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

  ifstream in(bed_filename);
  if (!in)
    throw runtime_error("could not open file: " + bed_filename);

  merge(in, out, cutoff);

  return EXIT_SUCCESS;
}


int
main(int argc, const char **argv) {

  try {
    const string prog_name = strip_path(argv[0]);

    const string main_help_message =
      "Uage: " + prog_name + " [COMMAND] [PARAMETERS]\n\n"
      "Available commands: \n"
      "  regression  Calculates multi-factor differential methylation scores.\n"
      "  adjust      Adjusts the p-value of each site based on the p-value of "
      "its neighbors.\n"
      "  merge       Combines significantly differentially methylated CpGs into"
      " DMRs.\n";

    if (argc == 1) {
      cerr << "Analysis of differential methylation in"
           << " multi-factor bisulfite sequencing experiments" << endl
           << main_help_message;
      return EXIT_SUCCESS;
    }

    // first argument is name of command
    const string command_name = argv[1];

    if (command_name == "regression")
      return run_regression(argc - 1, argv + 1);
    else if (command_name == "adjust")
      return run_adjust(argc - 1, argv + 1);
    else if (command_name == "merge")
      return run_merge(argc - 1, argv + 1);
    else {
      cerr << "ERROR: invalid command name: \""
           << command_name << "\"" << endl
           << main_help_message;
      return EXIT_SUCCESS;
    }
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
