/*
  Copyright (C) 2008 Cold Spring Harbor Laboratory
  Authors: Andrew D. Smith

  This file is part of rmap.

  rmap is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  rmap is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with rmap; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include "HMM.hpp"

#include <iomanip>
#include <numeric>
#include <limits>
#include <cmath>

#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>


using std::vector;
using std::pair;
using std::setw;
using std::max;
using std::cerr;
using std::endl;
using std::string;
using std::setprecision;


// struct geometric 
// {
// 		double prob;

// 		geometric(const double p) : prob(p) {}
// 		double operator()(const double &val) const;
// 		void fit(const vector<size_t> &vals,
// 				 const vector<double> &p);
// 		string tostring() const;
// };

// string
// geometric::tostring() const 
// {
// 		std::ostringstream os;
// 		os << setprecision(4) << prob;
// 		return os.str();
// }

// double 
// geometric::operator()(const double &val) const 
// {
// 		return val * gsl_sf_log(1 - prob) + gsl_sf_log(prob);
// }

// void
// geometric::fit(const vector<size_t> &vals, const vector<double> &p) 
// {
// 		const double p_total = std::accumulate(p.begin(), p.end(), 0.0);
// 		const double val_total = inner_product(vals.begin(), vals.end(), 
// 											   p.begin(), 0.0);
// 		prob = p_total / val_total;
// }

betadistribution::betadistribution(double a, double b)
{
		alpha = a;
		beta = b;
		lnbeta_helper = gsl_sf_lnbeta(a, b);
}

betadistribution& 
betadistribution::set_alpha(double a)
{
		alpha = a;
		lnbeta_helper = gsl_sf_lnbeta(alpha, beta);
		return *this;
}

betadistribution&
betadistribution::set_beta(double b)
{
		beta = b;
		lnbeta_helper = gsl_sf_lnbeta(alpha, beta);
		return *this;
}

double
betadistribution::operator()(double v) const
{
		return (alpha - 1) * gsl_sf_log(v)
				+ (beta -1 ) * gsl_sf_log(1 - v)
				- lnbeta_helper;
}

inline static double 
sign(double x) {
  return (x >= 0) ? 1.0 : -1.0;
}

inline static double
invpsi(const double tolerance, const double x)
{
		double L = 1.0, Y = std::exp(x);
		while (L > tolerance) {
				Y += L*sign(x - gsl_sf_psi(Y));
				L /= 2.0;
		}
		return Y;
}

static double
movement(const double curr, const double prev)
{
		return std::abs(curr - prev)/std::max(std::fabs(curr), std::fabs(prev));
}

void
betadistribution::fit(const vector<double> &vals,
					  const vector<double> &probs)
{
		const double probs_total = std::accumulate(probs.begin(), probs.end(), 0.0);
		vector<double> log_vals(vals.size());
		std::transform(vals.begin(), vals.end(), log_vals.begin(), gsl_sf_log);
		const double alpha_rhs = inner_product(log_vals.begin(), log_vals.end(), 
											   probs.begin(), 0.0)/probs_total;
		for (size_t i = 0; i < log_vals.size(); ++i)
				log_vals[i] = 1 - log_vals[i];
		const double beta_rhs = inner_product(log_vals.begin(), log_vals.end(), 
											  probs.begin(), 0.0)/probs_total;
		double prev_alpha = 0.0, prev_beta = 0.0;
		alpha = beta = 0.01;
		while (movement(alpha, prev_alpha) > tolerance &&
			   movement(beta, prev_beta) > tolerance)
		{
				prev_alpha = alpha;
				prev_beta = beta;
				alpha = invpsi(tolerance, gsl_sf_psi(prev_alpha + prev_beta) + alpha_rhs);
				beta = invpsi(tolerance, gsl_sf_psi(prev_alpha + prev_beta) + beta_rhs);
		}
		lnbeta_helper = gsl_sf_lnbeta(alpha, beta);
}


string
betadistribution::tostring() const
{
		std::ostringstream oss;
		oss << setprecision(4) << alpha << " " << setprecision(4) << beta;
		return oss.str();
}

/***********************/
/********* HMM *********/
/***********************/

inline double
HMM::log_sum_log(const double p, const double q) const {
		if (p == 0) {return q;}
		else if (q == 0) {return p;}
		const double larger = (p > q) ? p : q;
		const double smaller = (p > q) ? q : p;
		return larger + log(1.0 + exp(smaller - larger));
}


double
HMM::log_sum_log_vec(const vector<double> &vals, size_t limit) const {
		const vector<double>::const_iterator x = 
				std::max_element(vals.begin(), vals.begin() + limit);
		const double max_val = *x;
		const size_t max_idx = x - vals.begin();
		double sum = 1.0;
		for (size_t i = 0; i < limit; ++i) {
				if (i != max_idx) {
						sum += exp(vals[i] - max_val);
#ifdef DEBUG
						assert(finite(sum));
#endif
				}
		}
		return max_val + log(sum);
}


double 
HMM::forward_algorithm(const vector<double> &values, 
					   const size_t start, 
					   const size_t end,
					   const vector<double> &log_start_trans,
					   const vector< vector<double> > &log_trans,
					   const vector<double> &log_end_trans,
					   const vector<distro_type> &distros,
					   vector< vector<double> > &f) const
{
// 		f[start].first = fg_distro(vals[start]) + lp_sf;
// 		f[start].second = bg_distro(vals[start]) + lp_sb;
// 		for (size_t i = start + 1; i < end; ++i) {
// 				// assert(finite(fg_distro.log_likelihood(vals[i])));
// 				const size_t k = i - 1;
// 				f[i].first = (fg_distro(vals[i]) +
// 							  log_sum_log(f[k].first + lp_ff, f[k].second + lp_bf));
// 				f[i].second = (bg_distro(vals[i]) + 
// 							   log_sum_log(f[k].first + lp_fb, f[k].second + lp_bb));
// 		}
// 		return log_sum_log(f[end - 1].first + lp_ft, f[end - 1].second + lp_bt);

		for (size_t i = 0; i < distros.size(); ++i)
				f[start][i] =  log_start_trans[i] + distros[i](values[start]);
		for (size_t i = start + 1; i <  end; ++i)
		{
				const size_t k = i - 1;
				for (size_t j = 0; j < distros.size(); ++j)
				{
						// each element is like f[k].first + lp_ff
						vector<double> log_p(distros.size());
						for (size_t ii = 0; ii < distros.size(); ++ii)
								log_p[ii] = f[k][ii] + log_trans[ii][j];
						 
						f[i][j] = distros[j](values[i]) + 
								log_sum_log_vec(log_p, log_p.size());
				}
		}

		vector<double> log_p(distros.size());
		for (size_t ii = 0; ii < distros.size(); ++ii)
				log_p[ii] = f[end - 1][ii] + log_end_trans[ii];
		
		return log_sum_log_vec(log_p, log_p.size());
}

double 
HMM::backward_algorithm(const vector<double> &values, 
						const size_t start, 
						const size_t end,
						const vector<double> &log_start_trans,
						const vector< vector<double> > &log_trans,
						const vector<double> &log_end_trans,
						const vector<distro_type> &distros,
						vector< vector<double> > &b) const
{

		for (size_t i = 0; i < distros.size(); ++i)
				b[end - 1][i] = log_end_trans[i];
		for (size_t k = end - 1; k >  start; --k)
		{
				const size_t i = k - 1;
				for (size_t j = 0; j < distros.size(); ++j)
				{
						// each element is like f[k].first + lp_ff
						vector<double> log_p(distros.size());
						for (size_t ii = 0; ii < distros.size(); ++ii)
								log_p[ii] = b[k][ii] + distros[ii](values[k]) + log_trans[j][ii];
						 
						b[i][j] = log_sum_log_vec(log_p, log_p.size());
				}
		}

		vector<double> log_p(distros.size());
		for (size_t ii = 0; ii < distros.size(); ++ii)
				log_p[ii] = b[start][ii] + distros[ii](values[start]) + log_start_trans[ii];
		
		return log_sum_log_vec(log_p, log_p.size());
}


// double
// HMM::backward_algorithm(const vector<size_t > &vals,
// 								 const size_t start, const size_t end,
// 								 const double lp_sf, const double lp_sb,
// 								 const double lp_ff, const double lp_fb, 
// 								 const double lp_ft,
// 								 const double lp_bf, const double lp_bb, 
// 								 const double lp_bt,
// 								 const geometric &fg_distro,
// 								 const geometric &bg_distro,
// 								 vector<pair<double, double> > &b) const {
// 		b[end - 1].first = lp_ft;
// 		b[end - 1].second = lp_bt;
// 		for (size_t k = end - 1; k > start; --k) {
// 				size_t i = k - 1;
// 				const double fg_a = fg_distro(vals[k]) + b[k].first;
// 				const double bg_a = bg_distro(vals[k]) + b[k].second;
// 				b[i].first = log_sum_log(fg_a + lp_ff, bg_a + lp_fb);
// 				b[i].second = log_sum_log(fg_a + lp_bf, bg_a + lp_bb);
// 		}
// 		return log_sum_log(b[start].first + fg_distro(vals[start]) + lp_sf,
// 						   b[start].second + bg_distro(vals[start]) + lp_sb);
// }


void
HMM::estimate_emissions(const vector< vector<double> > &forward,
						const vector< vector<double> > &backward,
						vector< vector<double> > &probs) const
{
		for (size_t i = 0; i < forward.size(); ++i)
		{
				// the log_likelihood of being state j at time i
				vector<double> log_gamma(probs.size());
				for (size_t j = 0; j < probs.size(); ++j)
						log_gamma[j] = forward[i][j] + backward[i][j];

				const double denom = log_sum_log_vec(log_gamma, log_gamma.size());

				for (size_t j = 0; j < probs.size(); ++j)
						probs[i][j] = exp(log_gamma[j] - denom);
		}
}

void
HMM::estimate_transitions(const vector<value_type> &values,
						  const size_t start,
						  const size_t end,
						  const vector< vector<double> > &forward,
						  const vector< vector<double> > &backward,
						  const double total, 
						  const vector<distro_type> &distros,
						  const vector< vector<double> > &log_trans,
						  const vector<double> &log_end_trans,
						  vector< vector< vector<double> > > &vals) const
{
		for (size_t i = start + 1; i < end; ++i)
		{
				const size_t k = i - 1;
				for (size_t state_now = 0; state_now < distros.size(); ++state_now)
						for (size_t state_next = 0; state_next < distros.size(); ++state_next)
								vals[k][state_now][state_next]
										= forward[k][state_now] +
										log_trans[state_now][state_next] +
										distros[state_next](values[i]) +
										backward[i][state_next] -
										total;
		}
}


// void
// HMM::estimate_transitions(const vector<size_t > &vals,
// 								   const size_t start, const size_t end,
// 								   const vector<pair<double, double> > &f,
// 								   const vector<pair<double, double> > &b,
// 								   const double total,
// 								   const geometric &fg_distro,
// 								   const geometric &bg_distro,
// 								   const double lp_ff, const double lp_fb,
// 								   const double lp_bf, const double lp_bb,
// 								   const double lp_ft, const double lp_bt,
// 								   vector<double> &ff_vals,
// 								   vector<double> &fb_vals,
// 								   vector<double> &bf_vals,
// 								   vector<double> &bb_vals) const
// {
// 		for (size_t i = start + 1; i < end; ++i) {
// 				const size_t k = i - 1;
// 				const double b_first = b[i].first;
// 				const double b_second = b[i].second;
    
// 				const double ff = f[k].first + fg_distro(vals[i]) - total;
// 				const double bb = f[k].second + bg_distro(vals[i]) - total;
    
// 				ff_vals[k] = ff + lp_ff + b_first;
// 				fb_vals[k] = ff + lp_fb + b_second;
    
// 				bf_vals[k] = bb + lp_bf + b_first;
// 				bb_vals[k] = bb + lp_bb + b_second;
// 		}
// }


double 
HMM::single_iteration(const vector<double> &values,
					  const vector<size_t> &reset_points,
					  vector< vector<double> > &forward,
					  vector< vector<double> > &backward,
					  vector<double> &start_trans,
					  vector< vector<double> > &trans,
					  vector<double> &end_trans,
					  vector<distro_type> &distros)
{
// 		vector<double> log_fg_expected;
// 		vector<double> log_bg_expected;
  
		double total_score = 0;
  
		vector<double> log_start_trans(start_trans.size());
		std::transform(start_trans.begin(),
					   start_trans.end(),
					   log_start_trans.begin(),
					   log);

		vector< vector<double> > log_trans(trans.size(), vector<double>(trans.size()));
		for (size_t i = 0; i < trans.size(); ++i)
				std::transform(trans[i].begin(),
							   trans[i].end(),
							   log_trans[i].begin(),
							   log);

		vector<double> log_end_trans(end_trans.size());
		std::transform(end_trans.begin(),
					   end_trans.end(),
					   log_end_trans.begin(),
					   log);

		// for estimating transitions
		vector< vector< vector<double> > >
				vals(trans.size(),
					 vector< vector<double> >(trans.size(),
											  vector<double>(values.size(), 0)));
		
// 		vector<double> ff_vals(values.size(), 0);
// 		vector<double> fb_vals(values.size(), 0);
// 		vector<double> bf_vals(values.size(), 0);
// 		vector<double> bb_vals(values.size(), 0);
  
		for (size_t i = 0; i < reset_points.size() - 1; ++i)
		{
				const double score = forward_algorithm(values, 
													   reset_points[i], 
													   reset_points[i + 1],
													   log_start_trans,
													   log_trans,
													   log_end_trans,
													   distros,
													   forward);
				const double backward_score = 
						backward_algorithm(values, 
										   reset_points[i], 
										   reset_points[i + 1],
										   log_start_trans,
										   log_trans,
										   log_end_trans,
										   distros,
										   backward);								
    
				if (DEBUG && (fabs(score - backward_score)/
							  max(score, backward_score)) > 1e-10)
						cerr << "fabs(score - backward_score)/"
							 << "max(score, backward_score) > 1e-10" << endl;
    
				estimate_transitions(values,
									 reset_points[i], 
									 reset_points[i + 1],
									 forward, backward,
									 score, 
									 distros,
									 log_trans,
									 log_end_trans,
									 vals);

				total_score += score;
		}

		// Subtracting 1 from the limit of the summation because the final
		// term has no meaning since there is no transition to be counted
		// from the final observation (they all must go to terminal state)
// 		const double p_ff_new_estimate = exp(log_sum_log_vec(ff_vals, values.size() - 1));
// 		const double p_fb_new_estimate = exp(log_sum_log_vec(fb_vals, values.size() - 1));
// 		const double p_bf_new_estimate = exp(log_sum_log_vec(bf_vals, values.size() - 1));
// 		const double p_bb_new_estimate = exp(log_sum_log_vec(bb_vals, values.size() - 1));

		vector< vector<double> > trans_new_est(trans.size(), vector<double>(trans.size()));
		for (size_t i = 0; i < trans.size(); ++i)
				for (size_t j = 0; j < trans.size(); ++j)
						trans_new_est[i][j] = exp(log_sum_log_vec(vals[i][j], values.size() - 1));

  
// 		double denom = (p_ff_new_estimate + p_fb_new_estimate);
// 		p_ff = p_ff_new_estimate/denom - p_ft/2.0;
// 		p_fb = p_fb_new_estimate/denom - p_ft/2.0;
  
// 		if (p_ff < MIN_PROB) {
// 				if (DEBUG)
// 						cerr << "p_ff < MIN_PROB" << endl;
// 				p_ff = MIN_PROB;
// 		}
  
// 		if (p_fb < MIN_PROB) {
// 				if (DEBUG)
// 						cerr << "p_fb < MIN_PROB" << endl;
// 				p_fb = MIN_PROB;
// 		}
  
// 		denom = (p_bf_new_estimate + p_bb_new_estimate);
// 		p_bf = p_bf_new_estimate/denom - p_bt/2.0;
// 		p_bb = p_bb_new_estimate/denom - p_bt/2.0;
  
// 		if (p_bf < MIN_PROB) {
// 				if (DEBUG)
// 						cerr << "p_bf < MIN_PROB" << endl;
// 				p_bf = MIN_PROB;
// 		}

// 		if (p_bb < MIN_PROB) {
// 				if (DEBUG)
// 						cerr << "p_bb < MIN_PROB" << endl;
// 				p_bb = MIN_PROB;
// 		}

		for (size_t i = 0; i < trans.size(); ++i)
		{
				double denom = std::accumulate(trans_new_est[i].begin(),
											   trans_new_est[i].end(),
											   0.0);
				for (size_t j = 0; j < trans.size(); ++j)
				{
						trans[i][j] = trans_new_est[i][j] / denom - end_trans[i] / 2.0;
						if (trans[i][j] < MIN_PROB)
						{
								if (DEBUG)
										cerr << "transition probability is smaller than MIN_PROB"
											 << endl;
								trans[i][j] = MIN_PROB;
						}
				}
		}
		
		// for estimating emissions
// 		vector<double> fg_probs(values.size());
// 		vector<double> bg_probs(values.size());
// 		estimate_emissions(forward, backward, fg_probs, bg_probs);
  
// 		fg_distro.fit(values, fg_probs);
// 		bg_distro.fit(values, bg_probs);

		vector< vector<double> > probs(trans.size(), vector<double>(values.size()) );
		estimate_emissions(forward, backward, probs);

		for (size_t i = 0; i < distros.size(); ++i)
				distros[i].fit(values, probs[i]);
		
		return total_score;
}


double
HMM::BaumWelchTraining(const vector<value_type> &values,
				  const vector<size_t> &reset_points,
				  vector<double> &start_trans,
				  vector<vector<double> > &trans, 
				  vector<double> &end_trans,
				  vector<distro_type> &distros) const
{
		assert(distros.size() >= 2);
		assert(start_trans.size() == distros.size());
		assert(end_trans.size() == distros.size());
		assert(trans.size() == distros.size());

		for (size_t i = 0; i < trans.size(); ++i)
				assert(trans[i].size() == distros.size());
  
		vector< vector<double>  > forward(values.size(), vector<double>(distros.size(), 0));
		vector< vector<double>  > backward(values.size(), vector<double>(distros.size(), 0));
  
		if (VERBOSE)
		{
				cerr << "ITR | "
					 << "DELTA | "
					 << "PARAMS | "
					 << "sizes"
					 << endl;
		}
  
		double prev_total = - std::numeric_limits<double>::max();

		for (size_t i = 0; i < max_iterations; ++i)
		{
				vector<double> start_trans_est = start_trans;
				vector< vector<double> > trans_est = trans;
				vector<double> end_trans_est = end_trans;
    
				double total = single_iteration(values, 
												reset_points,
												forward, backward,
												start_trans_est,
												trans_est,
												end_trans_est,
												distros);
    
				if (VERBOSE)
				{
						cerr << i + 1 << " | ";
						cerr << total << "\t"
							 <<	prev_total << "\t"
								(total - prev_total)/std::fabs(total) << " | ";
						for (size_t i = 0; i < distros.size(); ++i)
								cerr << distros[i].tostring() << "\t";
						cerr << "| ";
						for (size_t i = 0; i < trans.size(); ++i)
								cerr << 1 /  trans[i][i] << "\t";
						cerr << endl;
				}

				if ((total - prev_total) < tolerance) {
						// if (fabs((total - prev_total)/total) < tolerance) {
//     if (fabs(total - prev_total)// /abs(max(total, prev_total))
// 	< tolerance) {
						if (VERBOSE)
								cerr << "CONVERGED" << endl << endl;
						break;
				}
    
				start_trans = start_trans_est;
				trans = trans_est;
				end_trans = end_trans_est;

				prev_total = total;
		}
		return prev_total;
}


// void
// HMM::PosteriorScores(const vector<size_t > &values,
// 					 const vector<size_t> &reset_points,
// 					 const vector<double> &start_trans,
// 					 const vector<vector<double> > &trans, 
// 					 const vector<double> &end_trans,
// 					 const double fg_lambda,
// 					 const double bg_lambda,
// 					 const vector<bool> &classes,
// 					 vector<double> &llr_scores) const {

// 		const geometric fg_distro(fg_lambda);
// 		const geometric bg_distro(bg_lambda);

// 		assert(start_trans.size() >= 2);
// 		assert(end_trans.size() >= 2);
// 		assert(trans.size() >= 2);
// 		for (size_t i = 0; i < trans.size(); ++i)
// 				assert(trans[i].size() >= 2);
  
// 		return PosteriorScores(values, reset_points,
// 							   start_trans[0], start_trans[1],
// 							   trans[0][0], trans[0][1], end_trans[0],
// 							   trans[1][0], trans[1][1], end_trans[1],
// 							   fg_distro, bg_distro, classes, llr_scores);
// }



// void
// HMM::PosteriorScores(const vector<size_t > &values,
// 					 const vector<size_t> &reset_points,
// 					 double p_sf, double p_sb,
// 					 double p_ff, double p_fb, double p_ft,
// 					 double p_bf, double p_bb, double p_bt,
// 					 const geometric &fg_distro,
// 					 const geometric &bg_distro,
// 					 const vector<bool> &classes,
// 					 vector<double> &llr_scores) const {

// 		double total_score = 0;
  
// 		const double lp_sf = log(p_sf);
// 		const double lp_sb = log(p_sb);
// 		const double lp_ff = log(p_ff);
// 		const double lp_fb = log(p_fb);
// 		const double lp_ft = log(p_ft);
// 		const double lp_bf = log(p_bf);
// 		const double lp_bb = log(p_bb);
// 		const double lp_bt = log(p_bt);
  
// 		assert(finite(lp_sf) && finite(lp_sb) && 
// 			   finite(lp_ff) && finite(lp_fb) && finite(lp_ft) && 
// 			   finite(lp_bf) && finite(lp_bb) && finite(lp_bt));

// 		vector<pair<double, double> > forward(values.size(), pair<double, double>(0, 0));
// 		vector<pair<double, double> > backward(values.size(), pair<double, double>(0, 0));

// 		for (size_t i = 0; i < reset_points.size() - 1; ++i) {
// 				const double score = forward_algorithm(values, 
// 													   reset_points[i],
// 													   reset_points[i + 1],
// 													   lp_sf, lp_sb,
// 													   lp_ff, lp_fb, lp_ft,
// 													   lp_bf, lp_bb, lp_bt,
// 													   fg_distro, bg_distro, forward);
    
// 				const double backward_score = 
// 						backward_algorithm(values,
// 										   reset_points[i],
// 										   reset_points[i + 1],
// 										   lp_sf, lp_sb,
// 										   lp_ff, lp_fb, lp_ft,
// 										   lp_bf, lp_bb, lp_bt,
// 										   fg_distro, bg_distro, backward);
    
// 				if (DEBUG && (fabs(score - backward_score)/
// 							  max(score, backward_score)) > 1e-10)
// 						cerr << "fabs(score - backward_score)/"
// 							 << "max(score, backward_score) > 1e-10" << endl;

// 				total_score += score;
// 		}
  
// 		llr_scores.resize(values.size());
// 		for (size_t i = 0; i < values.size(); ++i) {
// 				const double fg_state = forward[i].first + backward[i].first;
// 				const double bg_state = forward[i].second + backward[i].second;
// 				if (classes[i])
// 						llr_scores[i] = (fg_state - bg_state);
// 				else 
// 						llr_scores[i] = (bg_state - fg_state);
// 		}
// }



// void
// HMM::PosteriorScores(const vector<size_t > &values,
// 					 const vector<size_t> &reset_points,
// 					 const vector<double> &start_trans,
// 					 const vector<vector<double> > &trans, 
// 					 const vector<double> &end_trans,
// 					 const double fg_lambda,
// 					 const double bg_lambda,
// 					 const bool fg_class,
// 					 vector<double> &llr_scores) const {
  
// 		const geometric fg_distro(fg_lambda);
// 		const geometric bg_distro(bg_lambda);


// 		assert(start_trans.size() >= 2);
// 		assert(end_trans.size() >= 2);
// 		assert(trans.size() >= 2);
// 		for (size_t i = 0; i < trans.size(); ++i)
// 				assert(trans[i].size() >= 2);
  
// 		return PosteriorScores(values, reset_points,
// 							   start_trans[0], start_trans[1],
// 							   trans[0][0], trans[0][1], end_trans[0],
// 							   trans[1][0], trans[1][1], end_trans[1],
// 							   fg_distro, bg_distro, fg_class, llr_scores);
// }




// void
// HMM::PosteriorScores(const vector<size_t > &values,
// 					 const vector<size_t> &reset_points,
// 					 double p_sf, double p_sb,
// 					 double p_ff, double p_fb, double p_ft,
// 					 double p_bf, double p_bb, double p_bt,
// 					 const geometric &fg_distro,
// 					 const geometric &bg_distro,
// 					 const bool fg_class,
// 					 vector<double> &llr_scores) const {
  
// 		double total_score = 0;
  
// 		const double lp_sf = log(p_sf);
// 		const double lp_sb = log(p_sb);
// 		const double lp_ff = log(p_ff);
// 		const double lp_fb = log(p_fb);
// 		const double lp_ft = log(p_ft);
// 		const double lp_bf = log(p_bf);
// 		const double lp_bb = log(p_bb);
// 		const double lp_bt = log(p_bt);
  
// 		assert(finite(lp_sf) && finite(lp_sb) && 
// 			   finite(lp_ff) && finite(lp_fb) && finite(lp_ft) && 
// 			   finite(lp_bf) && finite(lp_bb) && finite(lp_bt));
  
// 		vector<pair<double, double> > forward(values.size(), pair<double, double>(0, 0));
// 		vector<pair<double, double> > backward(values.size(), pair<double, double>(0, 0));
  
// 		for (size_t i = 0; i < reset_points.size() - 1; ++i) {
// 				const double score = 
// 						forward_algorithm(values, reset_points[i], reset_points[i + 1],
// 										  lp_sf, lp_sb, lp_ff, lp_fb, lp_ft, lp_bf, lp_bb, lp_bt,
// 										  fg_distro, bg_distro, forward);
    
// 				const double backward_score = 
// 						backward_algorithm(values, reset_points[i], reset_points[i + 1],
// 										   lp_sf, lp_sb, lp_ff, lp_fb, lp_ft,
// 										   lp_bf, lp_bb, lp_bt,
// 										   fg_distro, bg_distro, backward);
    
// 				if (DEBUG && (fabs(score - backward_score)/
// 							  max(score, backward_score)) > 1e-10)
// 						cerr << "fabs(score - backward_score)/"
// 							 << "max(score, backward_score) > 1e-10" << endl;
// 				total_score += score;
// 		}
  
// 		llr_scores.resize(values.size());
// 		for (size_t i = 0; i < values.size(); ++i) {
// 				const double fg_state = forward[i].first + backward[i].first;
// 				const double bg_state = forward[i].second + backward[i].second;
// 				if (fg_class)
// 						llr_scores[i] = exp(fg_state - log_sum_log(fg_state, bg_state));
// 				else
// 						llr_scores[i] = exp(bg_state - log_sum_log(fg_state, bg_state));
// 				//     if (fg_class)
// 				//       llr_scores[i] = (fg_state - bg_state);
// 				//     else
// 				//       llr_scores[i] = (bg_state - fg_state);
// 		}
// }





// void
// HMM::TransitionPosteriors(const vector<size_t > &values,
// 						  const vector<size_t> &reset_points,
// 						  const vector<double> &start_trans,
// 						  const vector<vector<double> > &trans, 
// 						  const vector<double> &end_trans,
// 						  const double fg_lambda,
// 						  const double bg_lambda,
// 						  const size_t transition,
// 						  vector<double> &llr_scores) const {
  
// 		const geometric fg_distro(fg_lambda);
// 		const geometric bg_distro(bg_lambda);

// 		assert(start_trans.size() >= 2);
// 		assert(end_trans.size() >= 2);
// 		assert(trans.size() >= 2);
// 		for (size_t i = 0; i < trans.size(); ++i)
// 				assert(trans[i].size() >= 2);
  
// 		return TransitionPosteriors(values, reset_points,
// 									start_trans[0], start_trans[1],
// 									trans[0][0], trans[0][1], end_trans[0],
// 									trans[1][0], trans[1][1], end_trans[1],
// 									fg_distro, bg_distro, transition, llr_scores);
// }

void
HMM::TransitionPosteriors(const vector<value_type> &values,
						  const vector<size_t> &reset_points,
						  const vector<double> &start_trans, 
						  const vector<vector<double> > &trans, 
						  const vector<double> &end_trans, 
						  const vector<distro_type> &distros,
						  const size_t state_first,
						  const size_t state_second,
						  vector<double> &scores) const
{
		double total_score = 0;
  
// 		const double lp_sf = log(p_sf);
// 		const double lp_sb = log(p_sb);
// 		const double lp_ff = log(p_ff);
// 		const double lp_fb = log(p_fb);
// 		const double lp_ft = log(p_ft);
// 		const double lp_bf = log(p_bf);
// 		const double lp_bb = log(p_bb);
// 		const double lp_bt = log(p_bt);
		vector<double> log_start_trans(start_trans.size());
		std::transform(start_trans.begin(),
					   start_trans.end(),
					   log_start_trans.begin(),
					   log);

		vector< vector<double> > log_trans(trans.size(), vector<double>(trans.size()));
		for (size_t i = 0; i < trans.size(); ++i)
				std::transform(trans[i].begin(),
							   trans[i].end(),
							   log_trans[i].begin(),
							   log);

		vector<double> log_end_trans(end_trans.size());
		std::transform(end_trans.begin(),
					   end_trans.end(),
					   log_end_trans.begin(),
					   log);
  
// 		assert(finite(lp_sf) && finite(lp_sb) && 
// 			   finite(lp_ff) && finite(lp_fb) && finite(lp_ft) && 
// 			   finite(lp_bf) && finite(lp_bb) && finite(lp_bt));

		vector< vector<double>  > forward(values.size(), vector<double>(distros.size(), 0));
		vector< vector<double>  > backward(values.size(), vector<double>(distros.size(), 0));
  
		for (size_t i = 0; i < reset_points.size() - 1; ++i) {
				const double score = forward_algorithm(values, 
													   reset_points[i], 
													   reset_points[i + 1],
													   log_start_trans,
													   log_trans,
													   log_end_trans,
													   distros,
													   forward);
    
				const double backward_score =
						backward_algorithm(values, 
										   reset_points[i], 
										   reset_points[i + 1],
										   log_start_trans,
										   log_trans,
										   log_end_trans,
										   distros,
										   backward);
				
				if (DEBUG && (fabs(score - backward_score)/
							  max(score, backward_score)) > 1e-10)
						cerr << "fabs(score - backward_score)/"
							 << "max(score, backward_score) > 1e-10" << endl;

				total_score += score;
		}
  
		scores.resize(values.size());
		size_t j = 0;
		for (size_t i = 0; i < values.size(); ++i)
		{
				if (i == reset_points[j])
				{
						++j;
						scores[i] = 0;
				} else {
// 						const double fg_to_fg_state = forward[i - 1].first + lp_ff + // transition
// 								// emission for value i + 1
// 								fg_distro(values[i]) + backward[i].first;
// 						const double fg_to_bg_state = forward[i - 1].first + lp_fb + 
// 								bg_distro(values[i]) + backward[i].second;
// 						const double bg_to_fg_state = forward[i - 1].second + lp_bf + 
// 								fg_distro(values[i]) + backward[i].first;
// 						const double bg_to_bg_state = forward[i - 1].second + lp_bb + 
// 								bg_distro(values[i]) + backward[i].second;
// 						const double denom = log_sum_log(log_sum_log(fg_to_fg_state, fg_to_bg_state),
// 														 log_sum_log(bg_to_fg_state, bg_to_bg_state));
// 						double numerator = fg_to_fg_state;
// 						if (transition == 1)
// 								numerator = fg_to_bg_state;
// 						if (transition == 2)
// 								numerator = bg_to_fg_state;
// 						if (transition == 3)
// 								numerator = bg_to_bg_state;
// 						scores[i] = exp(numerator - denom);
						double denom = 0;
						double numerator = 0;
						
						for (size_t j = 0; j < distros.size(); ++j)
						{
								const double log_emission_prob = distros[j](values[i]);
								for (size_t k = 0; k < distros.size(); ++k)
								{
										const double sts = forward[i - 1][k]
												+ log_trans[k][j]
												+ log_emission_prob
												+ backward[i][j];
										denom = log_sum_log(denom, sts);
										
										if (k == state_first && j == state_second)
												numerator = sts;
								}
						}
						scores[i] = exp(numerator - denom);
				}
		}
}


// double
// HMM::PosteriorDecoding(const vector<size_t > &values,
// 					   const vector<size_t> &reset_points,
// 					   const vector<double> &start_trans,
// 					   const vector<vector<double> > &trans, 
// 					   const vector<double> &end_trans,
// 					   const double fg_lambda,
// 					   const double bg_lambda,
// 					   vector<bool> &classes,
// 					   vector<double> &llr_scores) const {

  
// 		const geometric fg_distro(fg_lambda);
// 		const geometric bg_distro(bg_lambda);

// 		assert(start_trans.size() >= 2);
// 		assert(end_trans.size() >= 2);
// 		assert(trans.size() >= 2);
// 		for (size_t i = 0; i < trans.size(); ++i)
// 				assert(trans[i].size() >= 2);
  
// 		return PosteriorDecoding(values, reset_points,
// 								 start_trans[0], start_trans[1],
// 								 trans[0][0], trans[0][1], end_trans[0],
// 								 trans[1][0], trans[1][1], end_trans[1],
// 								 fg_distro, bg_distro, classes, llr_scores);
// }


double
HMM::PosteriorDecoding(const vector<value_type> &values,
					   const vector<size_t> &reset_points,
					   const vector<double> &start_trans,
					   const vector< vector<double> > &trans,
					   const vector<double> &end_trans,
					   const vector<distro_type> &distros,
					   vector<state_type> &classes,
					   vector<double> &llr_scores) const 
{
		assert(distros.size() >= 2);
		assert(start_trans.size() == distros.size());
		assert(end_trans.size() == distros.size());
		assert(trans.size() == distros.size());
		for (size_t i = 0; i < trans.size(); ++i)
				assert(trans[i].size() == distros.size());

		double total_score = 0;
  
// 		const double lp_sf = log(p_sf);
// 		const double lp_sb = log(p_sb);
// 		const double lp_ff = log(p_ff);
// 		const double lp_fb = log(p_fb);
// 		const double lp_ft = log(p_ft);
// 		const double lp_bf = log(p_bf);
// 		const double lp_bb = log(p_bb);
// 		const double lp_bt = log(p_bt);
  
		vector<double> log_start_trans(start_trans.size());
		std::transform(start_trans.begin(),
					   start_trans.end(),
					   log_start_trans.begin(),
					   log);

		vector< vector<double> > log_trans(trans.size(), vector<double>(trans.size()));
		for (size_t i = 0; i < trans.size(); ++i)
				std::transform(trans[i].begin(),
							   trans[i].end(),
							   log_trans[i].begin(),
							   log);

		vector<double> log_end_trans(end_trans.size());
		std::transform(end_trans.begin(),
					   end_trans.end(),
					   log_end_trans.begin(),
					   log);
  
// 		assert(finite(lp_sf) && finite(lp_sb) && 
// 			   finite(lp_ff) && finite(lp_fb) && finite(lp_ft) && 
// 			   finite(lp_bf) && finite(lp_bb) && finite(lp_bt));

// 		vector<pair<double, double> > forward(values.size(), pair<double, double>(0, 0));
// 		vector<pair<double, double> > backward(values.size(), pair<double, double>(0, 0));

		vector< vector<double>  > forward(values.size(), vector<double>(distros.size(), 0));
		vector< vector<double>  > backward (values.size(), vector<double>(distros.size(), 0));


		for (size_t i = 0; i < reset_points.size() - 1; ++i)
		{
				const double score = forward_algorithm(values, 
													   reset_points[i], 
													   reset_points[i + 1],
													   log_start_trans,
													   log_trans,
													   log_end_trans,
													   distros,
													   forward);
				const double backward_score = 
						backward_algorithm(values, 
										   reset_points[i], 
										   reset_points[i + 1],
										   log_start_trans,
										   log_trans,
										   log_end_trans,
										   distros,
										   backward);								
    
				if (DEBUG && (fabs(score - backward_score)/
							  max(score, backward_score)) > 1e-10)
						cerr << "fabs(score - backward_score)/"
							 << "max(score, backward_score) > 1e-10" << endl;

				total_score += score;
		}
  
		classes.resize(values.size());
		llr_scores.resize(values.size());
		for (size_t i = 0; i < values.size(); ++i)
		{
				vector<double> state_posteriors(distros.size());
//				std::transform(forward[i].begin(), forward[i].end(),
//							   backward[i].begin(),
//							   state_posteriors.begin(),
//							   op_sum);
				for (size_t j = 0; j < distros.size(); ++j)
						state_posteriors[j] = forward[i][j] + backward[i][j]; 

				classes[i] = static_cast<state_type>(std::max_element(state_posteriors.begin(),
																	  state_posteriors.end())
													 - state_posteriors.begin());
				const double sum = log_sum_log_vec(state_posteriors, state_posteriors.size());
				llr_scores[i] = exp(state_posteriors[classes[i]] - sum);
				
// 				const double fg_state = forward[i].first + backward[i].first;
// 				const double bg_state = forward[i].second + backward[i].second;
// 				classes[i] = static_cast<bool>(fg_state > bg_state);
// 				llr_scores[i] = (fg_state - bg_state);
		}
  
		return total_score;
}


/*************************************************************
 *
 * Functions for Viterbi training and decoding.
 *
 *************************************************************/


// double
// HMM::ViterbiDecoding(const vector<size_t > &values,
// 					 const vector<size_t> &reset_points,
// 					 const vector<double> &start_trans,
// 					 const vector<vector<double> > &trans, 
// 					 const vector<double> &end_trans,
// 					 const double fg_lambda,
// 					 const double bg_lambda,
// 					 vector<bool> &classes) const {
  
// 		const geometric fg_distro(fg_lambda);
// 		const geometric bg_distro(bg_lambda);
  
// 		assert(start_trans.size() >= 2);
// 		assert(end_trans.size() >= 2);
// 		assert(trans.size() >= 2);
// 		for (size_t i = 0; i < trans.size(); ++i)
// 				assert(trans[i].size() >= 2);
  
// 		return ViterbiDecoding(values, reset_points,
// 							   start_trans[0], start_trans[1],
// 							   trans[0][0], trans[0][1], end_trans[0],
// 							   trans[1][0], trans[1][1], end_trans[1],
// 							   fg_distro, bg_distro, classes);
// }


double
HMM::ViterbiDecoding(const vector<value_type> &values,
					 const vector<size_t> &reset_points,
					 const vector<double> &start_trans,
					 const vector< vector<double> > &trans,
					 const vector<double> &end_trans,
					 const vector<distro_type> &distros,
					 vector<state_type> &ml_classes) const 
{
// 		const double lp_sf = log(p_sf);
// 		const double lp_sb = log(p_sb);
// 		const double lp_ff = log(p_ff);
// 		const double lp_fb = log(p_fb);
// 		const double lp_ft = log(p_ft);
// 		const double lp_bf = log(p_bf);
// 		const double lp_bb = log(p_bb);
// 		const double lp_bt = log(p_bt);
  		vector<double> log_start_trans(start_trans.size());
		std::transform(start_trans.begin(),
					   start_trans.end(),
					   log_start_trans.begin(),
					   log);

		vector< vector<double> > log_trans(trans.size(), vector<double>(trans.size()));
		for (size_t i = 0; i < trans.size(); ++i)
				std::transform(trans[i].begin(),
							   trans[i].end(),
							   log_trans[i].begin(),
							   log);

		vector<double> log_end_trans(end_trans.size());
		std::transform(end_trans.begin(),
					   end_trans.end(),
					   log_end_trans.begin(),
					   log);

		
		// ml_classes = vector<bool>(values.size());
		double total = 0;
		for (size_t i = 0; i < reset_points.size() - 1; ++i)
		{
				const size_t lim = reset_points[i + 1] - reset_points[i];
    
				vector< vector<double> > delta(lim, vector<double>(distros.size(), 0));
				vector< vector<size_t> > trace(lim, vector<size_t>(distros.size(), 0));
				
// 				v.front().first = fg_distro(values[reset_points[i]]) + lp_sf;
// 				v.front().second = bg_distro(values[reset_points[i]]) + lp_sb;
				
				// beginning of sequence
				for (size_t j = 0; j < distros.size(); ++j)
						delta[0][j] = log_end_trans[j] + distros[j](values[reset_points[i]]);

				// middle of sequence
                for (size_t j = 1; j < lim; ++j)
				{
						for (size_t k = 0; k < distros.size(); ++k)
						{
								// dynamic programming
								double max_prob = - std::numeric_limits<double>::max();
								for (size_t state_prev = 0; state_prev < distros.size(); ++state_prev)
								{
										const double log_trans_prob
												= delta[j - 1][state_prev] + log_trans[state_prev][k];
										if (log_trans_prob > max_prob)
										{
												max_prob = log_trans_prob;
												trace[j][k] = state_prev;
										}
								}

								const double log_emission_prob = distros[k](values[reset_points[i] + j]);
								delta[j][k] = max_prob + log_emission_prob;
						}
				}

				// end of sequence
				for (size_t j = 0; j < distros.size(); ++j)
						delta.back()[j] += log_end_trans[j];

// 				v.back().first += lp_ft;
// 				v.back().second += lp_bt;
					 
				vector<state_type> inner_ml_classes;
					 
				// do the traceback
				size_t opt_state = static_cast<size_t>(std::max_element(delta.back().begin(),
																		delta.back().end())
													   - delta.back().begin());
				inner_ml_classes.push_back(opt_state);
    
				for (size_t j = trace.size() - 1; j > 0; --j)
				{
// 						const size_t k = j - 1;
// 						if (prev == 0)
// 						{
// 								inner_ml_classes.push_back(true);
// 								prev = trace[k].first;
// 						} else {
// 								inner_ml_classes.push_back(false);
// 								prev = trace[k].second;
// 						}
						opt_state = trace[j][opt_state];
						inner_ml_classes.push_back(opt_state);
						
				}
				reverse(inner_ml_classes.begin(), inner_ml_classes.end());
				ml_classes.insert(ml_classes.end(), inner_ml_classes.begin(), 
								  inner_ml_classes.end());
					 
				total += *std::max_element(delta.back().begin(), delta.back().end());
		}
		return total;
}
		
