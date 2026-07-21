

#include "util.h"
#include <random>


using namespace Rcpp;


//' Logit and inverse logit functions.
//'
//'@param p Numeric vector of proportion(s)
//' @export
//'
//[[Rcpp::export]]
NumericVector logit(NumericVector p) {
    NumericVector out(p.size());
    for (uint32 i = 0; i < p.size(); i++) {
        logit_cpp(p[i], out[i]);
    }
    return out;
}
//' @rdname logit
//'
//'@param a Numeric vector of values that aren't necessarily proportions
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector inv_logit(NumericVector a){
    NumericVector out(a.size());
    for (uint32 i = 0; i < a.size(); i++) {
        inv_logit_cpp(a[i], out[i]);
    }
    return out;
}



// Truncated normal for testing
//[[Rcpp::export]]
arma::vec trunc_rnorm(const uint32& n, const double& mu, const double& sigma) {

    trunc_normal_distribution distr(mu, sigma);

    pcg32 eng;
    seed_pcg(eng);

    arma::vec out(n, arma::fill::none);
    for (double & x : out) x = distr(eng);

    return out;

}
// Same as above, but it only returns the min and max
//[[Rcpp::export]]
arma::vec trunc_rnorm_range(const uint32& n, const double& mu, const double& sigma) {

    trunc_normal_distribution distr(mu, sigma);

    pcg32 eng;
    seed_pcg(eng);

    arma::vec out = {100, -100};
    if (mu > 100) out.at(0) = mu;
    double x;
    for (uint32 i = 0; i < n; i++) {
        x = distr(eng);
        if (x < out.at(0)) out.at(0) = x;
        if (x > out.at(1)) out.at(1) = x;
    }

    return out;

}



/*
===============================================================================*
===============================================================================*
 Brière-2 thermal performance curve (TPC)
===============================================================================*
===============================================================================*
 */




//' Brière-2 thermal performance curve (TPC)
//'
//' Note that this TPC does not perform well when ctmin is positive, but
//' there are negative temperatures.
//' For example, if you run `briere2_tpc(c(-10, 0, 10), ctmin = 5,
//' ctmax = 30, a = 1, b = 0.2, TRUE)`, you'll see that a temperature of `-10`
//' gives a greater performance value than `10`.
//'
//' @param temp Numeric vector of temperatures
//' @param ctmin Single numeric for parameter `ctmin`.
//' @param ctmax Single numeric for parameter `ctmax`.
//' @param a Single numeric for parameter `a`.
//' @param b Single numeric for parameter `b`.
//' @param scale Single logical for whether to scale to make max value 1.
//'     Defaults to `FALSE`.
//'
//' @returns A numeric vector for measure of performance for each in `temp`
//'
//' @export
//'
//[[Rcpp::export]]
arma::vec briere2_tpc(const arma::vec& temp,
                          const double& ctmin,
                          const double& ctmax,
                          const double& a,
                          const double& b,
                          const bool& scale = false) {
    arma::vec out = briere2_tpc_cpp(temp, ctmin, ctmax, a, b, scale);
    return out;
}





//' Derivative of Brière-2 thermal performance curve (TPC) with respect to time
//'
//'
//' @inheritParams briere2_tpc
//'
//' @returns A numeric vector for first derivative of measures of
//' performance for each in `temp`
//'
//' @export
//'
//[[Rcpp::export]]
arma::vec briere2_tpc_deriv(const arma::vec& temp,
                            const double& ctmin,
                            const double& ctmax,
                            const double& a,
                            const double& b) {

    arma::vec out(temp.n_elem, arma::fill::none);
    double a_ctmax_T_b, a_ctmax_T_b1, temp_ctmin;
    for (uint32 i = 0; i < temp.n_elem; i++) {
        const double& T(temp.at(i));
        if (T >= ctmax || T <= ctmin) {
            out.at(i) = 0;
        } else {
            a_ctmax_T_b = a * std::pow(ctmax - T, b); // a (ctmax - T)^b
            a_ctmax_T_b1 = a_ctmax_T_b / (ctmax - T); // a (ctmax - T)^(b-1)
            temp_ctmin = T - ctmin;
            out.at(i) = a_ctmax_T_b * T + a_ctmax_T_b * temp_ctmin -
                b * a_ctmax_T_b1 * T * temp_ctmin;
        }
    }
    /*
     a (ctmax - T)^b T + a (ctmax - T)^b (T - ctmin) -
        a b (ctmax - T)^(b - 1) T (T - ctmin)
     */

    return out;
}
















/*
===============================================================================*
===============================================================================*
 Simulate data
===============================================================================*
===============================================================================*
 */


//' Simulate data with gamma distributed error (or normal approximation)
//'
//' @param temp Numeric vector giving the observed temperatures for each
//'     performance value.
//' @param n_reps Single integer for the number of experimental reps per
//'     temperature.
//' @param obs_cv Single numeric for the observation error coefficient of variation.
//'     Must be >= 0.
//'     Observation error is simulated using a gamma distribution, except for
//'     when performance values would fall below zero, in which case it uses
//'     a truncated normal approximation.
//' @inheritParams briere2_tpc
//' @param scale_tpc Single logical for whether to scale to make max value 1.
//'     Defaults to `FALSE`.
//'
//' @export
//'
//' @returns A dataframe of temperatures (column `"temp"`) and
//'     performance values (column `"y"`) with gamma (or normal approximation)
//'     error.
//'
//[[Rcpp::export]]
DataFrame sim_gamma_data(const arma::vec& temp,
                         const int& n_reps,
                         const double& obs_cv,
                         const double& ctmin,
                         const double& ctmax,
                         const double& a,
                         const double& b,
                         const bool& scale_tpc = false) {

    if (n_reps < 1) stop("n_reps must be >= 1");
    if (obs_cv <= 0) stop("obs_cv must be > 0");
    // if (ctmin >= ctmax) stop("ctmin must be < ctmax");
    if (a < 0) stop("a must be >= 0");
    if (b < 0) stop("b must be >= 0");

    pcg32 eng;
    seed_pcg(eng);
    trunc_normal_distribution tnorm_d;
    typedef std::gamma_distribution<double> gamma_distribution;
    typedef std::gamma_distribution<double>::param_type gamma_params;
    gamma_distribution gamma_d;

    arma::vec true_y = briere2_tpc_cpp(temp, ctmin, ctmax, a, b, scale_tpc);

    uint32 n_temps = temp.n_elem;
    uint32 n_rows = (uint32)n_reps * n_temps;

    arma::vec out_temp(n_rows, arma::fill::none);
    arma::vec out_y(n_rows, arma::fill::none);

    double obs_cv2 = obs_cv * obs_cv;
    double shape = 1 / obs_cv2;
    double scale, mu, sigma;

    uint32 k = 0;
    for (uint32 i = 0; i < n_temps; i++) {
        scale = true_y.at(i) * obs_cv2;
        if (scale <= 0) {
            // approximate with truncated normal if scale <= 0:
            mu = shape * scale;
            sigma = std::sqrt(shape * scale * scale);
            tnorm_d.reset(mu, sigma);
            for (uint32 j = 0; j < n_reps; j++) {
                out_temp.at(k) = temp.at(i);
                out_y.at(k) = tnorm_d(eng);
                k++;
            }
        } else {
            // otherwise, use gamma:
            gamma_d.param(gamma_params(shape, scale));
            for (uint32 j = 0; j < n_reps; j++) {
                out_temp.at(k) = temp.at(i);
                out_y.at(k) = R::rgamma(shape, scale);
                k++;
            }
        }
    }

    DataFrame out = DataFrame::create(_["temp"] = out_temp,
                                      _["y"] = out_y);

    return out;

}
