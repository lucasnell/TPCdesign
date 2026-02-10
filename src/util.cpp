
#include <RcppArmadillo.h>
#include <algorithm>
#include <math.h>
#include "TPCdesign_types.h"


using namespace Rcpp;






inline void logit_cpp(const double& p, double& out) {
    out = std::log(p / (1-p));
    return;
}
inline void inv_logit_cpp(const double& a, double& out) {
    out = 1 / (1 + std::exp(-a));
    return;
}

inline double logit_cpp(const double& p) {
    double a = std::log(p / (1-p));
    return a;
}
inline double inv_logit_cpp(const double& a) {
    double p = 1 / (1 + std::exp(-a));
    return p;
}



//' Logit and inverse logit functions.
//'
//'
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
//' @describeIn logit
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




/*
===============================================================================*
===============================================================================*
 Brière-2 thermal performance curve (TPC)
===============================================================================*
===============================================================================*
 */


inline NumericVector briere2_tpc_cpp(NumericVector temp,
                                     const double& ctmin,
                                     const double& ctmax,
                                     const double& a,
                                     const double& b,
                                     const bool& scale) {

    NumericVector out(temp.size());
    double out_max = 0;
    for (uint32 i = 0; i < temp.size(); i++) {
        out[i] = a * temp[i] * (temp[i] - ctmin) *
            std::pow(std::max(ctmax - temp[i], 0.0), b);
        if (out[i] > out_max) out_max = out[i];
    }
    if (scale) {
        if (out_max <= 0) out_max = 1;
        for (double& x : out) x /= out_max;
    }
    return out;
}




//' Brière-2 thermal performance curve (TPC)
//'
//' Note that this TPC does not perform well when ctmin is positive, but
//' there are negative temperatures.
//' For example, if you run `briere2_tpc(c(-10, 0, 10), ctmin = 5,
//' ctmax = 30, a = 1, b = 0.2, TRUE)`, you'll see that a temperature of `-10`
//' gives a greater performance value than `10`.
//'
//' @param temp Numeric vector of temperatures
//' @param ctmin Single numeric for parameter ctmin
//' @param ctmax Single numeric for parameter ctmax
//' @param a Single numeric for parameter a
//' @param b Single numeric for parameter b
//' @param scale Single logical for whether to scale to make max value 1.
//'
//' @returns A numeric vector for measure of performance for each in `temp`
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector briere2_tpc(NumericVector temp,
                          const double& ctmin,
                          const double& ctmax,
                          const double& a,
                          const double& b,
                          const bool& scale) {
    NumericVector out = briere2_tpc_cpp(temp, ctmin, ctmax, a, b, scale);
    return out;
}







/*
 ===============================================================================*
 ===============================================================================*
 Make vector of temperatures for input of parameters
 ===============================================================================*
 ===============================================================================*
 */

//' Make vector of temperatures for input of parameters
//'
//' @param params Numeric vector of parameters.
//'     The length of this vector must be equal to the number of desired
//'     temperatures plus 1. A vector of length `<2` returns an error.
//'     The last item in this vector must be the logit-transformed proportion
//'     the last temperature is from the midpoint
//'     (i.e., `mean(c(temp_min, temp_max))`)
//'     to the maximum possible temperature (i.e., `temp_max`).
//'     The first item is the logit-transformed weight affecting
//'     the difference between `temp_min` and the first sampled temperature.
//'     The remaining items are the logit-transformed weights affecting
//'     the distance between the remaining sampled temperatures and
//'     the previous ones.
//' @inheritParams design_temps
//'
//' @returns A numeric vector for the temperatures to sample.
//'
//' @export
//'
//'
//[[Rcpp::export]]
NumericVector make_temps(NumericVector params,
                         const double& temp_min,
                         const double& temp_max) {

    if (params.size() < 2U) stop("params must be of length 2 or longer");

    uint32 n_temps = params.size() - 1;

    double temp_diff = temp_max - temp_min;
    double temp_end = 0.5 * (1 + inv_logit_cpp(params[n_temps])) * temp_diff + temp_min;
    temp_diff = temp_end - temp_min;

    NumericVector cs_probs(n_temps);
    for (uint32 i = 0; i < n_temps; i++) {
        cs_probs[i] = inv_logit_cpp(params[i]);
        if (i > 0) cs_probs[i] += cs_probs[i-1U];
    }

    NumericVector temp(n_temps);
    double denom = (cs_probs[n_temps-1U] > 0) ? cs_probs[n_temps-1U] : 1;
    for (uint32 i = 0; i < n_temps; i++) {
        temp[i] = (cs_probs[i] / denom) * temp_diff + temp_min;
    }

    return temp;

}




/*
 ===============================================================================*
 ===============================================================================*
 Objective function
 ===============================================================================*
 ===============================================================================*
 */


//' Objective function with RMSE returned
//'
//' @param params Numeric vector of length 4 containing (in order)
//'     untransformed `ctmin`, untransformed `ctmax`, log-transformed `a`,
//'     and log-transformed `b`.
//' @param y Numeric vector giving the observed performance values for each
//'     temperature. It's recommended to scale this vector to have a max value
//'     of 1 by dividing by its max value.
//'     Must be the same length as `temp`.
//' @param temp Numeric vector giving the observed temperatures for each
//'     performance value. Must be the same length as `y`.
//' @param scale_tpc Single logical for whether to scale TPC to make max value 1.
//'
//' @returns A single numeric giving the RMSE between observed performance
//'     values and those predicted based on the input parameters.
//'
//' @noRd
//'
//'
//[[Rcpp::export]]
double rmse_objective(NumericVector params,
                      NumericVector y,
                      NumericVector temp,
                      const bool& scale_tpc) {

    if (params.size() != 4) stop("params must be length 4");
    if (y.size() != temp.size()) stop("y must be same length as temp");

    double ctmin = params[0];
    double ctmax = params[1];
    double a = std::exp(params[2]);
    double b = std::exp(params[3]);

    NumericVector y_pred = briere2_tpc_cpp(temp, ctmin, ctmax, a, b, scale_tpc);

    double out = 0;
    for (uint32 i = 0; i < y.size(); i++) {
        out += ((y[i] - y_pred[i]) * (y[i] - y_pred[i]));
    }
    out /= static_cast<double>(y.size());
    out = std::sqrt(out);
    if (Rcpp::traits::is_infinite<REALSXP>(out) ||
        NumericVector::is_na(out)) out = 1e10;

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
//' @inheritParams design_temps
//'
//' @export
//'
//' @returns A dataframe of temperatures (column `"temp"`) and
//'     performance values (column `"y"`) with gamma (or normal approximation)
//'     error.
//'
//[[Rcpp::export]]
DataFrame sim_gamma_data(NumericVector temp,
                         const int& n_reps,
                         const double& obs_cv,
                         const double& ctmin,
                         const double& ctmax,
                         const double& a,
                         const double& b,
                         const bool& scale_tpc = true) {

    if (n_reps < 1) stop("n_reps must be >= 1");
    if (obs_cv < 0) stop("obs_cv must be >= 0");
    // if (ctmin >= ctmax) stop("ctmin must be < ctmax");
    if (a < 0) stop("a must be >= 0");
    if (b < 0) stop("b must be >= 0");

    NumericVector true_y = briere2_tpc_cpp(temp, ctmin, ctmax, a, b, scale_tpc);

    uint32 n_temps = temp.size();
    uint32 n_rows = (uint32)n_reps * n_temps;

    NumericVector out_temp(n_rows);
    NumericVector out_y(n_rows);

    double obs_cv2 = obs_cv * obs_cv;
    double shape = 1 / obs_cv2;
    double scale, mu, sigma;

    uint32 k = 0;
    for (uint32 i = 0; i < n_temps; i++) {
        scale = true_y[i] * obs_cv2;
        for (uint32 j = 0; j < n_reps; j++) {
            out_temp[k] = temp[i];
            if (scale <= 0) {
                // approximate with normal if scale <= 0:
                mu = shape * scale;
                sigma = std::sqrt(shape * scale * scale);
                out_y[k] = R::rnorm(mu, sigma);
            } else {
                // otherwise, use gamma:
                out_y[k] = R::rgamma(shape, scale);
            }
            k++;
        }
    }

    DataFrame out = DataFrame::create(_["temp"] = out_temp,
                                      _["y"] = out_y);

    return out;

}
