
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




//' Brière-2 thermal performance curve (TPC)
//'
//' @param temp Numeric vector of temperatures
//' @param CTmin Single numeric for parameter CTmin
//' @param CTmax Single numeric for parameter CTmax
//' @param a Single numeric for parameter a
//' @param b Single numeric for parameter b
//' @param scale Single logical for whether to scale to make max value 1.
//'     Defaults to `TRUE`
//'
//' @returns A numeric vector for measure of performance for each in `temp`
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector briere2_tpc(NumericVector temp,
                          const double& CTmin,
                          const double& CTmax,
                          const double& a,
                          const double& b,
                          const bool& scale = true) {

    NumericVector out(temp.size());
    double out_max = 0;
    for (uint32 i = 0; i < temp.size(); i++) {
        out[i] = a * temp[i] * std::max(temp[i] - CTmin, 0.0) *
            std::pow(std::max(CTmax - temp[i], 0.0), b);
        if (out[i] > out_max) out_max = out[i];
    }
    if (scale) {
        if (out_max <= 0) out_max = 1;
        for (double& x : out) x /= out_max;
    }
    return out;
}



//' Make vector of temperatures for input of parameters
//'
//' @param temp Numeric vector of parameters. The first must be the
//'     logit-transformed proportion the first temperature to be sampled is
//'     from `temp_min` to `(temp_min + temp_max) / 2`.
//'     The second is the same for the proportion the last temperature to
//'     be sampled is from the midpoint to `temp_max`.
//'     The third and fourth (if present) are the `log(shape1)` and
//'     `log(shape2)`, respectively, for a beta distribution.
//'     This distribution is used to generate quantiles to make the temperatures
//'     not evenly distributed through the temperature space.
//' @param temp_min Single numeric for the minimum possible temperature allowed.
//' @param temp_max Single numeric for the maximum possible temperature allowed.
//' @param n_temps Single integer for the number of temperatures to sample.
//'
//' @returns A numeric vector for the temperatures to sample.
//'
//'
//[[Rcpp::export]]
NumericVector make_temps(NumericVector params,
                         const double& temp_min,
                         const double& temp_max,
                         const int& n_temps) {

    if (params.size() < 2U) stop("params must be of length 2 or longer");
    if (params.size() % 2U != 0U) stop("params must be a length divisible by 2");
    if (n_temps < 1) stop("n_temps must be at least 1");

    // max possible for temp_begin is just before midway between min and max temps
    // min possible for temp_end is just after midway between min and max temps
    // this is so that they do not affect each other but cannot overlap
    double temp_begin = inv_logit_cpp(params[0]) * (temp_max - temp_min) / 2 + temp_min;
    double temp_end = (inv_logit_cpp(params[1]) * (temp_max - temp_min) +
                       (temp_max + temp_min)) / 2;

    NumericVector temp_ps(n_temps);
    if (n_temps > 1) {
        double d = 1 / static_cast<double>(n_temps - (int)1);
        for (uint32 i = 0; i < (uint32)n_temps; i++) temp_ps[i] = d * (double)i;
    } else temp_ps[0] = 0.5;


    if (params.size() == 4) {
        double s1 = std::exp(params[2]);
        double s2 = std::exp(params[3]);
        for (double& p : temp_ps) p = R::qbeta(p, s1, s2, true, false);
    } else if (params.size() > 4) {
        // n_distrs <- (length(params) - 2L) %/% 2L
        // distr_list <- lapply(1:n_distrs, function(i) {
        //     s1 <- exp(params[[i * 2L + 1L]])
        //     s2 <- exp(params[[i * 2L + 2L]])
        //     dst_beta(s1, s2)
        // })
        // mix_distr <- mix(distr_list)
        // temp_ps <- eval_quantile(mix_distr, temp_ps)
        warning("more than one distribution not supported");
    }

    NumericVector design_temps = temp_begin + temp_ps * (temp_end - temp_begin);

    return design_temps;

}



//[[Rcpp::export]]
NumericVector make_temps2(NumericVector params,
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

    NumericVector temps(n_temps);
    for (uint32 i = 0; i < n_temps; i++) {
        temps[i] = (cs_probs[i] / cs_probs[n_temps-1U]) * temp_diff + temp_min;
    }

    return temps;

}
