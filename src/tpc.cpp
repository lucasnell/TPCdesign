

#include "tpc.h"


using namespace Rcpp;




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




//' Thermal optimum for Brière-2 thermal performance curve (TPC)
//'
//' Note that parameter `a` is not necessary for this calculation.
//'
//'
//' @param ctmin Numeric vector for parameter `ctmin`.
//' @param ctmax Numeric vector for parameter `ctmax`.
//' @param b Numeric vector for parameter `b`.
//'
//' @returns A numeric vector of optimum temperatures.
//'
//' @export
//'
//[[Rcpp::export]]
arma::vec briere2_tpc_Topt(arma::vec ctmin,
                           arma::vec ctmax,
                           arma::vec b) {

    uint32 n = ctmin.n_elem;
    if (ctmax.n_elem > n) n = ctmax.n_elem;
    if (b.n_elem > n) n = b.n_elem;

    if (n > 1) {
        if (ctmin.n_elem == 1U) {
            double tmp = ctmin.at(0);
            ctmin.set_size(n);
            ctmin.fill(tmp);
        }
        if (ctmax.n_elem == 1U) {
            double tmp = ctmax.at(0);
            ctmax.set_size(n);
            ctmax.fill(tmp);
        }
        if (b.n_elem == 1U) {
            double tmp = b.at(0);
            b.set_size(n);
            b.fill(tmp);
        }
    }

    if (ctmin.n_elem != n) stop("length(ctmin) must be 1 or equal to lengths of b and ctmax");
    if (ctmax.n_elem != n) stop("length(ctmax) must be 1 or equal to lengths of ctmin and b");
    if (b.n_elem != n) stop("length(b) must be 1 or equal to lengths of ctmin and ctmax");

    arma::vec out(n, arma::fill::none);

    double ccbc2;
    for (uint32 i = 0; i < n; i++) {
        const double& ctmin_i(ctmin.at(i));
        const double& ctmax_i(ctmax.at(i));
        const double& b_i(b.at(i));
        ccbc2 = -2 * ctmax_i - ctmin_i - b_i * ctmin_i;
        out.at(i) = (2 * ctmax_i + ctmin_i + b_i * ctmin_i +
            std::sqrt(-4 * (2 + b_i) * ctmax_i * ctmin_i + ccbc2 * ccbc2)) /
                (2 * (2 + b_i));
        // (2 ctmax + ctmin + b ctmin +
        //     Sqrt[-4 (2 + b) ctmax ctmin + (-2 ctmax - ctmin - b ctmin)^2]) /
        // (2 (2 + b))
    }


    return out;
}

