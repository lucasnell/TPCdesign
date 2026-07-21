
#include "TPCdesign_types.h"

#include <algorithm>
#include <math.h>


using namespace Rcpp;




//[[Rcpp::export]]
double utility_briere2D_cpp(const arma::mat& d, SEXP B) {

    // Because pace passes a list sometimes:
    arma::mat theta;
    if (TYPEOF(B) == VECSXP) {
        List lB(B);
        NumericMatrix mB(lB[0]);
        theta = arma::mat(mB.begin(), mB.nrow(), mB.ncol(), false);
    } else if (TYPEOF(B) == REALSXP) {
        NumericMatrix mB(B);
        theta = arma::mat(mB.begin(), mB.nrow(), mB.ncol(), false);
    } else stop("B must be list or numeric");

    const arma::vec temps(d.col(0));
    const uint32 n_temps(temps.n_elem);

    const uint32 n_draws = theta.n_rows;

    double logdet_sum = 0;

    uint32 n_pars = 4;

    double val, sign;
    bool ok;
    arma::mat G(n_temps, n_pars, arma::fill::none);
    arma::mat I(n_pars, n_pars, arma::fill::none);

    double warm_gap, warm_gap_b, cold_gap;

    for (uint32 i = 0; i < n_draws; i++) {

        const double& ctmin(theta.at(i, 0));
        const double& ctmax(theta.at(i, 1));
        const double& b(theta.at(i, 2));
        const double& a(theta.at(i, 3));

        for (uint32 j = 0; j < n_temps; j++) {
            const double& t(temps.at(j));
            if ((t > ctmin) & (t < ctmax)) {
                // inside [ctmin, ctmax] range
                warm_gap = ctmax - t;
                warm_gap_b = std::pow(warm_gap, b);
                cold_gap = t - ctmin;
                // dy / dctmin:
                G.at(j,0) = -a * t * warm_gap_b;
                // dy / dctmax:
                G.at(j,1) = a * t * cold_gap * b * std::pow(warm_gap, b - 1.0);
                // dy / db:
                G.at(j,2) = a * t * cold_gap * warm_gap_b * std::log(warm_gap);
                // dy / da:
                G.at(j,3) = t * cold_gap * warm_gap_b;
            } else {
                // NOT inside [ctmin, ctmax] range
                G.row(j).zeros();
            }
        }


        I = G.t() * G + 1e-10 * arma::eye(arma::size(I));

        // ----------------------------------------

        ok = arma::log_det(val, sign, I);
        logdet_sum += val;
    }


    double logdet_mean = logdet_sum / static_cast<double>(n_draws);

    return logdet_mean;
}
