# ifndef __TPCDESIGN_UTIL_H
# define __TPCDESIGN_UTIL_H


#include "TPCdesign_types.h"

#include <algorithm>
#include <math.h>

#include "pcg.h"

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






inline arma::vec briere2_tpc_cpp(const arma::vec& temp,
                                 const double& ctmin,
                                 const double& ctmax,
                                 const double& a,
                                 const double& b,
                                 const bool& scale) {

    arma::vec out(temp.n_elem, arma::fill::none);
    double out_max = 0;
    for (uint32 i = 0; i < temp.n_elem; i++) {
        out.at(i) = a * temp.at(i) * std::max(temp.at(i) - ctmin, 0.0) *
            std::pow(std::max(ctmax - temp.at(i), 0.0), b);
        if (out.at(i) > out_max) out_max = out.at(i);
    }
    if (scale) {
        if (out_max <= 0) out_max = 1;
        for (double& x : out) x /= out_max;
    }
    return out;
}







/*
 Normal distribution truncated above zero.
 Used for generating performance measures bc we never want them to be < 0.
 */
class trunc_normal_distribution {

    double mu;
    double sigma;
    double a_bar;
    double p;
    /*
     This method does weird things numerically if mu < 0 and its magnitude is
     greater than 5 * sigma. This is because this situation results in
     a distribution that essentially never produces anything > 0.
     The bool below checks for this situation and will make this always
     produce 0 instead of trying to generate numbers.
     */
    bool always_zero; // make output always zero if

public:

    inline trunc_normal_distribution()
        : mu(0),
          sigma(1),
          a_bar((0 - mu) / sigma),
          p(1),
          always_zero(false) {}
    inline trunc_normal_distribution(const double& mu_, const double& sigma_)
        : mu(mu_), sigma(sigma_), a_bar((0 - mu) / sigma),
          p(R::pnorm5(a_bar, 0, 1, 1, 0)),
          always_zero(mu < (-5 * sigma)) {}
    inline trunc_normal_distribution(const trunc_normal_distribution& other)
        : mu(other.mu), sigma(other.sigma), a_bar((0 - mu) / sigma),
          p(R::pnorm5(a_bar, 0, 1, 1, 0)),
          always_zero(mu < (-5 * sigma)) {}

    inline trunc_normal_distribution& operator=(const trunc_normal_distribution& other) {
        mu = other.mu;
        sigma = other.sigma;
        a_bar = (0 - mu) / sigma;
        p = R::pnorm5(a_bar, 0, 1, 1, 0);
        always_zero = mu < (-5 * sigma);
        return *this;
    }

    inline void reset(const double& mu_, const double& sigma_) {
        mu = mu_;
        sigma = sigma_;
        a_bar = (0 - mu) / sigma;
        p = R::pnorm5(a_bar, 0, 1, 1, 0);
        always_zero = mu < (-5 * sigma);
        return;
    }

    inline double operator()(pcg32& eng) {

        if (always_zero) return 0;
        if (sigma == 0) return mu;

        double u = runif_ab(eng, p, 1);

        double x = R::qnorm5(u, 0, 1, 1, 0);
        x = x * sigma + mu;

        return x;
    }


};






/*
 =======================================================================
 =======================================================================
 Manual log determinant
 =======================================================================
 =======================================================================
 */

/*
 I created this function so that the result of the log determinant
 didn't change when the BLAS library did.
 I included the tolerance based singularity check to avoid non-zero determinant
 outcomes as a result of rounding.
 In R on my system, this issue would arise with
 `as.numeric(determinant(matrix(1:9/10, 3), logarithm = TRUE)$modulus)`,
 which would result in `-40.06104`.
 The matrix's third column is a linear combination of the first two
 (col3 = 2*col2 - col1), so the determinant should be exactly 0 and the
 log-determinant should  be `-Inf`.
 This function returns the correct answer.
 */

// ------------------------------------------------------------
// Manually compute log|det(A)| via LU decomposition.
// ------------------------------------------------------------
inline double log_det_cpp(const arma::mat& A) {

    if (A.n_rows != A.n_cols) {
        stop("Matrix must be square.");
    }

    arma::mat L, U, P;
    arma::lu(L, U, P, A);   // A = P^T * L * U  (Armadillo's LU with partial pivoting)

    int n = A.n_rows;

    // ---- Tolerance-based singularity check ----
    // Instead of comparing u_ii to exactly 0, compare it to the
    // scale of the matrix. A common, robust choice mirrors what
    // LAPACK/Armadillo use internally for rank/condition checks:
    //   tol = n * machine_epsilon * max(|U_ii|)
    arma::vec abs_diag = arma::abs(U.diag());
    double max_diag = abs_diag.is_empty() ? 0.0 : abs_diag.max();
    double tol = n * std::numeric_limits<double>::epsilon() * max_diag;


    // log|det(A)| = sum(log|U_ii|)   since det(L) = 1 (unit diagonal)
    // det(A) = det(P^T) * det(L) * det(U) = det(P^T) * det(U)
    double log_det = 0.0;
    // int sign = 1;

    for (int i = 0; i < n; ++i) {
        double u_ii = U(i, i);

        if (std::abs(u_ii) <= tol) {
            // Effectively singular -> determinant is 0
            return -arma::datum::inf;
        }

        log_det += std::log(std::abs(u_ii));
        // if (u_ii < 0) sign *= -1;
    }

    // // Account for the sign of the permutation matrix P.
    // // det(P) = (-1)^(number of row swaps). We recover this from
    // // how far P deviates from the identity via its determinant.
    // double det_P = arma::det(P);   // will be exactly +1 or -1
    // if (det_P < 0) sign *= -1;

    return log_det;
}


#endif
