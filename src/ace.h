# ifndef __TPCDESIGN_ACE_H
# define __TPCDESIGN_ACE_H


#include "TPCdesign_types.h"

#include <algorithm>
#include <math.h>


using namespace Rcpp;





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
