# ifndef __TPCDESIGN_SIM_DATA_H
# define __TPCDESIGN_SIM_DATA_H


#include "TPCdesign_types.h"

#include <algorithm>
#include <math.h>

#include "tpc.h"
#include "pcg.h"

using namespace Rcpp;





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





#endif
