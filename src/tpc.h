# ifndef __TPCDESIGN_TPC_H
# define __TPCDESIGN_TPC_H


#include "TPCdesign_types.h"

#include <algorithm>
#include <math.h>


using namespace Rcpp;





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






#endif
