
library(acebayes)
library(lhs)


# utility_briere2D <- function(d, B) {
#     temps  <- d[, 1]
#     theta  <- B
#     ndraws <- nrow(theta)
#     logdets <- numeric(ndraws)
#
#     for (i in 1:ndraws) {
#         ctmin <- theta[i, "ctmin"]
#         ctmax <- theta[i, "ctmax"]
#         b <- theta[i, "b"]
#         a <- theta[i, "a"]
#
#         inside <- (temps > ctmin) & (temps < ctmax)
#
#         warm_gap <- pmax(ctmax - temps, 0)
#         cold_gap <- pmax(temps - ctmin, 0)
#
#         dy_dctmin <- -a * temps * warm_gap^b
#         dy_dctmax <- a * temps * cold_gap * b * warm_gap^(b - 1)
#         dy_db     <- a * temps * cold_gap * warm_gap^b * log(warm_gap)
#
#         G <- cbind(dy_dctmin, dy_dctmax, dy_db)
#         G[!inside, ] <- 0
#
#         I <- t(G) %*% G
#         logdets[i] <- as.numeric(determinant(I + diag(1e-10, nrow(I)), logarithm = TRUE)$modulus)
#     }
#
#     mean(logdets)
# }

Rcpp::cppFunction(
'double utility_briere2D_cpp(const arma::mat& d, SEXP B) {

    // Because pace passes a list sometimes:
    arma::mat theta;
    if (TYPEOF(B) == VECSXP) {
        Rcpp::List lB(B);
        Rcpp::NumericMatrix mB(lB[0]);
        theta = arma::mat(mB.begin(), mB.nrow(), mB.ncol(), false);
    } else if (TYPEOF(B) == REALSXP) {
        Rcpp::NumericMatrix mB(B);
        theta = arma::mat(mB.begin(), mB.nrow(), mB.ncol(), false);
        // const arma::mat theta(B);
    } else stop("B must be list or numeric");

    typedef arma::uword uint32;

    const arma::vec temps(d.col(0));
    const uint32 ntemps(temps.n_elem);

    const uint32 ndraws = theta.n_rows;

    double logdet_sum = 0;

    double val, sign;
    bool ok;
    arma::mat G(ntemps, 3U, arma::fill::none);
    arma::mat I(3U, 3U, arma::fill::none);

    double warm_gap, warm_gap_b, cold_gap;

    for (uint32 i = 0; i < ndraws; i++) {

            const double& ctmin(theta.at(i, 0));
            const double& ctmax(theta.at(i, 1));
            const double& b(theta.at(i, 2));
            const double& a(theta.at(i, 3));

            for (uint32 j = 0; j < ntemps; j++) {
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
                } else {
                    // NOT inside [ctmin, ctmax] range
                    G.at(j,0) = 0.0;
                    G.at(j,1) = 0.0;
                    G.at(j,2) = 0.0;
                }
            }


            I = G.t() * G + 1e-10 * arma::eye(arma::size(I));

            // ----------------------------------------

            ok = arma::log_det(val, sign, I);
            logdet_sum += val;
        }


    double logdet_mean = logdet_sum / static_cast<double>(ndraws);

    return logdet_mean;
}
',
depends = "RcppArmadillo")



# Fills in gaps to hedge against curve mismatch
gap_filler <- function(n_fill, opt_temps, min_temp, max_temp) {

    if (n_fill <= 0) return(opt_temps)

    points <- sort(opt_temps)

    # Include the domain boundaries as candidate gap edges
    points <- c(min_temp, points, max_temp)

    for (i in 1:n_fill) {
        gaps <- diff(points)
        largest_gap_idx <- which.max(gaps)
        new_point <- mean(points[largest_gap_idx:(largest_gap_idx + 1)])
        points <- sort(c(points, new_point))
    }

    # Drop the boundary markers if they weren't part of the original optimal set
    final <- setdiff(points, c(min_temp, max_temp))
    # unless lower/upper were themselves chosen as optimal points, add them back
    if (min_temp %in% opt_temps) final <- c(min_temp, final)
    if (max_temp %in% opt_temps) final <- c(max_temp, upper)

    final <- sort(round(final, 2))
    return(final)
}





ace_design_temps <- function(n_temps, ctmin, ctmax, a, b,
                             opt_temp_p = 0.7,
                             n_draws = 500L, n_starts = 10L,
                             n_threads = 1L) {

    # n_temps = 7L; ctmin = 5; ctmax = 45; a = 1; b = 0.2; n_draws = 500; n_starts = 10L

    n_optimal <- as.integer(round(opt_temp_p * n_temps))
    n_fill <- n_temps - n_optimal  # insurance against shape mismatch

    # genus-level relatives was 2.50°C for Tmin, 1.50°C for Tmax,
    # and 0.26 for log(b)
    prior_ctmin <- ctmin + c(-1,1) * 2.5
    prior_ctmax <- ctmax + c(-1,1) * 1.5
    prior_lb <- log(b) + c(-1,1) * 0.26


    theta_draws <- cbind(ctmin = runif(n_draws, prior_ctmin[1], prior_ctmin[2]),
                         ctmax = runif(n_draws, prior_ctmax[1], prior_ctmax[2]),
                         b = exp(runif(n_draws, prior_lb[1], prior_lb[2])),
                         a = rep(1, n_draws))


    min_temp <- min(prior_ctmin)
    max_temp <- max(prior_ctmax)
    start_temp_list <- lapply(1:n_starts, function(i) {
        d <- randomLHS(n = n_optimal, k = 1) * (max_temp - min_temp) + min_temp
        colnames(d) <- "temp"
        return(d)
    })

    design_robust <- pace(utility = utility_briere2D_cpp,
                          start.d = start_temp_list,
                          B = rep(list(theta = theta_draws), 2),
                          lower = min_temp - 5,
                          upper = max_temp + 5,
                          deterministic = TRUE,
                          mc.cores = n_threads)

    opt_temps <- sort(round(design_robust$d, 2))

    final_temps <- gap_filler(n_fill, opt_temps, min_temp, max_temp)

    return(final_temps)

}

