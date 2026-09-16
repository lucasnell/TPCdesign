#' TPCdesign: Design Informative Thermal Performance Experiments
#'
#' Tools to choose experimental temperatures that are maximally informative
#' for fitting thermal performance curves (TPCs). Given rough prior guesses
#' at a Brière-2 TPC's parameters, [design_temps()] uses an approximate
#' coordinate-exchange (ACE) algorithm to search for the set of temperatures
#' that minimizes the expected uncertainty in fitted TPC parameters, and
#' [sim_gamma_data()] simulates data from that curve for testing a design.
#'
#' @section Main functions:
#' \itemize{
#'   \item [design_temps()]: Choose a set of temperatures to sample from.
#'   \item [sim_gamma_data()]: Simulate performance data with observation
#'     error at a set of temperatures.
#'   \item [briere2_tpc()]: Evaluate the Brière-2 TPC.
#'   \item [briere2_tpc_Topt()]: Compute the thermal optimum of a Brière-2 TPC.
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib TPCdesign, .registration = TRUE
## usethis namespace: end
NULL
