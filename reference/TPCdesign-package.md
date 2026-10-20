# TPCdesign: Design Informative Thermal Performance Experiments

Tools to choose experimental temperatures that are maximally informative
for fitting thermal performance curves (TPCs). Given rough prior guesses
at a Brière-2 TPC's parameters,
[`design_temps()`](https://lucasnell.github.io/TPCdesign/reference/design_temps.md)
uses an approximate coordinate-exchange (ACE) algorithm to search for
the set of temperatures that minimizes the expected uncertainty in
fitted TPC parameters, and
[`sim_gamma_data()`](https://lucasnell.github.io/TPCdesign/reference/sim_gamma_data.md)
simulates data from that curve for testing a design.

## Main functions

- [`design_temps()`](https://lucasnell.github.io/TPCdesign/reference/design_temps.md):
  Choose a set of temperatures to sample from.

- [`sim_gamma_data()`](https://lucasnell.github.io/TPCdesign/reference/sim_gamma_data.md):
  Simulate performance data with observation error at a set of
  temperatures.

- [`briere2_tpc()`](https://lucasnell.github.io/TPCdesign/reference/briere2_tpc.md):
  Evaluate the Brière-2 TPC.

- [`briere2_tpc_Topt()`](https://lucasnell.github.io/TPCdesign/reference/briere2_tpc_Topt.md):
  Compute the thermal optimum of a Brière-2 TPC.

## See also

Useful links:

- <https://github.com/lucasnell/TPCdesign>

- <https://lucasnell.github.io/TPCdesign/>

- Report bugs at <https://github.com/lucasnell/TPCdesign/issues>

## Author

**Maintainer**: Lucas A. Nell <lucnell@gmail.com>
([ORCID](https://orcid.org/0000-0003-3209-0517))

Authors:

- Lucas A. Nell <lucnell@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-3209-0517))

- Kelsey Lyberger <Kelsey.Lyberger@asu.edu>
  ([ORCID](https://orcid.org/0000-0003-4903-380X))
