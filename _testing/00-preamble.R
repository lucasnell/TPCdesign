library(tidyverse)
library(nls.multstart)
library(TPCdesign)

options("mc.cores" = max(1L, parallel::detectCores()-2L))


set_theme(theme_classic() +
              theme(strip.background = element_blank()))
