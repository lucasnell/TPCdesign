suppressPackageStartupMessages({
    library(tidyverse)
    library(patchwork)
    library(nls.multstart)
    library(TPCdesign)
    library(ggtext)
    library(scico)
})

options("mc.cores" = max(1L, parallel::detectCores()-2L))


set_theme(theme_classic() +
              theme(strip.background = element_blank()))

#' This only gets run if my local .Rprofile has been run and if it's an
#' interactive session:
if (interactive() && exists(".LAN_USER")) {
    setHook(packageEvent("grDevices", "onLoad"),
            function(...) grDevices::quartz.options(width = 4, height = 4,
                                                    pointsize = 10))
    options("device" = "quartz")
    grDevices::graphics.off()
}
