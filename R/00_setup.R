# 00_setup.R -- packages, colours and a plain theme shared by every figure
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

col_placebo <- "#2a78d6"   # blue
col_hu      <- "#eb6834"   # orange
col_third   <- "#1baf7a"   # aqua, used when a third series is needed
col_ink     <- "#333333"
col_muted   <- "#8a8a86"

theme_pmx <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#e9e8e4", linewidth = 0.3),
    axis.title       = element_text(colour = col_ink),
    axis.text        = element_text(colour = "#52514e"),
    strip.text       = element_text(face = "bold", hjust = 0, colour = col_ink),
    plot.title       = element_text(face = "bold", size = 12, colour = col_ink),
    plot.subtitle    = element_text(size = 10, colour = "#52514e"),
    legend.position  = "top",
    legend.title     = element_blank()
  )
theme_set(theme_pmx)

save_fig <- function(plot, file, width = 8, height = 4.5) {
  dir.create("figures", showWarnings = FALSE)
  ggsave(file.path("figures", file), plot,
         width = width, height = height, dpi = 200, bg = "white")
}
