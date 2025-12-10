
# install.packages(c("ggplot2", "hexSticker", "cowplot")) # if needed
library(ggplot2)
library(hexSticker)
library(cowplot)

set.seed(1)

## Shared data ----
y <- seq(0, 10, length.out = 600)
x_left  <- sin(2 * pi * y / 2.0) * 0.6
x_right <- sin(2 * pi * y / 2.0 + pi) * 0.6

helix_left  <- data.frame(x = x_left,  y = y)
helix_right <- data.frame(x = x_right, y = y)

y_rungs <- seq(0.3, 9.7, by = 0.35)
rungs <- data.frame(
  x1 = sin(2 * pi * y_rungs / 2.0) * 0.6,
  x2 = sin(2 * pi * y_rungs / 2.0 + pi) * 0.6,
  y  = y_rungs
)

mk_plot <- function(DNAcol, segmentcol) {
  cpg_y <- y_rungs[sample(c(TRUE, FALSE), length(y_rungs), replace = TRUE)]
  cpg <- data.frame(
    x = sin(2 * pi * cpg_y / 2.0) * 0.6 + 0.08,
    y = cpg_y
  )

  ggplot() +
    # two strands (thinner)
    geom_path(data = helix_left,  aes(x, y),
              color = DNAcol, linewidth = 0.35, lineend = "round") +
    geom_path(data = helix_right, aes(x, y),
              color = DNAcol, linewidth = 0.35, lineend = "round") +
    # rungs (even thinner)
    geom_segment(data = rungs, aes(x = x1, xend = x2, y = y, yend = y),
                 color = segmentcol, linewidth = 0.25, lineend = "round") +
    # methylation dots (smaller)
    geom_point(data = cpg, aes(x, y), fill = "red", color = "black", pch =21, size = .8, alpha = 1) +
    theme_void() +
    coord_fixed(xlim = c(-1.1, 1.1), ylim = c(-0.5, 10.5), expand = FALSE)
}

p1 <- mk_plot("#0072B2", "#56B4E9")
p2 <- mk_plot("#0d7554", "#10c434")
p3 <- mk_plot("#fcba03", "#fce549")

p <- cowplot::plot_grid(p1, p2, p3, nrow = 1)

dir.create("inst/figures", recursive = TRUE, showWarnings = FALSE)

sticker(
  p,
  package = "hyperVarMeth",
  p_size  = 20,
  p_color = "black",
  s_x     = 1, s_y = 0.8,
  s_width = 1.0, s_height = 1.0,
  h_fill  = "white",
  h_color = "#261f1f",
  dpi     = 300,                  # higher dpi keeps thin lines crisp
  filename = "inst/figures/logo.png",
  url = "https://github.com/alicebalard/hyperVarMeth", u_size = 3
)
