
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

## bug Fix
theme_sticker <- function (size = 1.2, ...) {
  center <- 1
  radius <- 1
  h <- radius
  w <- sqrt(3)/2 * radius
  m <- 1.02
  list(theme_transparent() + theme(plot.margin = margin(b = 0,
                                                        l = 0, unit = "lines"), strip.text = element_blank(),
                                   line = element_blank(), text = element_blank(), title = element_blank(),
                                   ...), coord_fixed(), scale_y_continuous(expand = c(0,
                                                                                      0), limits = c(center - h * m, center + h * m)), scale_x_continuous(expand = c(0,
                                                                                                                                                                     0), limits = c(center - w * m, center + w * m)))
}

sticker <- function (subplot, s_x = 0.8, s_y = 0.75, s_width = 0.4, s_height = 0.5,
          package, p_x = 1, p_y = 1.4, p_color = "#FFFFFF", p_family = "Aller_Rg",
          p_fontface = "plain", p_size = 8, h_size = 1.2, h_fill = "#1881C2",
          h_color = "#87B13F", spotlight = FALSE, l_x = 1, l_y = 0.5,
          l_width = 3, l_height = 3, l_alpha = 0.4, url = "", u_x = 1,
          u_y = 0.08, u_color = "black", u_family = "Aller_Rg", u_size = 1.5,
          u_angle = 30, white_around_sticker = FALSE, ..., filename = paste0(package,
                                                                             ".png"), asp = 1, dpi = 300)
{
  hex <- ggplot() + geom_hexagon(size = h_size, fill = h_fill,
                                 color = NA)
  if (inherits(subplot, "character")) {
    d <- data.frame(x = s_x, y = s_y, image = subplot)
    sticker <- hex + geom_image(aes(x = !!sym("x"), y = !!sym("y"),
                                    image = !!sym("image")), d, size = s_width, asp = asp)
  }
  else {
    sticker <- hex + geom_subview(subview = subplot, x = s_x,
                                  y = s_y, width = s_width, height = s_height)
  }
  sticker <- sticker + geom_hexagon(size = h_size, fill = NA,
                                    color = h_color)
  if (spotlight)
    sticker <- sticker + geom_subview(subview = spotlight(l_alpha),
                                      x = l_x, y = l_y, width = l_width, height = l_height)
  sticker <- sticker + geom_pkgname(package, p_x, p_y, color = p_color,
                                    family = p_family, fontface = p_fontface, size = p_size,
                                    ...)
  sticker <- sticker + geom_url(url, x = u_x, y = u_y, color = u_color,
                                family = u_family, size = u_size, angle = u_angle)
  if (white_around_sticker)
    sticker <- sticker + white_around_hex(size = h_size)
  sticker <- sticker + theme_sticker(size = h_size)
  save_sticker(filename, sticker, dpi = dpi)
  class(sticker) <- c("sticker", class(sticker))
  invisible(sticker)
}

sticker(
  p,
  package = "hyperVarMeth",
  p_size  = 15, spotlight = F,
  p_color = "black",
  s_x     = 1, s_y = 0.8,
  s_width = 1.0, s_height = 1.0,
  h_fill  = "white",
  h_color = "black",
  dpi     = 300,                  # higher dpi keeps thin lines crisp
  filename = "inst/figures/logo.png",
  url = "github.com/alicebalard/hyperVarMeth", u_size = 4, h_size = .5
)
