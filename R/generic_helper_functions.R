################################################################################
#' Run function f in environment e
#'
#' @param f a function
#' @param e the environment to run f in
#'
#' @family helper functions
#'
#' @keywords internal
.with_env <- function(f, e = parent.frame()) {
  stopifnot(is.function(f))
  environment(f) <- e
  f
}


################################################################################
#' Clean ggplot2 env's by generating a ggplot2 grob then converting to a ggplot2
#'
#' Useful to remove large objects before writing to disk with qs or rds
#'
#' @family helper
#' @importFrom ggplot2 ggplotGrob
#' @importFrom ggpubr as_ggplot
#'
#' @keywords internal
.grobify_ggplot <- function(p) {
  if ("ggplot" %in% class(p)) {
    p <- ggplot2::ggplotGrob(p)
    p <- ggpubr::as_ggplot(p)
  }
  return(p)
}

################################################################################
#' Remove objects from a ggplot plot_env
#'
#' Useful to remove large objects before writing to disk with qs or rds
#' This function was deprecated by p$plot_env <- rlang::new_environment()
#'
#' @family helper
#'
#' @keywords internal
.clean_ggplot_plot_env <- function(p,
                                   drop_classes = c(
                                     "SingleCellExperiment",
                                     "ggplot"
                                   )) {
  if ("ggplot" %in% class(p)) {
    env_classes <- lapply(p$plot_env, class)
    drop_idx <- purrr::map_lgl(env_classes, ~ any(drop_classes %in% .))
    drop_names <- names(env_classes[drop_idx])
    for (drop_name in drop_names) {
      p$plot_env[[eval(drop_name)]] <- NULL
    }
    p$plot_env$... <- NULL
  }
  return(p)
}



################################################################################
#' Remove unwanted variables from a quosure environment
#'
#' Useful to remove large objects before writing to disk with qs or rds
#'
#' @family helper
#'
#' @keywords internal
.clean_quosure_env <- function(quo,
                               drop_vars = c(
                                 "sce"
                               )) {
  if ("quosure" %in% class(quo)) {
    quo_env <- rlang::quo_get_env(quo)
    drop_idx <- purrr::map_lgl(names(quo_env), ~ any(drop_vars %in% .))
    drop_names <- names(quo_env)[drop_idx]
    for (drop_name in drop_names) {
      quo_env[[eval(drop_name)]] <- NULL
    }
    quo_env$... <- NULL
    quo <- rlang::quo_set_env(quo, quo_env)
  }
  return(quo)
}

################################################################################
#' Remove unwanted env variables from ggplot quosures
#'
#' Useful to remove large objects before writing to disk with qs or rds
#'
#' @family helper
#'
#' @keywords internal
.clean_ggplot_quosures <- function(p, ...) {
  fargs <- list(...)
  if ("ggplot" %in% class(p)) {
    # clean mapping quosures
    p$mapping <- lapply(p$mapping, .clean_quosure_env)
    # clean layer closures
    for (layer in seq_along(p$layers)) {
      p$layers[[layer]]$mapping <- lapply(
        p$layers[[layer]]$mapping, function(x) {
          do.call(.clean_quosure_env, c(list(quo = x), fargs))
        }
      )
    }
  }
  return(p)
}


################################################################################
#' Retrieve a palette of n_colours discrete colours from paletteer
#'
#' Useful to remove large objects before writing to disk with qs or rds
#'
#' @family helper
#'
#' @importFrom paletteer paletteer_d
#' @importFrom grDevices colorRampPalette
#'
#' @keywords internal
.get_d_palette <- function(pal_name = "ggsci::nrc_npg", n_colours = NULL) {
  palette_choice <- paletteer::paletteer_d("ggsci::nrc_npg")
  get_palette <- grDevices::colorRampPalette(palette_choice)
  if (n_colours <= length(palette_choice)) {
    pal_values <- palette_choice[1:n_colours]
  } else {
    pal_values <- get_palette(n_colours)
  }
  return(pal_values)
}


################################################################################
#' scale_y_continuous(labels=.fancy_scientific)
#'
#' @family helper
#'
#' @keywords internal
.fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # prevent 0 x xx
  l <- gsub("0e\\+00", "0", l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (3x10^+2 -> 3x10^2)
  l <- gsub("e\\+", "e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text = l)
}

################################################################################
#' Copy of archived uniftest::kolmogorov.unif.test
#'
#' @family helper
#'
#' @keywords internal
.uniftest <- function(x, nrepl = 2000, k = 0) {
  DNAME <- deparse(substitute(x))
  l <- 0
  n <- length(x)
  if (k == 1) {
    d <- max(c(1:n) / (n + 1) - x)
    for (i in 1:nrepl) {
      z <- runif(n)
      D <- max(c(1:n) / (n + 1) - z)
      if (D > d) {
        l <- l + 1
      }
    }
  }
  if (k == 0) {
    d <- max(abs(c(1:n) / (n + 1) - x))
    for (i in 1:nrepl) {
      z <- runif(n)
      D <- max(abs(c(1:n) / (n + 1) - z))
      if (D > d) {
        l <- l + 1
      }
    }
  }
  if (k == -1) {
    d <- max(x - c(1:n) / (n + 1))
    for (i in 1:nrepl) {
      z <- runif(n)
      D <- max(z - c(1:n) / (n + 1))
      if (D > d) {
        l <- l + 1
      }
    }
  }
  p.value <- l / nrepl
  RVAL <- list(
    statistic = c(D = d), p.value = p.value, method = "Kolmogorov-Smirnov test for uniformity",
    data.name = DNAME
  )
  class(RVAL) <- "htest"
  return(RVAL)
}
