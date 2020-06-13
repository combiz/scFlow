################################################################################
#' Run function f in environment e
#'
#' @param f a function
#' @param e the environment to run f in
#'
#' @family helper functions
#'
#' @keywords internal
.with_env <- function(f, e=parent.frame()) {
  stopifnot(is.function(f))
  environment(f) <- e
  f
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
    # clean plot_env
    #p$plot_env <- rlang::new_environment()
    # clean mapping quosures
    p$mapping <- lapply(p$mapping, .clean_quosure_env)
    # clean layer closures
    for (layer in seq_along(p$layers)) {
      p$layers[[layer]]$mapping <- lapply(
        p$layers[[layer]]$mapping, function(x)
        do.call(.clean_quosure_env, c(list(quo = x), fargs))
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
  l <- gsub("0e\\+00","0",l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (3x10^+2 -> 3x10^2)
  l <- gsub("e\\+","e",l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}
