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
