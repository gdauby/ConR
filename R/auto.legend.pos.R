#' @title Set Legend Position Automatically
#'
#' @param x a vector containing the (estimated) number of mature individuals of
#'   the species
#' @param y a vector containing the years for which the number of mature
#'   individuals was estimated
#' @param xlim a vector containing the names of the models to be fitted to
#'   species population data
#' @param ylim a vector containing the years for which the number of mature
#'   individuals should be predicted
#'
#' @author user 'chan1142' in stackoverflow
#'
#' @references https://stackoverflow.com/questions/7198178/automatically-determine-position-of-plot-legend
#' 
#' @noRd
#' 
#' 
auto.legend.pos <- function(x, y, xlim=range(x), ylim=range(y)) {
  countIt <- function(a, zero.only = TRUE) {
    tl <- sum(x <= xlim[1]*(1-a)+xlim[2]*a & y >= ylim[1]*a+ylim[2]*(1-a), na.rm = TRUE)
    tr <- sum(x >= xlim[1]*a+xlim[2]*(1-a) & y >= ylim[1]*a+ylim[2]*(1-a), na.rm = TRUE)
    bl <- sum(x <= xlim[1]*(1-a)+xlim[2]*a & y <= ylim[1]*(1-a)+ylim[2]*a, na.rm = TRUE)
    br <- sum(x >= xlim[1]*a+xlim[2]*(1-a) & y <= ylim[1]*(1-a)+ylim[2]*a, na.rm = TRUE)
    c(topleft=tl, topright=tr, bottomleft=bl, bottomright=br)
  }
  for (k in seq(0.5,0.05,by=-0.05)) {
    a <- countIt(k)
    if (sum(a==0)>0) break
  }
  names(a)[which(a==0)][1]
}