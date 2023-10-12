#' @title Quantification of Population Fluctuation
#'
#' @description Based on a time series of population sizes, this function try to 
#' tell appart extreme fluctuations from directional population changes. 
#'
#' @param x a vector (one species) or a data frame (multiple species/
#'   subpopulations) containing the population sizes (e.g. number of mature
#'   individuals) per year, from the oldest to the most recent estimate.
#' @param years a vector containing the years for which the population sizes is
#'   available (i.e. time series).
#' @param plot.test logical. Should the the results be plotted with the popuation data?
#'
#' @return A data frame containing the mean change of the population
#'   ("Magnitude.fluctuation"), the percentage of intervals presenting
#'   population increases followed by decreases ("Alternance.prop"), the result
#'   of a test for the presence of directional changes in the population size
#'   ("Direct.change"), and the Standard Error Estimate of the linear regression
#'   model fitted to the population trend ("Std.Error.Est.").
#'             
#' @details This is a basic function that quantifies the mean fluctuation of
#'   population sizes across time, aiming at the detection 'extreme
#'   fluctuations' as defined by IUCN (2019). Here we quantify flucutuations as
#'   the mean change in population size between consecutive years in the time
#'   series (e.g. if t[i]= 9 and t[i+1]= 90, change is 10). As defined in IUCN
#'   (2019), extreme fluctuations are generally characterizes by changes higher
#'   than 10 or an order of magnitude.
#'   
#'   (Detailed xxplanation of 'alternance.prop' pending)
#'   
#'   The evidence of directional changes is evaluated based on a linear
#'   regression model fitted to the population size data. The sign of the slope
#'   parameter of the regression is used to assess if population declining or
#'   increasing. The confidence interval around the slope estimate is used to
#'   evaluated if the declining or increasing tendencies are significant
#'   (significance level of 0.05). 
#'   
#'   The same linear regression model is used to obtain the standard error of
#'   the estimate (SEE) of the linear regression fitted to the population trend,
#'   which can be used as a measure of the temporal variability in population
#'   size (Cuervo & Møller 2017).
#' 
#' @author Renato A. Ferreira de Lima
#'
#' @references 
#'  IUCN 2019. Guidelines for Using the IUCN Red List Categories and Criteria.
#'  Version 14. Standards and Petitions Committee. Downloadable from:
#'  http://www.iucnredlist.org/documents/RedListGuidelines.pdf. 
#'  
#'  Cuervo, J.J. & Møller, A.P. (2017). Colonial, more widely distributed and
#'  less abundant bird species undergo wider population fluctuations independent
#'  of their population trend. PloS one, 12(3): e0173220.
#'  
#' @export pop.fluctuation
#'
#' @examples
#' 
#' ## Examples adapted from Figure 4.4 in IUCN (2019)
#' data("example_fluctuation")
#' pop.fluctuation(x = example_fluctuation)
#' 
pop.fluctuation <- function(x, 
                            years= NULL, 
                            plot.test = TRUE) {
  
  if(is.null(x))
    stop("Please provide at least two estimates of population sizes")
  
  if(is.null(years)) {
    
    anos <- as.numeric(gsub("[^0-9]", "", names(x)[grepl("[0-9]", names(x))]))
    
    if(is.null(anos)) 
      stop("Please provide at least two years with estimates of population sizes") 
    
    years <- anos
    warning("The years of the population sizes were not given and were taken from the input population data", call. = FALSE)
  }
  
  if(is.vector(x)) {
    
    if(is.null(names(x))) {
      
      x = as.data.frame(matrix(x, ncol = length(x), dimnames = list(NULL, years)),
                        stringsAsFactors = FALSE)
      
    } else {
      
      x = as.data.frame(matrix(x, ncol = length(x), dimnames = list(NULL, names(x))),
                        stringsAsFactors = FALSE)
      
    }
  }
  
  if(length(years) < 2)
    stop("At least two years are needed to perform the assessment")
  
  if(class(x[,1]) %in% c("factor", "character")) {
    
    nomes <- x[,1]
    x <- x[,-1]
    
  } else { nomes = NULL }
  
  if(!all(names(x) %in% years)) {
    x <- x[,names(x) %in% years]
  }
  
  if(plot.test) {
    
    panels = c(ceiling(sqrt(dim(x)[1])), ceiling(sqrt(dim(x)[1]))) 
    if(panels[1]>3) panels[1] <- 3
    if(panels[2]>3) panels[2] <- 3
    par(mfrow = panels, mgp = c(2.5,0.5,0), mar= c(3.5,4,0.5,1), las=1, tcl = -0.25)
    ylim <- range(x, na.rm=TRUE) + c(0,5)
    
  }  

  res = NULL
  for(j in 1:dim(x)[1]) {
    
    time <- years
    obs  <- x[j,]
    time <- time[match(names(obs), time, nomatch = 0)]
    time <- time[!is.na(obs)]
    obs  <- obs[!is.na(obs)]

    # Mean magnitude of the fluctuations
    fluct <- NULL
    for(i in 1:(length(obs)-1)) {
      
      t1 <- obs[i]
      t2 <- obs[i+1]
      if(t2 > t1) fluct[i] <- (t2-t1)/t1
      if(t2 < t1) fluct[i] <- (t1-t2)/t2
      
    }
    mean.fluct = round(mean(fluct, na.rm = TRUE), 2)
    
    #Standard error of the estimate
    mod <- stats::lm(obs ~ time)
    #mean.sq.errors <- mean((predict(mod) - obs)^2)
    see <- round(stats::sigma(mod), 2)
    
    #Coeficient of variation   
    #coef.var <- round(raster::cv(obs), 2)
    
    #Alternance in the direction of the change?
    diffs <- diff(obs)
    diffs[diffs>(-1) & diffs<1] <- 0 # very small changes in number os individuals considered as zero
    diffs <- diffs[diffs != 0]
    signs <- sign(diffs)
    
    #### CHECK HERE: HOW TO BETTER COMPARE IF CHANGES ARE ALWAYS IN THE OPPOSITE DIRECTION (ALTERNANCE) ####
    perf.signs <- rep(c(signs[1],-1*signs[1]), ceiling(length(signs)/2))
    perf.signs <- perf.signs[1:length(signs)]
    alter.prop <- round(100*sum(signs == perf.signs)/length(signs),2)
    
    #Any evidence of directional population change?
    ci <- stats::confint(mod)[2,]
    
    if(stats::coef(mod)[2] < 0)
      test <- if(ci[1]<0 & ci[2]<0) "signif.decline" else "non.signif.decline"
    
    if(stats::coef(mod)[2] >= 0)
      test <- if(ci[1]>=0 & ci[2]>0) "signif.increase" else "non.signif.increase"
    
    result <- as.character(c(mean.fluct, alter.prop, test, see))
    if(is.null(res)) res = result else res = rbind.data.frame(res, result, stringsAsFactors = FALSE)

    if(plot.test) {
      
      plot(obs ~ time,
           xlab = "Years", ylim = ylim, cex.lab = 1.2, pch = 19,cex = 0.9,
           type = "b", ylab = "Population size", log = "y")
      graphics::curve(stats::coef(mod)[1] + stats::coef(mod)[2]*x, add=TRUE, lwd=2,col=2)
      
      #leg.pos <- auto.legend.pos(obs, time, xlim = range(time), ylim = ylim)
      legend("bottomleft", paste(c("Mean fluct.= "),mean.fluct), bty = "n")
      if(!is.null(nomes)) legend("topright", nomes[j], bty = "n")
      
    }
    
  }

 if(!is.null(nomes)) row.names(res) = nomes  
 names(res) = c("Magnitude.fluctuation", "Alternance.prop", "Direct.change", "Std.Error.Est.")  
 return(res)
}