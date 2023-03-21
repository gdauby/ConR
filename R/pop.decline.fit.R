#' @title Fit Statistical Models to Population Reduction
#'
#' @description Fitting statistical models to the decline on the number of
#'   mature individuals across time "can be used to extrapolate population
#'   trends, so that a reduction of three generations can be calculated" (IUCN
#'   2019). This function provide a comparison of five different models and
#'   returns the predictions of the model with best fit to data.
#'
#' @param pop.size a vector containing the (estimated) number of mature
#'   individuals of the species
#' @param years a vector containing the years for which the population sizes is
#'   available
#' @param models a vector containing the names of the statistical models to be
#'   fitted to species population data
#' @param project.years a vector containing the years for which the number of
#'   mature individuals should be predicted
#' @param plot.fit logical. Should the fit of the best model be plotted against
#'   the population data?
#' @param max.count numerical. Maximum number of attempts to fit piece-wise
#'   models. Default to 50.
#' @param ... other parameters to be passed as arguments for function ```ICtab.mod.select```
#'
#' @details
#' By default, the function compares the fit of six statistical models
#' to the population trends, namely: linear, quadratic, exponential, logistic,
#' generalized logistic and piece-wise. But, as stated in IUCN (2019), the
#' model used to do the predictions makes an important difference. So, model
#' fit to data should not be the only or most important criteria to choose
#' among models. Users should preferably choose one or two of the models based
#' on the best available information of types of threat (i.e. patterns of
#' exploitation or habitat loss), life history and ecology of the taxon being
#' evaluated or any other processes that may contribute to population decline.
#' See IUCN (2019) for more details on the assumptions of each model.
#' 
#' The linear and exponential patterns of decline are fully described in IUCN
#' (2019) and are easy to be described statistically through a model (see
#' Figure 4.2, pg. 33 of IUCN 2019). But IUCN (2019) also recognizes the
#' existence of more "complex patterns of decline". To describe more complex
#' patterns, ```pop.decline.fit``` provides fits to logistic and piece-wise
#' patterns of decline. Despite the options of models provided by
#' ```pop.decline.fit```, depending on the numbers of observations or the patterns
#' of decline, many or none of the models may provide a good fit to data. This
#' reinforces the role of the user in choosing the more appropriate pattern
#' for the area or taxon considered.
#' 
#' For simplicity, the population size data provided is transformed into
#' proportions using the maximum population estimate provided. Therefore,
#' models are fit to proportional data, but the projections are provided in
#' proportions and in the original scale. As suggested in IUCN (2019), no
#' model fit is performed if only two estimates of population size are
#' provided.
#' 
#' Some more technical notes on model fitting and selection. Here, we use a
#' quadratic model as an equivalent to the accelerating model described in 
#' IUCN (2019), but note that the quadratic model can generate non-realistic 
#' projections depending on the population data or on the years chosen for the 
#' projection (see example). Fitting piece-wise models can be unstable (model 
#' fitting is quite sensitive to the start parameters) and may take a while to 
#' converge; so, it should preferably be used when  many years of population 
#' data are available. For simplicity, only piece-wise models with up to 3 
#' breaks and linear functions between breaks are provided. For time intervals > 80, 
#' the best model among the candidate models is chosen based on Akaike 
#' Information Criterion, or AIC; the corrected AIC or the AICc (Burnham and 
#' Anderson, 2004) is used for time intervals < 80.
#' 
#' 
#' @author Lima, R.A.F. & Dauby, G.
#'
#' @references IUCN 2019. Guidelines for Using the IUCN Red List Categories and
#'   Criteria. Version 14. Standards and Petitions Committee. Downloadable from:
#'   http://www.iucnredlist.org/documents/RedListGuidelines.pdf.
#'
#' @examples
#' ## Creating vectors with the population data and time intervals 
#' #(adapted from the IUCN 2019 workbook for Criterion A, available 
#' #at: https://www.iucnredlist.org/resources/criterion-a)
#' pop = c(10000, 9100, 8200, 7500, 7000)
#' yrs = c(1970, 1975, 1980, 1985, 1990)
#' 
#' ## Fitting data with different models and setting
#' pop.decline.fit(pop.size = pop, years = yrs)
#' pop.decline.fit(pop.size = pop, years = yrs, models = c("linear"))
#' pop.decline.fit(pop.size = pop, years = yrs, models = c("linear", "exponential","quadratic"))
#' pop.decline.fit(pop.size = pop, years = yrs, models = c("linear", "exponential"), project.years = c(1960, 2005))
#' pop.decline.fit(pop.size = pop, years = yrs, models = c("linear", "exponential"), project.years = c(1973, 2005))
#' pop.decline.fit(pop.size = pop, years = yrs, models = c("quadratic", "general_logistic"))
#' pop.decline.fit(pop.size = pop, years = yrs, models = c("general_logistic"))
#' pop.decline.fit(pop.size = pop, years = yrs, models = c("quadratic"), project.years = c(1960, 2050))
#' 
#' ## Another examples, with less observations (warning or no model fit below 3 observations)
#' pop.decline.fit(pop.size = c(10000, 8200, 6000), years = c(1970, 1985, 2000), models = "all", project.years = 2030)
#' pop.decline.fit(pop.size = c(10000, 6000), years = c(1970, 2000))
#' 
#' @importFrom nls.multstart nls_multstart
#' @importFrom segmented segmented seg.control
#' @importFrom stats na.omit
#'
#' @export pop.decline.fit
#' 
pop.decline.fit <- function(pop.size, 
                            years,
                            models = "all", 
                            project.years = NULL,
                            plot.fit = TRUE,
                            max.count = 50,
                            ...) {
  
  if(is.null(years)) {
    
    anos <- as.numeric(gsub("[^0-9]", "", names(pop.size)[grepl("[0-9]", names(pop.size))]))
    
    if(is.null(anos)) {
      
      stop("Please provide at least two years with estimates of population sizes")
      
    }
  } else {
    
    anos <- years
    
  }
  
  if(length(as.numeric(pop.size)) > length(as.numeric(anos)))
    pop.size <- pop.size[grepl(paste0(years, collapse = "|"), names(pop.size))]
  
  years <- anos
  
  if(length(as.numeric(pop.size)) < 3 | length(years) < 3)
    stop("Too few time intervals to fit trends of population reduction")
  
  if(!is.null(project.years)) {
    
    proj <- project.years[!project.years %in% years]
    years1 <- c(years, proj)
    pop.size1 <- c(as.numeric(pop.size), rep(NA, length(proj)))
    DATA <- cbind.data.frame(pop.size = as.numeric(pop.size1[order(years1)]), 
                             years = as.numeric(years1[order(years1)]),
                             stringsAsFactors = FALSE)
    
  } else {
    
    DATA <- cbind.data.frame(pop.size = as.numeric(pop.size), 
                             years = as.numeric(years),
                             stringsAsFactors = FALSE)
    
  }
  
  DATA$ps <- DATA$pop.size/max(as.double(pop.size), na.rm = TRUE)
  DATA$ys <- DATA$years - min(DATA$years, na.rm = TRUE)
  
  if(all(models == "all")) {
    
    nomes.mds <- c(
      "linear",
      "quadratic",
      "exponential",
      "logistic",
      "general_logistic",
      "piecewise")
    model.ls <- vector("list", length(nomes.mds))
    names(model.ls) <- nomes.mds
    
  } else {
    
    model.ls <- vector("list", length(models))
    names(model.ls) <- models
    
  }
  
  if(any("linear" %in% models) | all(models == "all")) { # linear model (Figure 4.2 panel c)
    
    f <- ps ~ a + b*ys
    sts <- as.numeric(stats::coef(stats::lm(ps ~ ys, data = stats::na.omit(DATA))))
    lin <- try(stats::nls(f, data = stats::na.omit(DATA),
                          start = list(a = sts[1], b = sts[2])), TRUE)
    
    if (class(lin) == "try-error") 
      lin = stats::lm(ps ~ ys, data = stats::na.omit(DATA))
    
    model.ls[["linear"]] = lin
    
  }
  
  if(any("quadratic" %in% models) | all(models == "all")) { # quadratic (accelerating) model (Figure 4.2 panel d)
    
    f <- ps ~ a + b*ys + c*(ys^2)
    sts <- as.numeric(stats::coef(stats::lm(ps ~ ys + I(ys^2), data = stats::na.omit(DATA))))
    quad <- try(stats::nls(f, data = stats::na.omit(DATA),
                           start = list(a = sts[1], b = sts[2], c = sts[3]),
                           stats::nls.control(maxiter = 500)), TRUE)
    
    if (class(quad) == "try-error")
      quad <- suppressWarnings( nls.multstart::nls_multstart(f, data = stats::na.omit(DATA),
                                                             start_lower = c(a=0.1, b=-0.5, c= -0.1),
                                                             start_upper = c(a=1, b=0.5, c= 0.1),
                                                             iter = 500, supp_errors = 'Y', convergence_count = 100,
                                                             na.action = stats::na.omit)
      )
    
    if (class(quad) == "try-error") 
      quad = stats::lm(ps ~ ys + I(ys^2), data = stats::na.omit(DATA))
    
    model.ls[["quadratic"]] = quad
    
  }
  
  if(any("exponential" %in% models) | all(models == "all")) { # (negative) exponential model (Figure 4.2 panel b)
    
    f <- ps ~ a * exp(b * ys)
    exp <- try(stats::nls(f, data = stats::na.omit(DATA), 
                          start = list(a= 1, b= -0.1)), TRUE) 
    
    if (class(exp) == "try-error") { 
      
      exp <- suppressWarnings( nls.multstart::nls_multstart(f, data = stats::na.omit(DATA),
                                                            start_lower = c(a=0.1, b=-0.1),
                                                            start_upper = c(a=1, b=0.01),
                                                            iter = 500, supp_errors = 'Y', convergence_count = 100,
                                                            na.action = stats::na.omit)
      )
      
    }
    
    model.ls[["exponential"]] = exp
    
  }
  
  if(any("logistic" %in% models) | all(models == "all")) { # logistic model (not formally in IUCN guidelines)
    
    f <- ps ~ exp(a + b*ys)/(1 + exp(a + b*ys))
    logis <- try(stats::nls(f, data = stats::na.omit(DATA), 
                            start = list(a=1,b=-0.1)), TRUE)
    
    if (class(logis) == "try-error") { 
      
      logis <- suppressWarnings( nls.multstart::nls_multstart(f, data = stats::na.omit(DATA),
                                                              start_lower = c(a=1, b=-0.5),
                                                              start_upper = c(a=5, b=0.1),
                                                              iter = 500, supp_errors = 'Y', convergence_count = 100,
                                                              na.action = na.omit)
      )
      
    }
    
    model.ls[["logistic"]] = logis
    
  }
  
  if(any("general_logistic" %in% models) | all(models == "all")) { # gen. logistic model (not formally in IUCN guidelines)
    
    f <- ps ~ a + (1 - a)/(1 + exp(-b*(ys - m))) # a: lower asymptote; k = upper asymptote (fixed to 1); b: growth rate; m = location
    
    gen.logis <- try(stats::nls(f, data = stats::na.omit(DATA), 
                                start=list(a=min(stats::na.omit(DATA)$ps), b=-0.2, m = stats::median(stats::na.omit(DATA)$ys))), TRUE)
    
    if (class(gen.logis) == "try-error") { 
      
      gen.logis <- suppressWarnings( nls.multstart::nls_multstart(f, data = stats::na.omit(DATA),
                                                                  start_lower = c(a=0, b=-0.5, m=max(stats::na.omit(DATA)$ys)),
                                                                  start_upper = c(a=1, b=-0.01,m=0),
                                                                  iter = 500, supp_errors = 'Y', convergence_count = 100,
                                                                  na.action = na.omit)
      )
      
    }
    
    model.ls[["general_logistic"]] = gen.logis
    
  }
  
  if(any("piecewise" %in% models) | all(models == "all")) { # piece-wise model (Figure 4.2 panel a)
    
    ys <- stats::na.omit(DATA)$ys
    breaks <- ys[which(ys >= min(ys)+1 & ys <= max(ys)-1)]
    mse <- numeric(length(breaks))
    
    for(j in seq_along(breaks)) {
      
      piecewise1 <- stats::lm(ps ~ ys*(ys < breaks[j]) + ys*(ys >= breaks[j]),
                              data = stats::na.omit(DATA))
      mse[j] <- as.numeric(summary(piecewise1)[6])
      
    }
    
    quebras <- breaks[order(mse)]
    md <- stats::lm(ps ~ ys, data = stats::na.omit(DATA)) 
    
    # 1 breakpoint
    piece1 <- try(segmented::segmented(md, seg.Z = ~ys, psi = quebras[1], 
                                       control = segmented::seg.control(display = FALSE), it.max = 100, n.boot = 50), TRUE)
    # 2 breakpoints
    piece2 <- try(segmented::segmented(md, seg.Z = ~ys, psi = c(quebras[1]/2, quebras[1]), 
                                       control = segmented::seg.control(display = FALSE), it.max = 100, n.boot = 50), TRUE)
    warn <- warnings(piece2)
    
    warn.patt <- "no residual degrees of freedom"
    if(class(piece2)[1] == "try-error" & !any(grepl(warn.patt, attributes(warn)$names, ignore.case = TRUE))) {
      
      counter <- 0
      
      while(class(piece2)[1] == "try-error" & counter < max.count){
        
        try(piece2 <- segmented::segmented.default(md, seg.Z = ~ys, psi = c(jitter(quebras[1]/2, 1), jitter(quebras[1], 1)),
                                           control = segmented::seg.control(display = FALSE, it.max = 100, n.boot = 50)), TRUE)
        counter <- sum(counter, 1)
        
      }
    }
    
    # 3 breakpoints
    piece3 <- try(segmented::segmented(md, seg.Z = ~ys, psi = c(jitter(quebras[1]/3, 1), jitter(quebras[1]/2, 1), jitter(quebras[1], 1)), 
                                       control = segmented::seg.control(display = FALSE), it.max = 100,  n.boot = 50), TRUE)
    warn <- warnings(piece3)
    
    if(class(piece3)[1] == "try-error" & !any(grepl(warn.patt, attributes(warn)$names, ignore.case = TRUE))) {
      
      counter <- 0
      
      while(class(piece3)[1] == "try-error" & counter < max.count){
        
        try(piece3 <- segmented::segmented(md, seg.Z = ~ys, psi = c(quebras[1]/3, quebras[1]/2, quebras[1]),
                                           control = segmented::seg.control(display = FALSE, it.max = 100, n.boot = 50)), TRUE)
        counter <- sum(counter, 1)
        
      }
    }
    
    # Chosing the best piece-wise model (without errors and with break-point estimates)
    piece.mds <- list(piece1, piece2, piece3)
    piece.mds <- piece.mds[!sapply(piece.mds, function(x) class(x)[1] %in% "try-error")]
    piece.mds <- piece.mds[!sapply(piece.mds, function(x) is.null(x$psi))]
    
    
    if(length(piece.mds) > 0) {
      
      if(length(piece.mds) == 1) {
        
        piece = piece.mds[[1]]
        
      } else {
        
        md.sel = ICtab.mod.select(piece.mds, parsimony = TRUE, ...)
        piece <- md.sel$best.model
        
      }
      
      model.ls[["piecewise"]] <- piece
      
    } else {
      
      model.ls[["piecewise"]] <- NA
      warning("No piece-wise model was fit to population data due to lack of convergence, probably caused by too few observations")
      
    }    
  }
  
  md.sel <- ICtab.mod.select(model.ls, ...)
  AICs <- md.sel$ICtab
  best <- md.sel$best.model
  best.name <- md.sel$best.model.name
  attributes(best)$best.model.name <- best.name
  
  DATA$est.prop <- predict(best, data.frame(ys = DATA$ys))
  DATA$est.prop[DATA$est.prop<0] <- 0
  DATA$predicted <- round(DATA$est.prop * max(DATA$pop.size, na.rm = TRUE), 1)
  
  if(plot.fit) {
    
    ylim <- range(DATA[grepl("ps|est.prop", names(DATA))], na.rm=TRUE)
    xlim <- range(DATA$ys, na.rm=TRUE)
    preds <- predict(best, data.frame(ys = seq(range(DATA$ys)[1], range(DATA$ys)[2], by = 1)))
    
    if(!is.null(project.years) & (any(min(preds) < min(ylim)) | any(max(preds) < max(ylim)))) 
      ylim <- range(c(ylim, preds), na.rm=TRUE)     
    
    par(mfrow= c(1, 1), mgp = c(2.8,0.6,0), mar= c(4,4,1,1))
    graphics::plot(DATA$ps ~ DATA$ys, pch=19, cex=1.2, #data = DATA, 
         ylim = ylim, xlim = xlim,
         xlab = "Years", ylab = "Population size (%)",
         xaxt = "n", yaxt= "n", cex.lab = 1.2)
    axis(1, at = DATA$ys, labels = DATA$years, tcl = -0.3)
    ats <- pretty(seq(min(ylim), max(ylim), 0.1))
    axis(2, at = ats, labels = ats*100, las = 1, tcl = -0.3)
    
    if(!is.null(project.years)) 
      points(est.prop ~ ys, cex=1.2, data = subset(DATA, is.na(DATA$ps)))
    
    range1 <- (range(stats::na.omit(DATA)$ys)[1] - 1)
    range2 <- (range(stats::na.omit(DATA)$ys)[2] + 1)
    
    leg.pos <- auto.legend.pos(DATA$ps, DATA$ys, xlim = xlim, ylim = ylim)
    
    if(best.name == "linear") { 
      
      #mod <- model.ls[["linear"]]
      mod <- best
      if(!is.null(project.years)) graphics::curve(stats::coef(mod)[1] + stats::coef(mod)[2]*x, 
                                        add= TRUE, lwd= 2, lty= 2, col = "#D55E00")
      graphics::curve(stats::coef(mod)[1] + stats::coef(mod)[2]*x, from = range1, to = range2,
            add= TRUE, lwd= 2, col = "#D55E00")
      legend(leg.pos, c("Linear model"), lwd= 2, col= "#D55E00", bty = "n")
      
    }
    
    if(best.name == "quadratic") { 
      
      #mod <- model.ls[["quadratic"]]
      mod <- best
      if(!is.null(project.years)) graphics::curve(stats::coef(mod)[1] + stats::coef(mod)[2]*x + stats::coef(mod)[3]*(x^2),
                                        add= TRUE, lwd=2, lty= 2, col = "#56B4E9")
      graphics::curve(stats::coef(mod)[1] + stats::coef(mod)[2]*x + stats::coef(mod)[3]*(x^2), from = range1, to = range2,
            add= TRUE, lwd=2, col = "#56B4E9")
      legend(leg.pos, c("Quadratic model"), lwd= 2, col= "#56B4E9", bty = "n")
      
    }
    
    if(best.name == "exponential") { 
      
      #mod <- model.ls[["exponential"]]
      mod <- best
      if(!is.null(project.years)) graphics::curve(stats::coef(mod)[1]*exp(stats::coef(mod)[2]*x), 
                                        add= TRUE, lwd=2, lty= 2, col = "#009E73")
      graphics::curve(stats::coef(mod)[1]*exp(stats::coef(mod)[2]*x), from = range1, to = range2,
            add= TRUE, lwd=2, col = "#009E73")
      legend(leg.pos, c("Exponential model"), lwd= 2, col= "#009E73", bty = "n")
      
    }
    
    if(best.name == "logistic") { 
      
      #mod <- model.ls1[["logistic"]]
      mod <- best
      if(!is.null(project.years)) graphics::curve(exp(stats::coef(mod)[1] + stats::coef(mod)[2]*x)/(1 + exp(stats::coef(mod)[1] + stats::coef(mod)[2]*x)),
                                        lwd=2, lty= 2, col = "#F0E442", add= TRUE)
      graphics::curve(exp(stats::coef(mod)[1] + stats::coef(mod)[2]*x)/(1 + exp(stats::coef(mod)[1] + stats::coef(mod)[2]*x)), from = range1, to = range2,
            lwd=2, col = "#F0E442", add= TRUE)
      legend(leg.pos, c("Logistic model"), lwd= 2, col= "#F0E442", bty = "n")
      
    }
    
    if(best.name == "general_logistic") { 
      
      #mod <- model.ls1[["general_logistic"]]
      mod <- best
      if(!is.null(project.years)) graphics::curve(stats::coef(mod)[1] + (1 - stats::coef(mod)[1])/(1 + exp(-stats::coef(mod)[2] * (x - stats::coef(mod)[3]))),
                                        lwd=2, lty= 2, col = "#0072B2", add= TRUE)
      graphics::curve(stats::coef(mod)[1] + (1 - stats::coef(mod)[1])/(1 + exp(-stats::coef(mod)[2] * (x - stats::coef(mod)[3]))), from = range1, to = range2,
            lwd=2, col = "#0072B2", add= TRUE)
      legend(leg.pos, c("Gen. logistic model"), lwd= 2, col= "#0072B2", bty = "n")
      
    }
    
    if(best.name == "piecewise") { 
      
      #mod <- model.ls1[["piecewise"]]
      mod <- best
      plot(mod, lwd=2, col= "#CC79A7", add=TRUE, dens.rug=FALSE, rug=FALSE)
      graphics::abline(v=mod$psi[,2], lty=3, col= "#CC79A7")
      legend(leg.pos, c("Piecewise model", "Estim. break(s)"), lwd= 2, lty = c(1,3), col= "#CC79A7", bty = "n")

    }
  }
  
  if(!is.null(project.years)) {
    
    DATA1 <- DATA[,c("pop.size", "years", "ps", "est.prop","predicted")] 
    names(DATA1) <- c("Observed", "Year", "Obs_proportion", "Est_proportion", "Predicted")
    
  } else {
    
    DATA1 <- DATA[,c("pop.size", "years", "ps")] 
    names(DATA1) <- c("Observed", "Year", "Obs_proportion")
    
  }
  
  if(length(models) > 1 | all(models == "all")) {
    
    res <- list(best.model = best, model.selection.result = AICs, data = DATA1)
    
  } else { 
    
    res <- list(best.model = best, data = DATA1)
    
  }
  
  return(res)
} 
