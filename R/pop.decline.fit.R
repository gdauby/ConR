#' @title Fit Statistical Models to Population Reduction
#'
#' @description Fitting statistical models to the decline on the number of mature
#'  individuals across time "can be used to extrapolate population trends, so that 
#'  a reduction of three generations can be calculated" (IUCN 2019). This function
#'  provided a comparison of five different models and returns the model with best
#'  fit to data.
#'
#' @param pop.size a vector containing the (estimated) number of mature individuals of the species 
#' @param years a vector containing the years for which the number of mature individuals was estimated 
#' @param models a vector containing the names of the models to be fitted to species population data 
#' @param project.years a vector containing the years for which the number of mature individuals should be predicted
#' @param plot logical. Should the fit of the best model be plotted against the population data?
#'
#' 
#' @details TO BE CHECKED By default, the function compare the fit of five statistical models to the population trends, 
#'  namely: linear, accelerating, exponential, logistic and piece-wise. Here, we use a quadratic model as
#'  a proxy AIC for time intervals > 30 and AICc for intervals smaller than 30.
#'  Note that model fits are performed in the proportional population estimates and on the period of interval,
#'  and not on raw population sizes and years of assessments. So, models should only be used to project 
#'  proportions 
#'  Fitting piece-wise models is very unstable (model fitting is sensitive to the start parameters) and may take a while to converge (use preferably for many years of assessment).
#'  We assumed that populations are stable or declining, so proportions are obtained in respect to the oldest population
#'  estimate.
#'  
#'  The model to be fitted should be based on the pattern of decline (which may be exponential, linear, accelerated, or a more complex pattern),
#'  which may be inferred from the type of threat. The assumed pattern of decline can make an important difference. Assessors should indicate the 
#'  basis on which they have decided the form of the decline function. The best information about the processes that contribute to changes in population 
#'  size should be used to decide what form of decline function to apply over the three-generation period. Specifically, if a model is fitted, the assumptions 
#'  of the model must be justified by characteristics of life history, habitat biology, pattern of exploitation or other threatening processes, etc. For example:
#' 
#' @author Lima, R.A.F. & Dauby, G.
#'
#' @references IUCN 2019. Guidelines for Using the IUCN Red List Categories and Criteria. Version 14. Standards and Petitions Committee. Downloadable from: http://www.iucnredlist.org/documents/RedListGuidelines.pdf.
#'
#' @export pop.decline.fit
#'
#' @examples
#' pop = c(10000, 9100, 8200, 7500, 7000)
#' yrs = c(1970, 1975, 1980, 1985, 1990)
#' pop.decline.fit(pop.size = pop, years = yrs, models = c("linear"))
#' pop.decline.fit(pop.size = pop, years = yrs, models = c("linear", "exponential"))
#' pop.decline.fit(pop.size = pop, years = yrs, models = c("linear", "exponential"), project.years = c(1960, 2005))
#' pop.decline.fit(pop.size = pop, years = yrs, models = c("linear", "exponential"), project.years = c(1960, 2005), plot = FALSE)
#' 
pop.decline.fit <- function(pop.size, 
                            years,
                            models = "all", 
                            project.years = NULL,
                            plot = TRUE) {
  
  if(is.null(years)) { 
    anos = as.numeric(gsub("[^0-9]", "", names(x)[grepl("[0-9]", names(x))]))
    if(is.null(anos)) { 
      stop("please provide at least two years with estimates of population sizes") 
    }
  } else {
    anos = years
  }
  
  if(length(as.numeric(pop.size)) > length(as.numeric(anos)))
    pop.size = pop.size[grepl(paste0(years, collapse = "|"), names(pop.size))]

  years = anos
  
  if(length(as.numeric(pop.size)) < 3 | length(years) < 3) {
    stop("too few time intervals to fit trends of population reduction")
  }  

  if(!is.null(project.years)) {
    proj = project.years[!project.years %in% years]
    years1 = c(years, proj)
    pop.size1 = c(as.numeric(pop.size), rep(NA, length(proj)))
    DATA = cbind.data.frame(pop.size = pop.size1[order(years1)], 
                            years = years1[order(years1)])
  } else {
    DATA = cbind.data.frame(pop.size = as.numeric(pop.size), 
                            years = as.numeric(years))
  }
  DATA$ps = DATA$pop.size/max(pop.size, na.rm = TRUE)
  DATA$ys = DATA$years - min(DATA$years, na.rm = TRUE)
  
  if(all(models == "all")) { 
    model.ls = vector("list", 5) 
    names(model.ls) = c("linear", "accelerating", "exponential", "logistic", "general_logistic", "piecewise") 
  } else { 
    model.ls = vector("list", length(models)) 
    names(model.ls) = models 
  }
  
  if(any("linear" %in% models) | all(models == "all")) { # linear model (Figure 4.2 panel c)
    f <- ps ~ a + b*ys
    sts = as.numeric(coef(lm(ps ~ ys, data = na.omit(DATA))))
    lin <- try(nls(f, data = na.omit(DATA),
               start = list(a = sts[1], b = sts[2])), TRUE)
    if (class(lin) == "try-error") 
      lin = lm(ps ~ ys, data = na.omit(DATA))
    model.ls[["linear"]] = lin
  }
  
  if(any("accelerating" %in% models) | all(models == "all")) { # quadratic (accelerating) model (Figure 4.2 panel d)
    f <- ps ~ a + b*ys + c*(ys^2)
    sts = as.numeric(coef(lm(ps ~ ys + I(ys^2), data = na.omit(DATA))))
    quad <- try(nls(f, data = na.omit(DATA),
                start = list(a = sts[1], b = sts[2], c = sts[3]),
                nls.control(maxiter = 500)), TRUE)
    if (class(quad) == "try-error") 
      quad <- nls_multstart(f, data = na.omit(DATA),
                           start_lower = c(a=0.1, b=-0.5, c= -0.1),
                           start_upper = c(a=1, b=0.5, c= 0.1),
                           iter = 500, supp_errors = 'Y', convergence_count = 100,
                           na.action = na.omit)
    if (class(quad) == "try-error") 
      quad = lm(ps ~ ys + I(ys^2), data = na.omit(DATA))
    model.ls[["accelerating"]] = quad
  }
  
  if(any("exponential" %in% models) | all(models == "all")) { # (negative) exponential model (Figure 4.2 panel b)
    f <- ps ~ a * exp(b * ys)
    exp <- try(nls(f, data = na.omit(DATA), 
                   start = list(a= 1, b= -0.1)), TRUE) 
    if (class(exp) == "try-error") { 
      exp <- nls_multstart(f, data = na.omit(DATA),
                           start_lower = c(a=0.1, b=-0.1),
                           start_upper = c(a=1, b=0.01),
                           iter = 500, supp_errors = 'Y', convergence_count = 100,
                           na.action = na.omit)
    }
    model.ls[["exponential"]] = exp
  }
  
  if(any("logistic" %in% models) | all(models == "all")) { # logistic model (not in IUCN guidelines)
    f <- ps ~ exp(a + b*ys)/(1 + exp(a + b*ys))
    logis <- try(nls(f, data = na.omit(DATA), 
                     start=list(a=1,b=-0.1)), TRUE) # logistic
    if (class(logis) == "try-error") { 
      logis <- nls_multstart(f, data = na.omit(DATA),
                             start_lower = c(a=1, b=-0.5),
                             start_upper = c(a=5, b=0.1),
                             iter = 500, supp_errors = 'Y', convergence_count = 100,
                             na.action = na.omit)
    }
    model.ls[["logistic"]] = logis
  }
  
  if(any("general_logistic" %in% models) | all(models == "all")) { # logistic model (not in IUCN guidelines)
    f <- ps ~ a + (1 - a)/(1 + exp(b*(ys - m)))
    # a: lower asymptote; k = upper asymptote (fixed to 1); b: growth rate; m = location

    gen.logis <- try(nls(f, data = na.omit(DATA), 
                     start=list(a=min(na.omit(DATA)$ps), b=0.2, m = median(na.omit(DATA)$ys))), TRUE) # logistic
    if (class(logis) == "try-error") { 
      gen.logis <- nls_multstart(f, data = na.omit(DATA),
                             start_lower = c(a=0, b=0.5, m=max(na.omit(DATA)$ys)),
                             start_upper = c(a=1, b=0.01,m=0),
                             iter = 500, supp_errors = 'Y', convergence_count = 100,
                             na.action = na.omit)
    }
    model.ls[["general_logistic"]] = gen.logis
  }
  
  if(any("piecewise" %in% models) | all(models == "all")) { # piece-wise model (Figure 4.2 panel a)
    ys = na.omit(DATA)$ys
    breaks <- ys[which(ys >= min(ys)+1 & ys <= max(ys)-1)]
    mse <- numeric(length(breaks))
    for(j in 1:length(breaks)) {
      piecewise1 <- lm(ps ~ ys*(ys < breaks[j]) + ys*(ys>=breaks[j]),
                       data = na.omit(DATA))
      mse[j] <- as.numeric(summary(piecewise1)[6])
    }
    quebras = breaks[order(mse)]
    ## Piecewise Regression METHOD 1
    #piece0 <- try(lm(ps ~ ys*(ys < quebras[1]) + ys*(ys > quebras[1])), TRUE)   
    ## Piecewise Regression METHOD 2
    md <- lm(ps ~ ys, data = na.omit(DATA)) 
    # 1 breakpoint
    piece1 <- try(segmented::segmented(md, seg.Z = ~ys, psi = quebras[1], 
                                       control=seg.control(display=FALSE), it.max = 100, n.boot = 50), TRUE)
    # 2 breakpoints
    piece2 <- try(segmented::segmented(md, seg.Z = ~ys, psi = c(quebras[1]/2, quebras[1]), 
                                       control=seg.control(display=FALSE), it.max = 100, n.boot = 50), TRUE)
    if(class(piece2)[1] == "try-error") {
      counter <- 0
      while(class(piece2)[1] == "try-error" | counter < 50){
        try(piece2 <- segmented::segmented(md, seg.Z = ~ys, psi = c(quebras[1]/2, quebras[1]),
                                           control = seg.control(display=FALSE, it.max = 100, n.boot = 50)), TRUE)
        counter <- sum(counter, 1)
      }
    }
    # 3 breakpoints
    piece3 <- try(segmented::segmented(md, seg.Z = ~ys, psi = c(quebras[1]/3, quebras[1]/2, quebras[1]), 
                                       control=seg.control(display=FALSE), it.max = 50,  n.boot = 50), TRUE)
    if(class(piece3)[1] == "try-error") {
      counter <- 0
      while(class(piece3)[1] == "try-error" | counter < 50){
        try(piece3 <- segmented::segmented(md, seg.Z = ~ys, psi = c(quebras[1]/3, quebras[1]/2, quebras[1]),
                                           control = seg.control(display=FALSE, it.max = 50, n.boot = 50)), TRUE)
        counter <- sum(counter, 1)
      }
    }
    # Chosing the best piece-wise model
    piece.mds = list(piece1, piece2, piece3)
    piece.mds = piece.mds[!sapply(piece.mds, function(x) class(x[1])) %in% "try-error"]
    
    if(length(na.omit(DATA)$ys) >= 30) {
      AICs = AICtab(piece.mds, sort = FALSE)
    } else {
      AICs = AICctab(piece.mds, sort = FALSE)
    }
    
    if(all(AICs$dAIC <= log(8))) {
      piece = piece.mds[which(AICs$df == min(AICs$df))]
    } else {
      piece = piece.mds[which(AICs$dAIC == min(AICs$dAIC))]
    }
    model.ls[["piecewise"]] = piece[[1]]
  }
  
  if(length(years) >= 30) {
    AICs = AICtab(model.ls, sort = FALSE, mnames = names(model.ls))
  } else {
    AICs = AICctab(model.ls, sort = FALSE, mnames = names(model.ls))
  }

  if(any(AICs$dAIC <= log(8))) {
    if(all(AICs$dAIC <= log(8))) {
      best = model.ls[which(AICs$df == min(AICs$df))][[1]]
      best.name = names(model.ls[which(AICs$df == min(AICs$df))])
    } else {
      id = which(AICs$dAIC <= log(8))
      id = id[which(AICs$df[id] == min(AICs$df[id]))]
      best = model.ls[id][[1]]
      best.name = names(model.ls[id])
    }  
  } else {
    best = model.ls[which(AICs$dAIC == min(AICs$dAIC))][[1]]
    best.name = names(model.ls[which(AICs$dAIC == min(AICs$dAIC))])
  }
  
  attributes(best)$best.model.name = best.name

  DATA$est.prop = predict(best, list(ys = DATA$ys))
  DATA$predicted = round(DATA$est.prop * max(DATA$pop.size, na.rm = TRUE), 1)
  
  if(plot == TRUE) {
    ylim = range(DATA[grepl("ps|est.prop", names(DATA))], na.rm=TRUE)
    xlim = range(DATA$ys, na.rm=TRUE)
    par(mfrow= c(1, 1), mgp = c(2.8,0.6,0), mar= c(4,4,1,1))
    plot(ps ~ ys, pch=19, cex=1.2, data = DATA, 
         ylim = ylim, xlim = xlim,
         xlab = "Years", ylab = "Population size (%)",
         xaxt = "n", yaxt= "n", cex.lab = 1.2)
    axis(1, at = DATA$ys, label = DATA$years, tcl = -0.3)
    ats = pretty(seq(min(ylim), max(ylim), 0.1))
    axis(2, at = ats, label = ats*100, las = 1, tcl = -0.3)
    if(!is.null(project.years)) points(est.prop ~ ys, cex=1.2, data = subset(DATA, is.na(DATA$ps)))
    
    range1 = (range(na.omit(DATA)$ys)[1] - 1)
    range2 = (range(na.omit(DATA)$ys)[2] + 1)
    if(best.name == "linear") { 
      mod = model.ls[["linear"]]
      if(!is.null(project.years)) curve(coef(mod)[1] + coef(mod)[2]*x, 
              add= TRUE, lwd= 2, lty= 2, col = "#D55E00")
      curve(coef(mod)[1] + coef(mod)[2]*x, from = range1, to = range2,
            add= TRUE, lwd= 2, col = "#D55E00")
      legend("bottomleft", c("Linear model"), lwd= 2, col= "#D55E00", bty = "n")
    }  
    if(best.name == "accelerating") { 
      mod = model.ls[["accelerating"]]
      if(!is.null(project.years)) curve(coef(mod)[1] + coef(mod)[2]*x + coef(mod)[3]*(x^2),
              add= TRUE, lwd=2, lty= 2, col = "#56B4E9")
      curve(coef(mod)[1] + coef(mod)[2]*x + coef(mod)[3]*(x^2), from = range1, to = range2,
            add= TRUE, lwd=2, col = "#56B4E9")
      legend("bottomleft", c("Accelerating (quadratic) model"), lwd= 2, col= "#56B4E9", bty = "n")
    }
    if(best.name == "exponential") { 
      mod = model.ls[["exponential"]]
      if(!is.null(project.years)) curve(coef(mod)[1]*exp(coef(mod)[2]*x), 
              add= TRUE, lwd=2, lty= 2, col = "#009E73")
      curve(coef(mod)[1]*exp(coef(mod)[2]*x), from = range1, to = range2,
            add= TRUE, lwd=2, col = "#009E73")
      legend("bottomleft", c("Exponential model"), lwd= 2, col= "#009E73", bty = "n")
    }    
    if(best.name == "logistic") { 
      mod = model.ls[["logistic"]]
      if(!is.null(project.years)) curve(exp(coef(mod)[1] + coef(mod)[2]*x)/(1 + exp(coef(mod)[1] + coef(mod)[2]*x)),
              lwd=2, lty= 2, col = "#F0E442", add= TRUE)
      curve(exp(coef(mod)[1] + coef(mod)[2]*x)/(1 + exp(coef(mod)[1] + coef(mod)[2]*x)), from = range1, to = range2,
            lwd=2, col = "#F0E442", add= TRUE)
      legend("bottomleft", c("Exponential model"), lwd= 2, col= "#F0E442", bty = "n")
    }
    if(best.name == "general_logistic") { 
      mod = model.ls[["general_logistic"]]
      if(!is.null(project.years)) curve(coef(mod)[1] + (1 - coef(mod)[1])/(1 + exp(coef(mod)[2] * (x - coef(mod)[3]))),
                                        lwd=2, lty= 2, col = "#0072B2", add= TRUE)
      curve(coef(mod)[1] + (1 - coef(mod)[1])/(1 + exp(coef(mod)[2] * (x - coef(mod)[3]))), from = range1, to = range2,
            lwd=2, col = "#0072B2", add= TRUE)
      legend("bottomleft", c("Generalized logistic model"), lwd= 2, col= "#0072B2", bty = "n")
    }
    if(best.name == "piecewise") { 
      mod = model.ls[["piecewise"]]
      plot(mod, lwd=2, col= "#CC79A7", add=TRUE, dens.rug=FALSE, rug=FALSE)
      abline(v=mod$psi[,2], lty=3, col= "#CC79A7")
      legend("bottomleft", c("Piece-wise model", "Est. breaks"), lwd= 2, lty = c(1,3), col= "#CC79A7", bty = "n")
      #curve((coef(tmp4)[1] + coef(tmp4)[3]) + (coef(tmp4)[2]+coef(tmp4)[5])*x, add=T, from=min(ys), to=quebras[1], lwd=2, col=4)
      #curve((coef(tmp4)[1] + coef(tmp4)[4]) + coef(tmp4)[2]*x, add=T, from=quebras[1], to=max(ys), lwd=2, col=4)
    }
  }

  if(!is.null(project.years)) {
    DATA1 = DATA[,c("pop.size", "years", "ps", "est.prop","predicted")] 
    names(DATA1) = c("Observed", "Year", "Obs_proportion", "Est_prportion", "Predicted")
  } else {
    DATA1 = DATA[,c("pop.size", "years", "ps")] 
    names(DATA1) = c("Observed", "Year", "Obs_proportion")
  }
  
  if(length(models) > 1 | all(models == "all")) {
    res = list(best.model = best, model.selection.result = AICs, data = DATA1)
  } else { 
    res = list(best.model = best, data = DATA1)
  }
  
  return(res)
}  