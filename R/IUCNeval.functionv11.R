

.Convex.Hull.Poly <- function(XY) {
  hpts <- chull(x =  XY[,1], y = XY[,2]) ; hpts <- c(hpts, hpts[1])
  coord <- matrix(NA, length(hpts), 2)
  POLY <- "POLYGON(("
  for (i in 1:length(hpts)){
    POLY <- paste(POLY,XY[hpts[i],1]," ",XY[hpts[i],2], sep="")
    if(i!=length(hpts)) POLY <- paste(POLY,", ", sep="")
    if(i==length(hpts)) POLY <- paste(POLY,"))", sep="")
    
    coord[i,1] <- XY[hpts[i],2]
    coord[i,2] <- XY[hpts[i],1]
    
  }
  p1 = readWKT(POLY)
  crs(p1) <- "+proj=longlat +datum=WGS84"
  makePoly(p1)
  return(p1)
}

## The functions ahull_to_SPLDF and alpha.hull.poly were originally posted in the website https://casoilresource.lawr.ucdavis.edu/software/r-advanced-statistical-package/working-spatial-data/converting-alpha-shapes-sp-objects/
## in a now broken link. It is also used in github functions written by David Bucklin, see https://github.com/dnbucklin/r_movement_homerange 

.ahull_to_SPLDF <- function(x, proj4string=NA)
{
  if(class(x) != 'ahull')
    stop('This function only works with `ahull` class objects')
  
  # convert ashape edges to DF
  x.ah.df <- as.data.frame(x$arcs)
  
  # convert each arc to a line segment
  l.list <- list()
  for(i in 1:nrow(x.ah.df))
  {
    # extract row i
    row_i <- x.ah.df[i,]
    
    # extract elements for arc()
    v <- c(row_i$v.x, row_i$v.y)
    theta <- row_i$theta
    r <- row_i$r
    cc <- c(row_i$c1, row_i$c2)
    # from arc()
    angles <- alphahull::anglesArc(v, theta)
    seqang <- seq(angles[1], angles[2], length = 100)
    x <- cc[1] + r * cos(seqang)
    y <- cc[2] + r * sin(seqang)
    
    # convert to line segment
    l.list[[i]] <- Line(cbind(x,y))
  }
  
  # promote to Lines class, then to SpatialLines class
  l <- Lines(l.list, ID=1)
  
  # copy over CRS data from original point data
  l.spl <- SpatialLines(list(l), proj4string=CRS(as.character(NA)))
  
  # promote to SpatialLinesDataFrame, required for export to GRASS / OGR
  l.spldf <- SpatialLinesDataFrame(l.spl, data=data.frame(id=1), match.ID=FALSE)
  
  return(l.spldf)
}

.alpha.hull.poly <- function(XY, alpha=1, buff=0.1){
  
   Used_data=unique(XY)
   Used_data <- apply(Used_data,2,jitter)
   if (any(rownames(installed.packages())=="alphahull")) {
     
     loadNamespace("alphahull")
     ahull.obj <- alphahull::ahull(Used_data[,c(1,2)], alpha = alpha)
     y.as.spldf <- .ahull_to_SPLDF(ahull.obj)
     y.as.spldf_buff <- gBuffer(y.as.spldf, width=buff)
     
     NZp <- slot(y.as.spldf_buff, "polygons")
     holes <- lapply(NZp, function(x) sapply(slot(x, "Polygons"), slot, 
                                             "hole"))
     res <- lapply(1:length(NZp), function(i) slot(NZp[[i]], 
                                                   "Polygons")[!holes[[i]]])
     IDs <- row.names(y.as.spldf_buff)
     NZfill <- SpatialPolygons(lapply(1:length(res), function(i) 
       Polygons(res[[i]], ID=IDs[i])), proj4string=CRS(proj4string(y.as.spldf_buff)))
     return(NZfill)
     crs(NZfill) <- "+proj=longlat +datum=WGS84"
   }else{
     stop("The package alphahull is required for this procedure, please install it")
   }
}


.crop.poly <- function(poly, crop){

  p1_owin <- spatstat::as.owin(poly)
  africa_owin <- spatstat::as.owin(crop)
  
  if(round(area.owin(union.owin(p1_owin,africa_owin)),3)!=round(area.owin(africa_owin),3)) {
    w <- spatstat::setminus.owin(p1_owin, africa_owin)
    w2 <- spatstat::setminus.owin(p1_owin, w)
    poly_masked <- as(w2, "SpatialPolygons")
    
    raster::crs(poly_masked) <- raster::crs(crop)
  }else{
    poly_masked <- poly
  }
  EOO <- round(areaPolygon(poly_masked)/1000000,1)  
  return(list(EOO, poly_masked))
}


.EOO.comp <-  function(XY, 
                       exclude.area=FALSE, 
                       buff_width=0.1, 
                       country_map=NULL, 
                       Name_Sp="tax", 
                       alpha.hull=FALSE, 
                       convex.hull=TRUE,
                       alpha=1, buff.alpha=0.1, method.less.than3="not comp") { # , verbose=TRUE
  
  ### Checking if the method of calculating EOO has been chosen
  if(!convex.hull & !alpha.hull) 
    stop("alpha.hull and convex.hull are both FALSE, choose one of them")
  
  if(nrow(unique(XY))>1) if(max(dist(XY[,2]), na.rm=T)>=180) 
    stop(paste("EOO for species", as.character(Name_Sp),"cannot be computed because occurrences spans more than 180 degrees longitude"))
  
  ## Check if there are less than 3 unique occurrences
   if(nrow(unique(XY))<3) {
     
     ## if there is only one occurrence, EOO is NA
     if(nrow(unique(XY))<2) {
       EOO <- NA
       message(paste("EOO for", as.character(Name_Sp),"is not computed because there is only 1 unique occurrence"))
      
     }else{
       if(method.less.than3=="arbitrary") {
         
         ## if there are two unique occurences, EOO is 1/10 of the distance between the two points
         projEAC=crs("+proj=cea +lon_0=Central Meridian+lat_ts=Standard Parallel+x_0=False Easting+y_0=False Northing +ellps=WGS84")
         coordEAC <- as.data.frame(matrix(unlist(rgdal::project(as.matrix(unique(XY)[,1:2]),proj=as.character(projEAC),inv=FALSE)), ncol=2))
         EOO <- as.numeric(dist(coordEAC/1000)*0.1*dist(coordEAC/1000))  #
       }
       
       if(method.less.than3=="not comp") {
         
         ## if there are two unique occurences, EOO is not computed neither
         message(paste("EOO for", as.character(Name_Sp),"is not computed because there is less than 3 unique occurrences"))
         EOO <- NA
       }
     }
     
     OUTPUT <- round(EOO, 0)
     names(OUTPUT) <- c("EOO")
     
  }else{
    
    ### Checking if all occurrences are on a straight line
    if(length(unique(XY[,1]))==1 || length(unique(XY[,2]))==1 || round(abs(cor(XY[,1],XY[,2])),6)==1) {
      ## If so, a straight line is built and a buffer of buff_width is added
      message(paste("Occurrences of",as.character(Name_Sp),"follow a straight line, thus EOO is based on an artificial polygon using buff_width"))
      hpts <- unique(XY[,c(2,1)])
      POLY <- "LINESTRING("
      for (Z in 1:dim(hpts)[1]){
        POLY <- paste(POLY,hpts[Z,1]," ", hpts[Z,2], sep="")
        if(Z!=dim(hpts)[1]) POLY <- paste(POLY,", ", sep="")
        if(Z==dim(hpts)[1]) POLY <- paste(POLY,")", sep="")
      }
      p1 = readWKT(POLY)
      crs(p1) <- "+proj=longlat +datum=WGS84"
      
      makeLine(p1) ### Add vertices to line
      
      p1 <- gBuffer(p1, width=buff_width) ### Add buffer to line
      
      ## If exclude.area is TRUE
      if(exclude.area) {
        croped.EOO <- .crop.poly(poly=p1, crop=country_map)
        p1 <- croped.EOO[[2]]
        EOO <- croped.EOO[[1]]
      }else{
        EOO <- round(areaPolygon(p1)/1000000,0)
      }
      
    }else{
      
      if(alpha.hull) p1 <- .alpha.hull.poly(cbind(XY[,2],XY[,1]), alpha=alpha, buff=buff.alpha)
      
      if(convex.hull) p1 <- .Convex.Hull.Poly(cbind(XY[,2],XY[,1]))
      
      if(exclude.area) {
        croped.EOO <- .crop.poly(poly=p1, crop=country_map)
        p1 <- croped.EOO[[2]]
      }
      
      ## If exclude.area is TRUE
      if(exclude.area) {EOO <- croped.EOO[[1]]}else{EOO <- round(areaPolygon(p1)/1000000,0)}
    }
    
    OUTPUT <- list(EOO, p1)
    names(OUTPUT) <- c("EOO", "spatial.polygon")
  }
  
  # if(verbose) cat(" ",paste(Name_Sp,"EOO comp."))
  
  return(OUTPUT)
}


EOO.computing <- function(XY, 
                          exclude.area=FALSE, 
                          country_map=NULL, 
                          export_shp=FALSE, 
                          write_shp=FALSE, 
                          alpha=1, 
                          buff.alpha=0.1, 
                          method.range="convex.hull", 
                          Name_Sp="species1", 
                          buff_width=0.1, 
                          method.less.than3="not comp",
                          write_results=TRUE, 
                          file.name="EOO.results", 
                          parallel=FALSE, 
                          NbeCores=2, 
                          show_progress=FALSE){ # , verbose=TRUE
  
  if(any(is.na(XY[,c(1:2)]))) {
    print(paste("Skipping",length(which(rowMeans(is.na(XY[,1:2]))>0)) ,"occurrences because of missing coordinates for",  # if(verbose) 
                paste(as.character(unique(XY[which(rowMeans(is.na(XY[,1:2]))>0),3])), collapse=" AND ") ))
    XY <- XY[which(!is.na(XY[,1])),]
    XY <- XY[which(!is.na(XY[,2])),]
  }
  
  XY <- as.data.frame(XY)
  
  if(exclude.area & is.null(country_map)) stop("exclude.area is TRUE but no country_map is provided")
  if(buff_width>80) stop("buff_width has unrealistic value")
  if(any(XY[,2]>180) || any(XY[,2]< -180)|| any(XY[,1]< -180) || any(XY[,1]>180)) stop("coordinates are outside of expected range")
  
  if(method.range=="convex.hull") {
    convex.hull=TRUE
    alpha.hull=FALSE
  }
  
  if(method.range=="alpha.hull") {
    convex.hull=FALSE
    alpha.hull=TRUE
  }
  
  if(ncol(XY)>2) {
    colnames(XY)[1:3] <- c("ddlat","ddlon","tax")
    XY$tax <- as.character(XY$tax)
    list_data <- split(XY, f = XY$tax)
  }else{
    colnames(XY)[1:2] <- c("ddlat","ddlon")
    list_data <- list(XY)
  }
  
  # OUTPUT <- lapply(list_data, function(x) .EOO.comp(x, Name_Sp=ifelse(ncol(XY)>2, as.character(unique(x$tax)), Name_Sp),
  #                                                   exclude.area=exclude.area, buff_width=buff_width, country_map=country_map,
  #                                                   alpha=alpha, buff.alpha=buff.alpha, alpha.hull=alpha.hull, convex.hull=convex.hull,
  #                                                   method.less.than3=method.less.than3)) # , verbose=verbose
  
  # if(show_progress) prog. <- "text"
  # if(!show_progress) prog. <- "none"
  # 
  # if(parallel) registerDoParallel(NbeCores)
  # 
  # OUTPUT <- plyr::llply(list_data, .fun=function(x) {
  #   .EOO.comp(x, Name_Sp=ifelse(ncol(XY)>2, as.character(unique(x$tax)), Name_Sp),
  #             exclude.area=exclude.area, buff_width=buff_width, country_map=country_map,
  #             alpha=alpha, buff.alpha=buff.alpha, alpha.hull=alpha.hull, convex.hull=convex.hull,
  #             method.less.than3=method.less.than3) # , verbose=verbose
  # }
  # , .progress = prog., .parallel=parallel, .paropts = "ConR")
  
  if(parallel) {
    if("doParallel" %in% 
       rownames(installed.packages()) == FALSE) {stop("Please install doParallel package")}
    
    library(doParallel)
    
    doParallel::registerDoParallel(NbeCores)
    message('doParallel running with ',
            NbeCores, ' cores - progressbar not valid')
    `%d%` <- foreach::`%dopar%`
  }else{
    `%d%` <- foreach::`%do%`
  }
  
  if(show_progress)  pb <- 
    utils::txtProgressBar(min = 0, max = length(list_data), style = 3)
  
  if(is.null(names(list_data))) names_ <- rep(Name_Sp, length(list_data)) else names_ <- names(list_data)
  
  # ,
  # .packages = c("rgeos", "geosphere", "raster")
  # ,
  # arg3 = rep(ifelse(is.null(country_map), FALSE, TRUE), length(list_data))
  
  output <- 
    foreach(x = 1:length(list_data), 
            .combine='c', 
            arg1 = rep(exclude.area, length(list_data)), 
            arg2 = rep(buff_width, length(list_data)),
            arg4 = names_,
            arg5 = rep(alpha.hull, length(list_data)),
            arg6 = rep(convex.hull, length(list_data)),
            arg7 = rep(alpha, length(list_data)),
            arg8 = rep(buff.alpha, length(list_data)),
            arg9 = rep(method.less.than3, length(list_data))) %d% {
            # source("./R/IUCNeval.functionv11.R")
            if(show_progress)  utils::setTxtProgressBar(pb, x)
              
              # if(arg3) countr <- country_map
              # if(!arg3) countr <- FALSE
              
              
              res <- 
                .EOO.comp(XY = list_data[[x]], 
                exclude.area = arg1,
                buff_width = arg2, 
                country_map = country_map,
                Name_Sp = arg4,
                alpha.hull = arg5, 
                convex.hull = arg6, 
                alpha = arg7,
                buff.alpha = arg8, 
                method.less.than3 = arg9)
                
              names(res)[1] <- paste0(names(res)[1], "_" ,x)
              if(length(res)>1) names(res)[2] <- paste0(names(res)[2], "_" ,x)
              
              res
              
              
            }
  
  if(parallel) doParallel::stopImplicitCluster()
  
  if(length(output)==1) names(output) <- Name_Sp
  
  if(write_shp) {
    dir.create(file.path(paste(getwd(),"/shapesIUCN", sep="")), showWarnings = FALSE)
    output_spatial <- unlist(output[grep("spatial", names(output))])
    id_spatial <- as.numeric(unlist(lapply(strsplit(names(output_spatial), "_"), function(x) x[[2]])))
    for (i in 1:length(output_spatial)) {
      
      # if(!is.na(output[[i]][[1]])) {
        
        if(length(list.files(paste(getwd(),"/shapesIUCN", sep="")))>0){
          if(length(grep(paste(names(output)[i],"_EOO_poly", sep=""), unique(sub("....$", '', list.files(paste(getwd(),"/shapesIUCN", sep=""))))))>0) 
            {
            FILES <- list.files(paste(getwd(),"/shapesIUCN", sep=""), full.names = TRUE)
            file.remove(FILES[grep(paste(names(output)[i],"_EOO_poly", sep=""), FILES)])
            }
                                                                      }
        NAME <- names_[id_spatial[i]]
        output_spatial[[i]]@polygons[[1]]@ID <- "1"
        ConvexHulls_poly_dataframe <- SpatialPolygonsDataFrame(output_spatial[[i]], data=as.data.frame(names(output_spatial[[i]])))
        colnames(ConvexHulls_poly_dataframe@data) <- paste(substr(names_[id_spatial[i]], 0, 3), collapse = '')
        writeOGR(ConvexHulls_poly_dataframe,"shapesIUCN", paste(names_[id_spatial[i]],"_EOO_poly", sep=""), driver="ESRI Shapefile", overwrite_layer = TRUE)
      # }
    }
  }
  
  Results_short <- data.frame(EOO = unlist(output[grep("EOO", names(output))]))
  row.names(Results_short) <- names_
  
  if(write_results) write.csv(Results_short, paste(getwd(),"/", file.name, ".csv", sep=""))
  
  if(!export_shp) output <- Results_short
  
  output
}


.subpop.comp <- function(XY, Resol_sub_pop,
                        projEAC=crs("+proj=cea +lon_0=Central Meridian+lat_ts=Standard Parallel+x_0=False Easting+y_0=False Northing +ellps=WGS84")) {
    XY <- XY[,c(2,1)]
    coordEAC <- as.data.frame(matrix(unlist(rgdal::project(as.matrix(XY),proj=as.character(projEAC),inv=FALSE)), ncol=2))
    rownames(coordEAC) <- seq(1,nrow(coordEAC),1)
    
    p2=readWKT(paste("POINT(",mean(unique(coordEAC)[1,1])," ", mean(unique(coordEAC)[1,2]),")", sep=""))
    p2_Buffered1 <- gBuffer(p2, width = Resol_sub_pop*1000, id=1)
    if(nrow(unique(coordEAC))>1){
      for (LL in 2:nrow(unique(coordEAC))){
        p2=readWKT(paste("POINT(",mean(unique(coordEAC)[LL,1])," ", mean(unique(coordEAC)[LL,2]),")", sep=""))
        p2_Buffered <- gBuffer(p2, width = Resol_sub_pop*1000, id=LL)
        p2_Buffered1 <- gUnion(p2_Buffered1, p2_Buffered)
      }
    }
    splited_pol <- lapply(p2_Buffered1@polygons, slot, "Polygons")[[1]]
    NbeSubPop <- length(splited_pol)
    
    SubPopPoly <- SpatialPolygons(Srl=list(p2_Buffered1@polygons[[1]]), pO=as.integer(1), 
                                 proj4string=crs("+proj=cea +lon_0=Central Meridian+lat_ts=Standard Parallel+x_0=False Easting+y_0=False Northing +ellps=WGS84"))
    
    SubPopPoly <- sp::spTransform(SubPopPoly, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    
    OUTPUT <- list(NbeSubPop, SubPopPoly)
    names(OUTPUT) <- c("Number of subpopulation","subpop.poly")
    return(OUTPUT)
}


subpop.comp <- function(XY, Resol_sub_pop=NULL) {
  
  if(is.null(Resol_sub_pop)) stop("Resol_sub_pop is missing, please provide a value")
  
  if(any(is.na(XY[,c(1,2)]))) {
    length(which(rowMeans(is.na(XY[,1:2]))>0))
    unique(XY[which(rowMeans(is.na(XY[,1:2]))>0),3])
    print(paste("Skipping",length(which(rowMeans(is.na(XY[,1:2]))>0)) ,"occurrences because of missing coordinates for", 
                paste(as.character(unique(XY[which(rowMeans(is.na(XY[,1:2]))>0),3])), collapse=" AND ") ))
    XY <- XY[which(!is.na(XY[,1])),]
    XY <- XY[which(!is.na(XY[,2])),]
  }
  
  if(any(XY[,1]>180) || any(XY[,1]< -180)|| any(XY[,2]< -180) || any(XY[,2]>180)) stop("coordinates are outside of expected range")
  
  colnames(XY)[1:3] <- c("ddlat","ddlon","tax")
  XY$tax <- as.character(XY$tax)
  list_data <- split(XY, f = XY$tax)
  
  OUTPUT <- lapply(list_data, function(x) .subpop.comp(x, Resol_sub_pop=Resol_sub_pop))
  if(length(OUTPUT)==1) OUTPUT <- OUTPUT[[1]]
  return(OUTPUT)
}


.AOO.estimation <- function(coordEAC, 
                            cell_size = 2, 
                            nbe_rep = 0, 
                            poly_borders = NULL) {
  
  if(is.null(poly_borders)) crs_proj <- "+proj=longlat +datum=WGS84 +no_defs"
  if(!is.null(poly_borders)) crs_proj <- raster::crs(poly_borders)
  
  Corners <- rbind(c(min(coordEAC[,1]), 
                     max(coordEAC[,1])), 
                   c(min(coordEAC[,2]), 
                     max(coordEAC[,2])))
  
  if(nbe_rep==0) {
    Occupied_cells <- vector(mode = "numeric", length = 4)
    decal <- c(0,1,2,3)
    for (h in decal) {
      ext <- 
        raster::extent(floor(Corners[1,1])-h*(cell_size*1000/4)-2*cell_size*1000, 
                           floor(Corners[1,2])+h*(cell_size*1000/4)+2*cell_size*1000, 
                           floor(Corners[2,1])-h*(cell_size*1000/4)-2*cell_size*1000, 
                           floor(Corners[2,2])+h*(cell_size*1000/4)+2*cell_size*1000)
      r <- 
        raster::raster(ext, resolution=cell_size*1000, crs=crs_proj)
      r2_AOO <- 
        raster::rasterize(coordEAC[,1:2], r)
      OCC <- 
        length(which(!is.na(raster::values(r2_AOO))))
      Occupied_cells[h+1] <- OCC
      
      ### If only one occupied cell, stop the production of raster
      if(OCC==1) break
    }
    # h <- decal[which.min(Occupied_cells)]
    # Occupied_cells <- min(Occupied_cells)
  }
  
  if(nbe_rep>0) {
    Occupied_cells <- vector(mode = "numeric", length = nbe_rep)
    
    for (h in 1:nbe_rep) {
      rd.1 <- runif(1)*cell_size*1000
      rd.2 <- runif(1)*cell_size*1000
      
      ext = raster::extent(floor(Corners[1,1])-rd.1-2*cell_size*1000, floor(Corners[1,2])+rd.1+2*cell_size*1000, 
                           floor(Corners[2,1])-rd.2-2*cell_size*1000, floor(Corners[2,2])+rd.2+2*cell_size*1000)
      r = raster::raster(ext, resolution=cell_size*1000, crs=crs_proj)
      # r
      r2_AOO <- raster::rasterize(coordEAC[,1:2], r)
      OCC <- length(which(!is.na(raster::values(r2_AOO))))
      Occupied_cells[h] <- OCC
      # rd.1.vec <- c(rd.1.vec, rd.1)
      # rd.2.vec <- c(rd.2.vec, rd.2)
      if(OCC==1) break
    }
    
  }
  
  Occupied_cells <- Occupied_cells[Occupied_cells>0]
  Occupied_cells <- min(Occupied_cells)
  
  AOO <- Occupied_cells*cell_size*cell_size  ### AOO
  return(AOO)
  
}



AOO.computing <- function(XY, 
                          Cell_size_AOO = 2,
                          nbe.rep.rast.AOO = 0,
                          show_progress=TRUE, 
                          parallel = FALSE, 
                          NbeCores=2) {
  
  if(!any(class(XY)=="data.frame")) XY <- as.data.frame(XY)
  if(any(XY[,2]>180) || any(XY[,2]< -180)|| any(XY[,1]< -180) || any(XY[,1]>180)) stop("coordinates are outside of expected range")
  
  ## Equal Area cylindrical projection used for AOO estimation
  projEAC <- crs("+proj=cea +lon_0=Central Meridian+lat_ts=Standard Parallel+x_0=False Easting+y_0=False Northing +ellps=WGS84")
  
  coordEAC <- 
    data.frame(matrix(unlist(rgdal::project(as.matrix(XY[,c(2,1)]), 
                                            proj=as.character(projEAC),inv=FALSE)), 
                      ncol=2),
               tax = XY[,3])
  
  ## if any missing coordinates
  if(any(is.na(coordEAC[,c(1:2)]))) {
    print(paste("Skipping",length(which(rowMeans(is.na(coordEAC[,1:2]))>0)),
                "occurrences because of missing coordinates for", 
                paste(as.character(unique(coordEAC[which(rowMeans(is.na(coordEAC[,1:2]))>0),3])), collapse=" AND ") ))
    coordEAC <- coordEAC[which(!is.na(coordEAC[,1])),]
    coordEAC <- coordEAC[which(!is.na(coordEAC[,2])),]
  }
  
  coordEAC$tax <- as.character(coordEAC$tax)
  list_data <- split(coordEAC, f = coordEAC$tax)
  
  # if(show_progress) prog. <- "text"
  # if(!show_progress) prog. <- "none"
  
  if(show_progress)  pb <- 
    utils::txtProgressBar(min = 0, max = length(list_data), style = 3)
  
  if(parallel) {
    if("doParallel" %in% 
       rownames(installed.packages()) == FALSE) {stop("Please install doParallel package")}
    
    library(doParallel)
    
    registerDoParallel(NbeCores)
    message('doParallel running with ',
            NbeCores, ' cores')
    `%d%` <- `%dopar%`
  }else{
    `%d%` <- `%do%`
  }
  
  output <- 
    foreach(x=1:length(list_data), .combine='c', 
            arg1 = rep(Cell_size_AOO, length(list_data)), 
            arg2 = rep(nbe.rep.rast.AOO, length(list_data))) %d% {
              # source("./R/IUCNeval.functionv11.R")
              if(show_progress)  utils::setTxtProgressBar(pb, x)
              # utils::setTxtProgressBar(pb, x)
              
              res <- .AOO.estimation(coordEAC = list_data[[x]], 
                                     cell_size = arg1, nbe_rep = arg2)
              
              res
            }
  
  
  # OUTPUT <- plyr::llply(list_data, .fun=function(x) {
  #   source("./R/IUCNeval.functionv11.R")
  #   # Cell_size_AOO=2
  #   # nbe.rep.rast.AOO=NULL
  #   .AOO.estimation(x, 
  #                   cell_size = Cell_size_AOO, 
  #                   nbe_rep = nbe.rep.rast.AOO, 
  #                   poly_borders = NULL) # , verbose=verbose
  # }
  # , .progress = prog., .parallel=parallel, 
  # .paropts = list(.packages = c("ConR"), 
  #                 .combine = c, 
  #                 Cell_size_AOO = Cell_size_AOO,
  #                 nbe_rep = nbe.rep.rast.AOO))
  
  
  if(parallel) stopImplicitCluster()
  
  res <- unlist(output)
  names(res) <- names(list_data)
  
  return(res)
}





.IUCN.comp <- function(DATA, 
                       poly_borders=NULL, 
                       Cell_size_AOO=2, 
                       Cell_size_locations=10, 
                       Resol_sub_pop=5, 
                       method_locations=c("fixed_grid"), 
                       Rel_cell_size=0.05, 
                       protec.areas=NULL, exclude.area=FALSE, 
                       method_protected_area="no_more_than_one",
                       ID_shape_PA="WDPA_PID", 
                       buff_width=0.1, 
                       NamesSp="species1", 
                       write_shp=FALSE, 
                       file_name=NULL, 
                       add.legend=TRUE, 
                       DrawMap=TRUE, 
                       map_pdf=FALSE, 
                       draw.poly.EOO=TRUE, 
                       SubPop=TRUE, 
                       MinMax,
                       alpha=1, 
                       buff.alpha=0.1, 
                       method.range="convex.hull", 
                       nbe.rep.rast.AOO=0, 
                       verbose=TRUE, 
                       showWarnings=TRUE) {
  
  ### Getting by default land map if poly_borders is not provided
  if(is.null(poly_borders)) {
    data('land', package='ConR', envir=environment()) 
    land <- get("land", envir=environment()) 
    # data(land, envir = environment())
    poly_borders=land
  }
  
  ### cropping poly_borders according to range of occurrences shapefile for producing lighter maps
  if(DrawMap) {
      full_poly_borders <- poly_borders
      if(!is.null(poly_borders)) poly_borders <- raster::crop(poly_borders, extent(MinMax)+30)
  }
  

  
  ## Equal Area cylindrical projection used for AOO estimation
  projEAC <- crs("+proj=cea +lon_0=Central Meridian+lat_ts=Standard Parallel+x_0=False Easting+y_0=False Northing +ellps=WGS84")
  
  ## Initialization of data frame for stocking results
  if(is.null(protec.areas)){
    Results <- as.data.frame(matrix(NA,9,1))
    colnames(Results) <- NamesSp
    rownames(Results) <- c("EOO","AOO","Nbe_unique_occ.", "Nbe_subPop", "Nbe_loc", "Category_CriteriaB", "Category_code",
                           "Category_AOO","Category_EOO")
  }else{
    Results <- as.data.frame(matrix(NA,11,1))
    colnames(Results) <- NamesSp
    rownames(Results) <- c("EOO","AOO","Nbe_unique_occ.", "Nbe_subPop", "Nbe_loc", "Nbe_loc_PA", 
                           "Category_CriteriaB", "Category_code",  "Ratio_occ_within_PA",
                           "Category_AOO","Category_EOO")
  }
  
  # if(verbose) cat(" ",paste("Evaluation of", as.character(NamesSp)))
  
  XY <- DATA[,c(2:1)]
  ### Projection into equal area cylindrical
  coordEAC <- as.data.frame(matrix(unlist(rgdal::project(as.matrix(XY),proj=as.character(projEAC),inv=FALSE)), ncol=2))
  rownames(coordEAC) <- seq(1,nrow(coordEAC),1)
  
  ##########################################################################################
  ##############################  Sub-populations estimations ############################## 
  
  if(SubPop) {
    subpop_stats <- subpop.comp(DATA, Resol_sub_pop=Resol_sub_pop)
    SubPopPoly <- subpop_stats[[2]]
    NbeSubPop <- subpop_stats[[1]]
  }
  
  
  ##########################################################################################
  ##############################  Estimations of number of Locations ####################### 
  
  ## range of lat and long
  Corners <- rbind(c(min(XY[,1]), max(XY[,1])), c(min(XY[,2]), max(XY[,2])))
  
  ## geographical distances for all pairs of occurrences
  if(nrow(coordEAC)>1) pairwise_dist <- dist(coordEAC,  upper = F)
  
  ## resolution definition
  if(any(method_locations=="fixed_grid")) Resolution <- Cell_size_locations*1000
  if(any(method_locations=="sliding scale")){
    if(nrow(coordEAC)>1) {Resolution <- max(pairwise_dist)*Rel_cell_size
    }else{
      Resolution <- 10000
    }
  }
  
  ## Convert resolution of grid in km into decimal degrees after projection, taking mean X and Y of occurrences
  border_to_center <- as.data.frame(matrix(NA, 2, 2))
  border_to_center[1,] <- c(mean(coordEAC[,1]), mean(coordEAC[,2]))
  border_to_center[2,] <- c( border_to_center[1,1]+Resolution,  border_to_center[1,2])
  DIST_circle <- matrix(unlist(rgdal::project(as.matrix(border_to_center),proj=as.character(projEAC),inv =T)), ncol=2)
  cell_size_deg <- abs((DIST_circle[1,1]-DIST_circle[2,1]))
  
  
  if(!is.null(protec.areas)) { ### Taking into account Protected Areas if provided
    DATA_SF <- as.data.frame(XY)
    colnames(DATA_SF) <- c("ddlon","ddlat")
    coordinates(DATA_SF) <-  ~ddlon+ddlat
    crs(DATA_SF) <- crs(protec.areas)
    Links_NatParks <- over(DATA_SF, protec.areas)
    
    if(length(which(!is.na(Links_NatParks[,1])))!=0){
      if(method_protected_area=="no_more_than_one"){
        ## if method is 'no_more_than_one' the number of location is the number of occupied protected areas
        LocNatParks <- 
          length(unique(Links_NatParks[which(!is.na(Links_NatParks[,1])), ID_shape_PA]))
      }else{ #### If method for accounting Protected areas should superimpose a grid
        if(length(which(!is.na(Links_NatParks[,1])))>1) {
          
          ### By default, rasters with different starting position are created
          ## by doing translation of 1/4 increment north and east
          LocNatParks <- c()
          decal <- c(0,1,2,3)
          
          for (h in decal) {
            ext = extent(floor(min(Corners[,1]))-3+h*(cell_size_deg/4), floor(max(Corners[,2]))+3+h*(cell_size_deg/4), 
                         floor(min(Corners[,1]))-3+h*(cell_size_deg/4), floor(max(Corners[,2]))+3+h*(cell_size_deg/4))
            r = raster(ext, resolution=cell_size_deg,crs=crs(poly_borders))
            r2_PA <- rasterize(XY[which(!is.na(Links_NatParks[,1])),], r)
            LOC <- length(which(!is.na(values(r2_PA))))
            LocNatParks <- c(LocNatParks, LOC)
            
            ### stop to produce rasters with different starting position if number of locations is 1 or >10
            if(LOC>10 || LOC==1) break 
          }
          
          ### The minimum number is chosen
          h <- decal[which.min(LocNatParks)]
          LocNatParks <- min(LocNatParks)
          
          ext = extent(floor(min(Corners[,1]))-3+h*(cell_size_deg/4), floor(max(Corners[,2]))+3+h*(cell_size_deg/4), 
                       floor(min(Corners[,1]))-3+h*(cell_size_deg/4), floor(max(Corners[,2]))+3+h*(cell_size_deg/4))
          r = raster(ext, resolution=cell_size_deg,crs=crs(poly_borders))
          r2_PA <- rasterize(XY[which(!is.na(Links_NatParks[,1])),], r)
          
        }else{LocNatParks <- 1
              ext = extent(floor(min(Corners[,1]))-3, floor(max(Corners[,2]))+3, 
                           floor(min(Corners[,1]))-3, floor(max(Corners[,2]))+3)
              r = raster(ext, resolution=cell_size_deg,crs=crs(poly_borders))
              r2_PA <- rasterize(XY[which(!is.na(Links_NatParks[,1])),], r)
        }
      }
    }else{LocNatParks <- 0}
    
    if(length(which(is.na(Links_NatParks[,1])))!=0){
      if(length(which(is.na(Links_NatParks[,1])))>1){
        pairwise_dist <- dist(coordEAC[which(is.na(Links_NatParks[,1])),],  upper = F)
        
        ### By default, rasters with different starting position are created
        ## by doing translation of 1/4 increment north and east
        LocOutNatParks <- c()
        decal <- c(0,1,2,3)
        for (h in decal) {
          ext = extent(floor(min(Corners[,1]))-3+h*(cell_size_deg/4), floor(max(Corners[,2]))+3+h*(cell_size_deg/4), 
                       floor(min(Corners[,1]))-3+h*(cell_size_deg/4), floor(max(Corners[,2]))+3+h*(cell_size_deg/4))
          r = raster(ext, resolution=cell_size_deg,crs=crs(poly_borders))
          r2 <- rasterize(XY[which(is.na(Links_NatParks[,1])),], r)
          LOC <- length(which(!is.na(values(r2))))
          LocOutNatParks <- c(LocOutNatParks, LOC)
          
          ### stop to produce rasters with different starting position if number of locations is 1 or >10
          if(LOC>10 || LOC==1) break
        }
        
        ### The minimum number is chosen
        h <- decal[which.min(LocOutNatParks)]
        LocOutNatParks <- min(LocOutNatParks)
        ext = extent(floor(min(Corners[,1]))-3+h*(cell_size_deg/4), floor(max(Corners[,2]))+3+h*(cell_size_deg/4), 
                     floor(min(Corners[,1]))-3+h*(cell_size_deg/4), floor(max(Corners[,2]))+3+h*(cell_size_deg/4))
        r = raster(ext, resolution=cell_size_deg,crs=crs(poly_borders))
        r2 <- rasterize(XY[which(is.na(Links_NatParks[,1])),], r)
      }else{LocOutNatParks <- 1
            ext = extent(floor(min(Corners[,1]))-3, floor(max(Corners[,2]))+3, 
                         floor(min(Corners[,1]))-3, floor(max(Corners[,2]))+3)
            r = raster(ext, resolution=cell_size_deg,crs=crs(poly_borders))
            r2 <- rasterize(XY[which(is.na(Links_NatParks[,1])),], r)
      }
    }else{LocOutNatParks <- 0}
    
  }else{
    if(nrow(coordEAC)>1) {
      
      ### If protected areas are not taken into account and the number of occurrence i > 1
      if(max(pairwise_dist) > Resolution) {
        
        ### By default, rasters with different starting position are created
        ## by doing translation of 1/4 increment north and east
        Locations <- c()
        decal <- c(0,1,2,3)
        for (h in decal) {
          ext = extent(floor(min(Corners[,1]))-3+h*(cell_size_deg/4), floor(max(Corners[,2]))+3+h*(cell_size_deg/4), 
                       floor(min(Corners[,1]))-3+h*(cell_size_deg/4), floor(max(Corners[,2]))+3+h*(cell_size_deg/4))
          r = raster(ext, resolution=cell_size_deg, crs=crs(poly_borders))
          r2 <- rasterize(XY, r)
          LOC <- length(which(!is.na(values(r2))))
          Locations <- c(Locations, LOC)
          
          ### stop to produce rasters with different starting position if number of locations is 1 or >10
          if(LOC>10 || LOC==1) break
        }
        
        ### The minimum number is chosen
        h <- decal[which.min(Locations)]
        Locations <- min(Locations)
        ext = extent(floor(min(Corners[,1]))-3+h*(cell_size_deg/4), floor(max(Corners[,2]))+3+h*(cell_size_deg/4), 
                     floor(min(Corners[,1]))-3+h*(cell_size_deg/4), floor(max(Corners[,2]))+3+h*(cell_size_deg/4))
        r = raster(ext, resolution=cell_size_deg,crs=crs(poly_borders))
        r2 <- rasterize(XY, r)
        
      }else{Locations <- 1
            ext = extent(floor(min(Corners[,1]))-3, floor(max(Corners[,2]))+3, 
                         floor(min(Corners[,1]))-3, floor(max(Corners[,2]))+3)
            r = raster(ext, resolution=cell_size_deg,crs=crs(poly_borders))
            r2 <- rasterize(XY, r)
      }
    }else{
      Locations <- 1
      ext = extent(floor(min(Corners[,1]))-3, floor(max(Corners[,2]))+3, 
                   floor(min(Corners[,1]))-3, floor(max(Corners[,2]))+3)
      r = raster(ext, resolution=cell_size_deg,crs=crs(poly_borders))
      r2 <- rasterize(XY, r)
    }
  }
  
  
   
  if(nrow(unique(XY))>2) { ### if more than 2 uniques occurrences
    
  ##########################################################################################
  ##############################  EOO estimation              ##############################
    
    EOO_ <- EOO.computing(XY[,c(2,1)], exclude.area=exclude.area, country_map=poly_borders, Name_Sp=NamesSp, 
                          buff_width=buff_width, export_shp=TRUE,
                          alpha=alpha, buff.alpha=buff.alpha, method.range=method.range, write_results=FALSE) # , verbose=FALSE
    p1 <- EOO_[[2]]
    EOO <- EOO_[[1]]
    
  ##########################################################################################
  ##############################  AOO estimation              ##############################
    
    AOO <- 
      .AOO.estimation(coordEAC, 
                      cell_size = Cell_size_AOO, 
                      nbe_rep = nbe.rep.rast.AOO, 
                      poly_borders = poly_borders)
    
    if(EOO<AOO) EOO <- AOO ### If EOO is < AOO, EOO is put equal to AOO
    
    ## recording results
    Results["EOO",1] <- as.numeric(EOO)
    Results["AOO",1] <- as.numeric(AOO)
    if(SubPop) Results["Nbe_subPop",1] <- NbeSubPop
    Results["Nbe_unique_occ.",1] <- nrow(unique(XY))
    
    if(!is.null(protec.areas)) Results["Nbe_loc",1] <- LocNatParks + LocOutNatParks
    if(is.null(protec.areas)) Results["Nbe_loc",1] <- Locations
    
    if(!is.null(protec.areas)) {Results["Nbe_loc_PA",1] <- LocNatParks}
    
    if(!is.null(protec.areas)) Results["Ratio_occ_within_PA",1] <- round(length(which(!is.na(Links_NatParks[,1])))/nrow(Links_NatParks)*100,1)
    Nbe_Loc <- as.numeric(Results["Nbe_loc",1])
    
    ##################
    ### Criteria B assessment following IUCN thresholds
    if(EOO<20000){
      Rank_EOO <- 3
      if(EOO<5000){
        Rank_EOO <- 2
        if(EOO<100){
          Rank_EOO <- 1
        }}}else(Rank_EOO <- 4)
    
    if(AOO<2000){
      Rank_AOO <- 3
      if(AOO<500){
        Rank_AOO <- 2
        if(AOO>10){
          Rank_AOO <-1
        }}}else{Rank_AOO <- 4}
    
    if(Nbe_Loc<=10){
      Rank_Loc <-3
      if(Nbe_Loc <=5){
        Rank_Loc <-2
        if(Nbe_Loc==1){
          Rank_Loc <- 1
        }}}else{Rank_Loc <- 4}
    
    Rank_B1a <- max(Rank_EOO, Rank_Loc)
    Rank_B2a <- max(Rank_AOO, Rank_Loc)
    Rank_CriteriaB <- min(Rank_B1a, Rank_B2a)
    
    if(Rank_CriteriaB==1) Cat <- "CR"
    if(Rank_CriteriaB==2) Cat <- "EN"
    if(Rank_CriteriaB==3 && Nbe_Loc>0 && Nbe_Loc<11) Cat <- "VU"

    if(Rank_CriteriaB>3 && Nbe_Loc>=0) Cat <- "LC or NT" ### 

    if(Rank_B1a>Rank_B2a) Cat_Code <- paste(Cat,"B2a")
    if(Rank_B1a<Rank_B2a) Cat_Code <- paste(Cat,"B1a")
    if(Rank_B1a==Rank_B2a) Cat_Code <- paste(Cat,"B1a+B2a")
    
    if(!is.null(protec.areas)) {
      if(as.numeric(Results["Ratio_occ_within_PA",1])==100){
        Results["Category_CriteriaB",1] <- "LC or NT"
        Results["Category_code",1] <- Cat_Code 
      }else{
        Results["Category_CriteriaB",1] <- Cat
        Results["Category_code",1] <- Cat_Code
      }
    }else{
      Results["Category_CriteriaB",1] <- Cat
      Results["Category_code",1] <- Cat_Code
    }
    
    if(Rank_B2a==1) Results["Category_AOO",1] <- "CR"
    if(Rank_B2a==2) Results["Category_AOO",1] <- "EN"
    if(Rank_B2a==3) Results["Category_AOO",1] <- "VU"
    if(Rank_B2a>3) Results["Category_AOO",1] <- "LC or NT"
    
    if(Rank_B1a==1) Results["Category_EOO",1] <- "CR"
    if(Rank_B1a==2) Results["Category_EOO",1] <- "EN"
    if(Rank_B1a==3) Results["Category_EOO",1] <- "VU"
    if(Rank_B1a>3) Results["Category_EOO",1] <- "LC or NT"
    
  }else{ ### if less than 3 uniques occurrences
    
    p1 <- NULL ## EOO shapefile is NULL
    
    
    if(nrow(coordEAC)==2) { ## if two uniques occurrences
      if(pairwise_dist<=Resolution) {
        AOO <- Cell_size_AOO*Cell_size_AOO ## 1 occupied cell if distance <= resolution
      }else{
        AOO <- 2*Cell_size_AOO*Cell_size_AOO ## 2 occupied cells if distance > resolution
      }
    }else{
      AOO <- Cell_size_AOO*Cell_size_AOO ## 1 occupied cell if one unique occurrence
    }
    
    Results["AOO",1] <- AOO
    if(SubPop) Results["Nbe_subPop",1] <- NbeSubPop
    Results["Nbe_unique_occ.",1] <- nrow(unique(XY))
    
    if(!is.null(protec.areas)) Results["Ratio_occ_within_PA",1] <- round(length(which(!is.na(Links_NatParks[,1])))/nrow(Links_NatParks)*100,2)
    
    if(is.null(protec.areas)) Results["Nbe_loc",1] <- Locations
    if(!is.null(protec.areas)) Results["Nbe_loc",1] <- LocNatParks + LocOutNatParks
    if(!is.null(protec.areas)) Results["Nbe_loc_PA",1] <- LocNatParks
    
    
    if(!is.null(protec.areas)){
      if(as.numeric(Results["Ratio_occ_within_PA",1])==100){
        ### If all occurences are found within protected areas, the species is considered as not threatened
        Results["Category_CriteriaB",1] <- "LC or NT"
      }else{
        if(as.numeric(Results["AOO",1]) < 10 & as.numeric(Results["Nbe_loc",1])==1) {
          Results["Category_AOO",1] <- Results["Category_CriteriaB",1] <- "CR"
          
        }else{
          if(as.numeric(Results["AOO",1]) < 500 & as.numeric(Results["Nbe_loc",1])<6) {
            Results["Category_AOO",1] <- Results["Category_CriteriaB",1] <- "EN"
            
          }else{
            if(as.numeric(Results["AOO",1]) < 2000 & as.numeric(Results["Nbe_loc",1])<11) {
              Results["Category_AOO",1] <- Results["Category_CriteriaB",1] <- "VU"
              
            }else{
              Results["Category_AOO",1] <- Results["Category_CriteriaB",1] <- "LC or NT"
            }
          }
        }
      }
    }else{
      if(as.numeric(Results["AOO",1]) < 10 & as.numeric(Results["Nbe_loc",1])==1) {
        Results["Category_AOO",1] <- Results["Category_CriteriaB",1] <- "CR"
      }else{
        if(as.numeric(Results["AOO",1]) < 500 & as.numeric(Results["Nbe_loc",1])<6) {
          Results["Category_AOO",1] <- Results["Category_CriteriaB",1] <- "EN"
          
        }else{
          if(as.numeric(Results["AOO",1]) < 2000 & as.numeric(Results["Nbe_loc",1])<11) {
            Results["Category_AOO",1] <- Results["Category_CriteriaB",1] <- "VU"
            
          }else{
            Results["Category_AOO",1] <- Results["Category_CriteriaB",1] <- "LC or NT"
          }
        }
      }
    }
    
    if(is.na(Results["Category_AOO",1])) {
      if(as.numeric(Results["AOO",1]) < 10 & as.numeric(Results["Nbe_loc",1])==1) {
        Results["Category_AOO",1] <- "CR"
        
      }else{
        if(as.numeric(Results["AOO",1]) < 500 & as.numeric(Results["Nbe_loc",1])<6) {
          Results["Category_AOO",1] <- "EN"
          
        }else{
          if(as.numeric(Results["AOO",1]) < 2000 & as.numeric(Results["Nbe_loc",1])<11) {
            Results["Category_AOO",1] <- "VU"
            
          }else{
            Results["Category_AOO",1] <- "LC or NT"
          }
        }
      }
    }
    
    Results["Category_code",1] <- paste(Results["Category_CriteriaB",1],"B2a")
    if(showWarnings) warning(paste("EOO statistic is not computed for", NamesSp,"because there is less than 3 records"))
  } ## End less than 3 records
  
  if(DrawMap) {
    
    ## pdf or png format initialization
    if(!map_pdf) {
      if(!is.null(file_name)) {
        NAME_FILE <- paste(file_name, gsub(" ",replacement = "_", as.character(NamesSp)) , sep="")
      }else{
        NAME_FILE <- paste("IUCN_", gsub(" ",replacement = "_", as.character(NamesSp)), sep="")
      }
      FILE_NAME <- ifelse(!is.null(file_name), file_name, "IUCN_")
      dir.create(file.path(paste(getwd(),paste("/",FILE_NAME,"_results_map", sep=""), sep="")), showWarnings = FALSE)
    }
    
    if(!map_pdf) png(paste(file.path(paste(getwd(),paste("/",FILE_NAME,"_results_map", sep=""), sep="")),"/",NAME_FILE,".png", sep=""), width=2000, height=2000)
    
    
    ### Layout of the map
    par(mar=c(10, 12, 10, 2), xpd=FALSE, las=1)
    if(add.legend & !any(colnames(DATA)=="coly")) nf <- layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), c(4,1.5), c(4,1.5))
    if(any(colnames(DATA)=="coly") & add.legend) nf <- layout(matrix(c(1,1,1,1,1,1,2,3,4), 3, 3, byrow = TRUE), c(2,1.5,1.5), c(4,1.5,1.5))
    
    
    ### Mapping 
    if(!is.null(protec.areas)){
      if(LocOutNatParks==0){
        plot(poly_borders, xlim=c(range(XY[,1])[1]-1, range(XY[,1])[2]+1), ylim=c(range(XY[,2])[1]-1, range(XY[,2])[2]+1), axes=FALSE, xlab="", ylab="")
      }else{
        r2_pol <- rasterToPolygons(r2, fun=NULL, n=4, na.rm=TRUE, digits=6, dissolve=FALSE)
        if(LocOutNatParks==1){
          plot(r2_pol, col=rgb(red=1, green=0, blue=0, alpha=0.2), 
               xlim=c(range(XY[,1])[1]-1, range(XY[,1])[2]+1), ylim=c(range(XY[,2])[1]-1, range(XY[,2])[2]+1))
        }else{
          plot(r2_pol, col=rgb(red=1, green=0, blue=0, alpha=0.2), 
               xlim=c(range(XY[,1])[1]-1, range(XY[,1])[2]+1), ylim=c(range(XY[,2])[1]-1, range(XY[,2])[2]+1))
        }
      }
    }else{
      r2_pol <- rasterToPolygons(r2, fun=NULL, n=4, na.rm=TRUE, digits=6, dissolve=FALSE)
      plot(r2_pol, col=rgb(red=1, green=0, blue=0, alpha=0.2), 
           xlim=c(range(XY[,1])[1]-1, range(XY[,1])[2]+1), ylim=c(range(XY[,2])[1]-1, range(XY[,2])[2]+1))
    }
    
    if(SubPop) plot(SubPopPoly, add=T, border="black", lwd=2, lty=1)
    
    if(!is.null(protec.areas)){
      if(LocNatParks>0){
        if(method_protected_area!="no_more_than_one"){
          r2_PA_pol <- rasterToPolygons(r2_PA, fun=NULL, n=4, na.rm=TRUE, digits=6, dissolve=FALSE)
          plot(r2_PA_pol, add=T, col=rgb(red=0, green=0, blue=1, alpha=0.2))                 
        }
      }
    }
    
    if(!is.null(p1) & draw.poly.EOO) plot(p1, add=T, col=rgb(red=0.2, green=0.2, blue=0.2, alpha=0.1))
    
    plot(poly_borders, axes=FALSE, lty=1, add=T, lwd=1)
    
    if(!is.null(protec.areas)) plot(protec.areas, add=T, col=rgb(red=0.2, green=0.2, blue=0.2, alpha=0.05), lty=2)
    
    if(!is.null(protec.areas)){
      colnames(XY) <- c("ddlon","ddlat")
      XY_sp <- XY[which(is.na(Links_NatParks[,1])),]
      if(nrow(XY_sp)>0){
        coordinates(XY_sp) <-  ~ddlon+ddlat
        plot(XY_sp, pch=19, cex=2, col="black", add=T)
      }
      XY_sp <- XY[which(!is.na(Links_NatParks[,1])),]
      if(nrow(XY_sp)>0){
        coordinates(XY_sp) <-  ~ddlon+ddlat
        plot(XY_sp, pch=19, cex=2, col="blue", add=T)            
      }
    }else{
      colnames(XY) <- c("ddlon","ddlat")
      XY_sp <- XY
      coordinates(XY_sp) <-  ~ddlon+ddlat
      plot(XY_sp, pch=19, cex=2, col="black", add=T)
    }
    
    axis(1, outer=FALSE, cex.axis=3, tick = FALSE, line=1.5)  #pos=min(range(XY[,2]))-2)
    axis(1, outer=FALSE,labels=FALSE, cex.axis=3, tick = TRUE, line=0)  #pos=min(range(XY[,2]))-2)
    axis(2, outer=FALSE, cex.axis=3, tick = FALSE, line=1.5)  #pos=min(range(XY[,2]))-2)
    axis(2, outer=FALSE,labels=FALSE, cex.axis=3, tick = TRUE, line=0)  #pos=min(range(XY[,2]))-2)
    box()
    
    if(Results["Nbe_loc",1]>1) {
      xlim <- par("xaxp")[1:2]
      xlim <- abs(xlim[1]-xlim[2])
      border_to_center <- as.data.frame(matrix(NA, 2, 2))
      border_to_center[,1] <- c(xlim/10, 0)
      border_to_center[,2] <- c(0,0)
      scaleBAR <- round(matrix(unlist(rgdal::project(as.matrix(border_to_center), proj=as.character(projEAC),inv =F)), ncol=2)/1000,0)[1,1]
    }else{
      scaleBAR <- Resolution/1000
    }
    scalebar(scaleBAR, type="bar", below="kilometres", cex=2.5)
    
    
    mtext(NamesSp, side=3, cex=3, line=3)
    if(any(colnames(DATA)=="higher.tax.rank")) mtext(DATA[which(DATA[,3]==NamesSp),"higher.tax.rank"][1], side=3, cex=3, line=0.4)
    
    if(add.legend) {
      par(mar=c(1,1,1,1), xpd=T)
      plot(1:10, 1:10, type="n", bty='n', xaxt='n', yaxt='n')
      if(is.null(protec.areas)){
        legend(1,10,  c(paste("EOO=", ifelse(!is.na(Results["EOO",1]), round(as.numeric(Results["EOO",1]),1), NA), "km2"),
                        paste("AOO (grid res.",Cell_size_AOO,"km)=", format(Results["AOO",1], scientific = 5),"km2"),
                        paste("Number of unique occurrences=", Results["Nbe_unique_occ.",1]),
                        paste("Number of sub-populations (radius",Resol_sub_pop,"km)=", Results["Nbe_subPop",1]),
                        paste("Number of locations (grid res.:",round(Resolution/1000,1)," km)","=", Results["Nbe_loc",1]),
                        paste("IUCN category according to criterion B:", Results["Category_CriteriaB",1])), cex=3.5,bg = grey(0.9))
      }
      if(!is.null(protec.areas)){
        legend(1,10,  c(paste("EOO=", ifelse(!is.na(Results["EOO",1]), round(as.numeric(Results["EOO",1]),1), NA),"km2"),
                        paste("AOO (grid res.",Cell_size_AOO,"km)=", format(Results["AOO",1], scientific = 5),"km2"),
                        paste("Number of unique occurrences=", Results["Nbe_unique_occ.",1]),
                        paste("Number of sub-populations (radius",Resol_sub_pop,"km)=",Results["Nbe_subPop",1]),
                        paste("Number of locations (grid res.:",round(Resolution/1000,1)," km)","=", Results["Nbe_loc",1]),
                        paste("Number of occupied protected areas=", Results["Nbe_loc_PA",1]),
                        paste("IUCN category according to criterion B:", Results["Category_CriteriaB",1]),
                        paste("Proportion of occurences within protected areas"), Results["Ratio_occ_within_PA",1]), cex=3.5,bg = grey(0.9))
      }
      par(mar=c(4,1,1,1))
      plot(full_poly_borders, lty=1, lwd=1,axes=FALSE)
      points(XY[,1],XY[,2], pch=8, cex=2, col="red") 
    }
    
    if(any(colnames(DATA)=="coly") & add.legend) {
      par(mar=c(12,6,1,2), las=2, yaxs="r", xpd=FALSE)
      subdata <- DATA[which(DATA[,"tax"]==NamesSp),"coly"]
      if((sum(subdata, na.rm=T))>0) {
        plot(table(subdata), col="grey", ylab=" ", xlab=" ", cex.lab=4, cex.axis=4, axes=F)
        axis(1, outer=FALSE, cex.axis=3, tick = FALSE, line=1.5)  #pos=min(range(XY[,2]))-2)
        axis(1, outer=FALSE,labels=FALSE, cex.axis=3, tick = TRUE, line=0)  #pos=min(range(XY[,2]))-2)
        axis(2, outer=FALSE, cex.axis=3, tick = FALSE, line=1.5)  #pos=min(range(XY[,2]))-2)
        axis(2, outer=FALSE,labels=FALSE, cex.axis=3, tick = TRUE, line=0)  #pos=min(range(XY[,2]))-2)            
      }
    }
    
    if(!map_pdf) dev.off()
  } # end draw map
  
  if(write_shp) {
    
    dir.create(file.path(paste(getwd(),"/shapesIUCN", sep="")), showWarnings = FALSE)
    
    if(!is.null(p1)) {
      if(length(list.files(paste(getwd(),"/shapesIUCN", sep="")))>0){
        if(length(grep(paste(NamesSp,"_EOO_poly", sep=""), unique(sub("....$", '', list.files(paste(getwd(),"/shapesIUCN", sep=""))))))>0) {
          FILES <- list.files(paste(getwd(),"/shapesIUCN", sep=""), full.names = TRUE)
          file.remove(FILES[grep(paste(NamesSp,"_EOO_poly", sep=""), FILES)])
        }
      }
      NAME <- names(p1)
      p1@polygons[[1]]@ID <- "1"
      ConvexHulls_poly_dataframe <- SpatialPolygonsDataFrame(p1, data=as.data.frame(names(p1)))
      colnames(ConvexHulls_poly_dataframe@data) <- paste(substr(unlist(strsplit(NamesSp, "[ ]")), 0, 3), collapse = '')
      writeOGR(ConvexHulls_poly_dataframe,"shapesIUCN",paste(NamesSp,"_EOO_poly", sep=""),driver="ESRI Shapefile")
    }
    
    if(SubPop) {
      if(length(list.files(paste(getwd(),"/shapesIUCN", sep="")))>0){
        if(length(grep(paste(NamesSp,"_subpop_poly", sep=""), unique(sub("....$", '', list.files(paste(getwd(),"/shapesIUCN", sep=""))))))>0) {
          FILES <- list.files(paste(getwd(),"/shapesIUCN", sep=""), full.names = TRUE)
          file.remove(FILES[grep(paste(NamesSp,"_subpop_poly", sep=""), FILES)])
        }
      }
      NAME <- names(SubPopPoly)
      SubPopPoly@polygons[[1]]@ID <- "1"
      ConvexHulls_poly_dataframe <- SpatialPolygonsDataFrame(SubPopPoly, data=as.data.frame(names(SubPopPoly)))
      colnames(ConvexHulls_poly_dataframe@data) <- paste(substr(unlist(strsplit(NamesSp, "[ ]")), 0, 3), collapse = '')
      writeOGR(ConvexHulls_poly_dataframe,"shapesIUCN",paste(NamesSp,"_subpop_poly", sep=""),driver="ESRI Shapefile")      
    }
  }
  
  
  if(SubPop) {
    OUTPUT <- list(Results, p1, SubPopPoly)
    names(OUTPUT) <- c("Results","spatialPoly_EOO","spatialPoly_subpop")
  } else{
    OUTPUT <- list(Results, p1)
    names(OUTPUT) <- c("Results","spatialPoly_EOO")    
  }

  return(OUTPUT)
}


IUCN.eval <- function (DATA, 
                       country_map = NULL, 
                       Cell_size_AOO = 2, 
                       Cell_size_locations = 10, 
                       Resol_sub_pop = 5, 
                       method_locations = "fixed_grid", 
                       Rel_cell_size = 0.05, 
                       DrawMap = FALSE, 
                       add.legend = TRUE, 
                       file_name = NULL, 
                       export_shp = FALSE, 
                       write_shp = FALSE, 
                       write_results=TRUE, 
                       protec.areas = NULL, 
                       map_pdf = FALSE, 
                       draw.poly.EOO = TRUE, 
                       exclude.area = FALSE, 
                       method_protected_area = "no_more_than_one", 
                       ID_shape_PA = "WDPA_PID", 
                       buff_width = 0.1, 
                       SubPop=TRUE, 
                       alpha=1, 
                       buff.alpha=0.1, 
                       method.range="convex.hull", 
                       nbe.rep.rast.AOO=0,
                       showWarnings=TRUE, # verbose=TRUE,
                       write_file_option="excel", 
                       parallel=F, 
                       NbeCores=2) {
  
  if(class(DATA)[1]=="spgeoIN") {
    DATA_2 <- cbind(DATA$species_coordinates, DATA$identifier)
    DATA <- DATA_2[,c(2,1,3)]
  }
  colnames(DATA)[1:3] <- c("ddlat","ddlon","tax")
  
  if(is_tibble(DATA)) DATA <- as.data.frame(DATA)
  
  if(any(is.na(DATA[,1:2]))) {
    length(which(rowMeans(is.na(DATA[,1:2]))>0))
    unique(DATA[which(rowMeans(is.na(DATA[,1:2]))>0),3])
    print(paste("Skipping",length(which(rowMeans(is.na(DATA[,1:2]))>0)) ,"occurrences because of missing coordinates for", 
                paste(as.character(unique(DATA[which(rowMeans(is.na(DATA[,1:2]))>0),3])), collapse=" AND ") ))
    DATA <- DATA[which(!is.na(DATA[,1])),]
    DATA <- DATA[which(!is.na(DATA[,2])),]
  }
  
  if(is.factor(DATA[,"tax"])) DATA[,"tax"] <- as.character(DATA[,"tax"])
  
  
  if(!is.numeric(DATA[,1]) || !is.numeric(DATA[,2]) ) {
    if(!is.double(DATA[,1]) || !is.double(DATA[,2])) stop("coordinates in DATA should be numeric")
  }
  
  if(any(DATA[,1]>180) || any(DATA[,1]< -180)|| any(DATA[,2]< -180) || any(DATA[,2]>180)) stop("coordinates are outside of expected range")
  if(!is.null(country_map)) if(!class(country_map)=="SpatialPolygonsDataFrame") stop("Country_map should be a spatialpolygondataframe")
  if(!is.null(protec.areas)) {
    if(!class(protec.areas)=="SpatialPolygonsDataFrame") stop("protec.areas should be a spatialpolygondataframe")
    if(!any(colnames(protec.areas@data) %in% ID_shape_PA)) stop("Check argument ID_shape_PA because selected ID field in the protected area shapefile does not exist")
  }
  
    if(is.null(country_map)) {
      data('land', package='ConR', envir=environment()) 
      land <- get("land", envir=environment()) 
      country_map <- land
    }
  
  if(!is.null(protec.areas)) {
    if(!identicalCRS(protec.areas, land)) crs(protec.areas) <-crs(land)
  }
  
  if(length(grep("[?]", DATA[,3]))>0) DATA[,3] <- gsub("[?]", "_", DATA[,3])
  if(length(grep("[/]", DATA[,3]))>0) DATA[,3] <- gsub("[/]", "_", DATA[,3])
  
  #####
  
  list_data <- split(DATA, f = DATA$tax)

  if(map_pdf){
    if(!is.null(file_name)) {
      NAME_FILE <- file_name
    }else{
      NAME_FILE <- "IUCN_"
    }
    FILE_NAME <- ifelse(!is.null(file_name), file_name, "IUCN_")
    dir.create(file.path(paste(getwd(),paste("/",FILE_NAME,"_results_map", sep=""), sep="")), showWarnings = FALSE)
    
    pdf(paste(paste(getwd(),paste("/",FILE_NAME,"_results_map", sep=""), sep=""),"/","results.pdf", sep=""), width=25, height=25)
  }
  
  if(parallel) {
    if("doParallel" %in% 
       rownames(installed.packages()) == FALSE) {stop("Please install doParallel package")}
    
    library(doParallel)
    
    doParallel::registerDoParallel(NbeCores)
    message('doParallel running with ',
            NbeCores, ' cores')
    `%d%` <- foreach::`%dopar%`
  }else{
    `%d%` <- foreach::`%do%`
  }
  
  pb <- 
    utils::txtProgressBar(min = 0, max = length(list_data), style = 3)
  
  # if(parallel) doParallel::registerDoParallel(NbeCores)

  Results <- 
    foreach::foreach(x = 1:length(list_data), 
                     # .combine='c', 
                     arg1 = names(list_data), 
                     arg2 = rep(DrawMap, length(list_data)),
                     arg3 = rep(exclude.area, length(list_data)),
                     arg4 = rep(write_shp, length(list_data)),
                     arg6 = rep(method_protected_area, length(list_data)),
                     arg7 = rep(Cell_size_AOO, length(list_data)),
                     arg8 = rep(Cell_size_locations, length(list_data)),
                     arg9 = rep(Resol_sub_pop, length(list_data)),
                     arg10 = rep(method_locations, length(list_data)),
                     arg11 = rep(Rel_cell_size, length(list_data)),
                     arg12 = rep(ifelse(is.null(file_name), NA, file_name), length(list_data)),
                     arg13 = rep(buff_width, length(list_data)),
                     arg14 = rep(map_pdf, length(list_data)),
                     arg15 = rep(ID_shape_PA, length(list_data)),
                     arg16 = rep(SubPop, length(list_data)),
                     arg17 = rep(ifelse(is.null(protec.areas), NA, protec.areas), length(list_data)),
                     arg18 = rep(add.legend, length(list_data)),
                     arg20 = rep(alpha, length(list_data)),
                     arg21 = rep(buff.alpha, length(list_data)),
                     arg22 = rep(method.range, length(list_data)),
                     arg23 = rep(nbe.rep.rast.AOO, length(list_data)),
                     arg24 = rep(showWarnings, length(list_data))) %d% {
                       # source("./R/IUCNeval.functionv11.R")
                       utils::setTxtProgressBar(pb, x)
                       cat(" ", x)
                       
                       if(is.na(arg12)) arg_file_name <- NULL
                       if(is.na(arg17)) arg_protec.areas <- NULL
                       res <- 
                         .IUCN.comp(DATA = list_data[[x]],
                                    NamesSp = arg1, 
                                    DrawMap = arg2, 
                                    exclude.area = arg3, 
                                    write_shp = arg4, 
                                    method_protected_area = arg6, 
                                    Cell_size_AOO = arg7, 
                                    Cell_size_locations = arg8, 
                                    Resol_sub_pop = arg9, 
                                    method_locations = arg10, 
                                    Rel_cell_size = arg11, 
                                    poly_borders = country_map, 
                                    file_name = arg_file_name, 
                                    buff_width = arg13, 
                                    map_pdf = arg14, 
                                    ID_shape_PA = arg15, 
                                    SubPop = arg16, 
                                    protec.areas = arg_protec.areas, 
                                    add.legend = arg18,
                                    alpha = arg20, 
                                    buff.alpha = arg21, 
                                    method.range = arg22, 
                                    nbe.rep.rast.AOO = arg23, 
                                    showWarnings = arg24,
                                    MinMax = c(min(list_data[[x]][,2]), max(list_data[[x]][,2]), min(list_data[[x]][,1]), max(list_data[[x]][,1])))
                       
                       # names(res)[1] <- paste0(names(res)[1], "_" ,x)
                       # if(length(res)>1) names(res)[2] <- paste0(names(res)[2], "_" ,x)
                       
                       res
                     }
  
  # Results <- 
  #   plyr::llply(list_data, .fun=function(x) {
  #   .IUCN.comp(x, NamesSp=as.character(unique(x$tax)), 
  #              DrawMap=DrawMap, 
  #              exclude.area=exclude.area,
  #              write_shp=write_shp, 
  #              poly_borders=country_map, 
  #              method_protected_area=method_protected_area, 
  #              Cell_size_AOO=Cell_size_AOO, 
  #              Cell_size_locations=Cell_size_locations, 
  #              Resol_sub_pop=Resol_sub_pop,
  #              method_locations=method_locations, 
  #              Rel_cell_size=Rel_cell_size,
  #              file_name=file_name, 
  #              buff_width=buff_width, 
  #              map_pdf=map_pdf,
  #              ID_shape_PA=ID_shape_PA, 
  #              SubPop=SubPop,protec.areas=protec.areas, 
  #              MinMax=c(min(DATA[,2]), max(DATA[,2]), min(DATA[,1]), max(DATA[,1])),
  #              alpha=alpha, 
  #              buff.alpha=buff.alpha, 
  #              method.range=method.range, 
  #              nbe.rep.rast.AOO=nbe.rep.rast.AOO, #verbose=TRUE, verbose=verbose, 
  #              showWarnings=showWarnings, 
  #              draw.poly.EOO=draw.poly.EOO)
  #               }
  #               , .progress = "text", .parallel=parallel, .paropts = "ConR")
  
  
  if(parallel) doParallel::stopImplicitCluster()
  
  if(map_pdf) dev.off()
  
  Results_short <- lapply(Results, `[`, 1)
  Results_short <- lapply(Results_short, FUN = function(x) t(x[[1]]))
  Results_short <- as.data.frame(do.call(rbind, Results_short), stringsAsFactors=FALSE)
  Results_short[,1:4] <- apply(Results_short[,1:4], MARGIN = 2, as.numeric)
  
  if(write_results) {
      if(!is.null(file_name)) {
        NAME_FILE <- file_name
    }else{
        NAME_FILE <- "IUCN_results"
    }
    
    if(write_file_option=="csv") write.csv(Results_short, paste(getwd(),"/", NAME_FILE, ".csv", sep=""))

    if(write_file_option=="excel") {
      Results_short <- data.frame(taxa=rownames(Results_short), Results_short)
      write_xlsx(Results_short, path = paste(getwd(),"/", NAME_FILE, ".xlsx", sep=""))
    }
  
  }

  if(!export_shp) {
    Results <- Results_short
    if(length(list_data)>5) {
      print("Number of species per category")
      print(table(Results[,"Category_CriteriaB"]))
      print("Ratio of species per category")
      print(round(table(Results[,"Category_CriteriaB"])/nrow(Results)*100,1))
    }
  }
  
  Results
}

### Internal function used within map.res function
.prop_threat <- function(Cell_count, threshold) {
  NbeRec <- nrow(Cell_count)
  if(NbeRec >= threshold) {
    NbeEsp <- length(unique(Cell_count$tax))
    NbeThreatened <- length(unique(Cell_count[which(Cell_count$Category_CriteriaB %in% c("CR","EN","VU")),"tax"]))
    PropThreatened <- round(NbeThreatened/NbeEsp*100,1)
  }else{
    NbeEsp <- NbeThreatened <- PropThreatened <- NA
  }
  c(NbeRec, NbeEsp, NbeThreatened, PropThreatened)
}


map.res <- function(Results, 
                    Occurrences, 
                    country_map=NULL, 
                    Resol=1, 
                    threshold=0, 
                    LatMin=NULL, 
                    LatMax=NULL, 
                    LongMin=NULL, 
                    LongMax=NULL, 
                    export_map=FALSE, 
                    file_name=NULL, 
                    export_data=FALSE) {
  
  if(nrow(Results)!=length(unique(as.character(Occurrences[[3]])))) stop("Results and Occurrences input files have different number of species")
  
  if(class(Occurrences)=="spgeoIN") {
    DATA_2 <- cbind(Occurrences$species_coordinates, Occurrences$identifier)
    Occurrences <- DATA_2[,c(2,1,3)]
    colnames(Occurrences)[1:3] <- c("ddlat","ddlon","tax")
  }else{
    colnames(Occurrences)[1:3] <- c("ddlat","ddlon","tax")
  }
  
  Results_full <- cbind(rownames(Results), Results)
  colnames(Results_full)[1] <- "tax"
  merged_data_criteriaB <- merge(Results_full, Occurrences, by.x="tax", by.y="tax")
  
  if(is.null(LatMin)) LatMin = min(merged_data_criteriaB[,"ddlat"])
  if(is.null(LatMax)) LatMax = max(merged_data_criteriaB[,"ddlat"])
  if(is.null(LongMin)) LongMin = min(merged_data_criteriaB[,"ddlon"])
  if(is.null(LongMax)) LongMax = max(merged_data_criteriaB[,"ddlon"])
  
  if(LatMin>=LatMax) stop("LatMin must be lower than LatMax")
  if(LongMin>=LongMax) stop("LongMin must be lower than LongMax")
  if(LongMin>180 || LongMin< -180 || LatMin>180 || LatMin< -180 || LatMax>180 || LatMax< -180 || LongMax>180 || LongMax< -180) stop("Latitude and longitude must be within [-180; 180] intervall")
  
  EXTENT <- extent(LongMin, LongMax, LatMin, LatMax)
  
  if(is.null(country_map))  {
    data('land', package='ConR', envir=environment()) 
    land <- get("land", envir=environment()) 
    # data(land, envir = environment())
    country_map=land
  }
  
  if(!is.null(country_map)) {
    cropped_country_map <- crop(country_map, EXTENT+20)
  }
  
  DATA_SF <- merged_data_criteriaB
  coordinates(DATA_SF) <-  ~ddlon+ddlat
  crs(DATA_SF) <- crs(country_map)
  DATA_SF$X <- floor(merged_data_criteriaB[,"ddlon"]/Resol)
  DATA_SF$Y <- floor(merged_data_criteriaB[,"ddlat"]/Resol)
  DATA_SF$Cell <- paste("M",DATA_SF$X,"x",DATA_SF$Y, sep="")
  DATA_SF@data <- cbind(DATA_SF@data, merged_data_criteriaB[,c("ddlat","ddlon")])
  
  print(paste("Number of cell with at least one occurrence is", length(unique(as.character(DATA_SF@data[,"Cell"])))))
  print(paste("Number of cell with number of occurrences higher or equal to",threshold,"is", length(which(table(DATA_SF@data[,"Cell"])>threshold))))
  if(length(which(table(DATA_SF@data[,"Cell"])>threshold))==0) stop("No cell left after filtering")
  
  counts <- by(DATA_SF@data, DATA_SF@data$Cell, function(d) c(d$X[1], d$Y[1], mean(d$ddlat), mean(d$ddlon), 
                                                              .prop_threat(d, threshold)))
  
  threatened_rec <- matrix(unlist(counts), nrow=8)
  rownames(threatened_rec) <- c("X", "Y","meanLat","meanLong", "NbeRec", "NbeEsp", "NbeThreatened", "PropThreatened")
  colnames(threatened_rec) <- names(counts)
  if(ncol(threatened_rec)<2) stop("All records are within one grid cell, decrease the resolution to have a relevant map")
  
  threatened_rec_cut <- as.data.frame(threatened_rec[, which(threatened_rec["Y",]*Resol> (EXTENT+Resol)[3])])
  if(ncol(threatened_rec_cut)<2) stop("All records are within one grid cell, decrease the resolution/threshold or modify extent")
  threatened_rec_cut <- as.data.frame(threatened_rec_cut[,which(threatened_rec_cut["Y",]*Resol< (EXTENT+Resol)[4])])
  if(ncol(threatened_rec_cut)<2) stop("All records are within one grid cell, decrease the resolution/threshold or modify extent")
  threatened_rec_cut <- as.data.frame(threatened_rec_cut[,which(threatened_rec_cut["X",]*Resol> (EXTENT+Resol)[1])])
  if(ncol(threatened_rec_cut)<2) stop("All records are within one grid cell, decrease the resolution/threshold or modify extent")
  threatened_rec_cut <- as.data.frame(threatened_rec_cut[,which(threatened_rec_cut["X",]*Resol< (EXTENT+Resol)[2])])
  if(ncol(threatened_rec_cut)<2) stop("All records are within one grid cell, decrease the resolution/threshold or modify extent")
  
  COORD <- t(rbind(threatened_rec_cut[1,],threatened_rec_cut[2,]))
  
  grid.list <- list( x= (Resol/2+Resol*seq(range(COORD[,1])[1],range(COORD[,1])[2],1)), y=(Resol/2+Resol*seq(range(COORD[,2])[1],range(COORD[,2])[2],1)))
  
  
  SelectedCells <- which(threatened_rec_cut["NbeRec",]>=threshold)
  
  if(export_map) {
    FILE_NAME <- ifelse(!is.null(file_name), file_name, "IUCN_")
    dir.create(file.path(paste(getwd(),paste("/",FILE_NAME,"_results_map", sep=""), sep="")), showWarnings = FALSE)    
  }
  
  BG <- 'white'
  Border <- 'black'
  COlor <- 'grey97' # rgb(0.1, 0.3, 0.1, alpha=0.1)
  
  if(dev.cur()>1) {dev.off()}
  if(dev.cur()>1) {dev.off()}
  
  if (!export_map) par(mfrow=c(2,2))
  if (export_map) par(mfrow=c(1,1))
  if (export_map) png(paste(file.path(paste(getwd(),paste("/",FILE_NAME,"_results_map", sep=""), sep="")),"/","number_of_records",".png", sep=""),width=20, height=20, units="cm",res=150)
  
  ################################### Number of records
  coltab<- two.colors(256, start="lightblue", end="red", middle="yellow")
  
  if (!export_map) par(mar=c(4,1,1,4), las=1, omi=c(0.5,1,0.5,0.3))
  if (export_map) par(mar=c(2,2,1,5), las=1, omi=c(0.3,0.4,0.3,0.1))
  
  if (export_map) plot(cropped_country_map, axes=T, lty=1,border=Border, col=COlor, xlim=c(LongMin,LongMax), ylim=c(LatMin,LatMax), cex.axis=1,lwd=1)

  if (!export_map) plot(cropped_country_map, axes=T, lty=1,border=Border, col=COlor, xlim=c(LongMin,LongMax), ylim=c(LatMin,LatMax), cex.axis=1,lwd=1, xaxt='n')
  
  mtext(text="Number of records",side=3, cex=1)
  VALUES <- as.numeric(threatened_rec_cut["NbeRec",SelectedCells])
  
  quilt.plot(COORD[SelectedCells,1]*Resol+Resol/2 , COORD[SelectedCells,2]*Resol+Resol/2, VALUES  ,grid=grid.list , cex.axis=1, 
             cex.lab=1, add.legend=FALSE, col=coltab, add=T)
  plot(cropped_country_map, add=T)
  if(min(VALUES)==max(VALUES)) Range <- c(min(VALUES), min(VALUES)+1)
  if(min(VALUES)!=max(VALUES)) Range <- range(VALUES)
  image.plot(zlim=Range,legend.only=TRUE, col=coltab, legend.shrink = 1 ,
             legend.width=1, cex.lab=2, axis.args=list(cex.axis = 1, col.lab = Border, col.axis = Border))
  box()
  if (export_map) dev.off()
  
  ################################### Number of species
  coltab<- two.colors(256, start="cyan", end="darkorange4", middle="gold")
  
  if (export_map) png(paste(file.path(paste(getwd(),paste("/",FILE_NAME,"_results_map", sep=""), sep="")),"/","species_richness",".png", sep=""),width=20, height=20, units="cm",res=150)
  if (export_map)   par(mar=c(2,2,1,5), las=1, omi=c(0.3,0.4,0.3,0.1))
  if (export_map) plot(cropped_country_map, axes=T, lty=1,border=Border, col=COlor, xlim=c(LongMin,LongMax), ylim=c(LatMin,LatMax), cex.axis=1,lwd=1)
  
  if (!export_map)   plot(cropped_country_map, axes=T, lty=1,border=Border, col=COlor, xlim=c(LongMin,LongMax), ylim=c(LatMin,LatMax), cex.axis=1,lwd=1, xaxt='n', yaxt='n')
  mtext(text="Species richness",side=3, cex=1)
  VALUES <- as.numeric(threatened_rec_cut["NbeEsp",SelectedCells])
  
  quilt.plot(COORD[SelectedCells,1]*Resol+Resol/2 , COORD[SelectedCells,2]*Resol+Resol/2, VALUES  ,grid=grid.list , cex.axis=1, 
             cex.lab=1, add.legend=FALSE, col=coltab, add=T)
  plot(cropped_country_map, add=T)
  if(min(VALUES)==max(VALUES)) Range <- c(min(VALUES), min(VALUES)+1)
  if(min(VALUES)!=max(VALUES)) Range <- range(VALUES)
  image.plot(zlim=Range,legend.only=TRUE, col=coltab, legend.shrink = 1 ,
             legend.width=1, cex.lab=2, axis.args=list(cex.axis = 1, col.lab = Border, col.axis = Border))
  box()
  if (export_map) dev.off()

  ################################### Number of threatened species
  coltab<- two.colors(256, start="slategray2", end="deeppink4", middle="burlywood2")
  
  if (export_map) png(paste(file.path(paste(getwd(),paste("/",FILE_NAME,"_results_map", sep=""), sep="")),"/","number_threatened_sp",".png", sep=""),width=20, height=20, units="cm",res=150)
  if (export_map)   par(mar=c(2,2,1,5), las=1, omi=c(0.3,0.4,0.3,0.1))

  if (export_map) plot(cropped_country_map, axes=T, lty=1,border=Border, col=COlor, xlim=c(LongMin,LongMax), ylim=c(LatMin,LatMax), cex.axis=1,lwd=1)

  if (!export_map)  plot(cropped_country_map, axes=T, lty=1,border=Border, col=COlor, xlim=c(LongMin,LongMax), ylim=c(LatMin,LatMax), cex.axis=1,lwd=1)
  if (export_map)   par(mar=c(2,2,1,5), las=1, omi=c(0.3,0.4,0.3,0.1))
  mtext(text="Number of threatened species",side=3, cex=1)
  VALUES <- as.numeric(threatened_rec_cut["NbeThreatened",SelectedCells])
  
  quilt.plot(COORD[SelectedCells,1]*Resol+Resol/2 , COORD[SelectedCells,2]*Resol+Resol/2, VALUES  ,grid=grid.list , cex.axis=1, 
             cex.lab=1, add.legend=FALSE, col=coltab, add=T)
  plot(cropped_country_map, add=T)
  if(min(VALUES)==max(VALUES)) Range <- c(min(VALUES), min(VALUES)+1)
  if(min(VALUES)!=max(VALUES)) Range <- range(VALUES)
  image.plot(zlim=Range,legend.only=TRUE, col=coltab, legend.shrink = 1 ,
             legend.width=1, cex.lab=2, axis.args=list(cex.axis = 1, col.lab = Border, col.axis = Border))
  box()
  if (export_map) dev.off()

  ################################### Proportion of threatened species
  coltab<- two.colors(256, start="darkslategray1", end="hotpink2", middle="khaki")
  
  if (export_map) png(paste(file.path(paste(getwd(),paste("/",FILE_NAME,"_results_map", sep=""), sep="")),"/","proportion_threatened_sp",".png", sep=""),width=20, height=20, units="cm",res=150)
  if (export_map) par(mar=c(2,2,1,5), las=1, omi=c(0.3,0.4,0.3,0.1))

  if (!export_map) plot(cropped_country_map, axes=T, lty=1,border=Border, col=COlor, xlim=c(LongMin,LongMax), ylim=c(LatMin,LatMax), cex.axis=1,lwd=1, yaxt='n')
  if (export_map) plot(cropped_country_map, axes=T, lty=1,border=Border, col=COlor, xlim=c(LongMin,LongMax), ylim=c(LatMin,LatMax), cex.axis=1,lwd=1)

  mtext(text="Proportion of threatened species",side=3, cex=1)
  VALUES <- as.numeric(threatened_rec_cut["PropThreatened",SelectedCells])
  
  quilt.plot(COORD[SelectedCells,1]*Resol+Resol/2 , COORD[SelectedCells,2]*Resol+Resol/2, VALUES  ,grid=grid.list , cex.axis=1, 
             cex.lab=1, add.legend=FALSE, col=coltab, add=T)
  plot(cropped_country_map, add=T)
  if(min(VALUES)==max(VALUES)) Range <- c(min(VALUES), min(VALUES)+1)
  if(min(VALUES)!=max(VALUES)) Range <- range(VALUES)
  image.plot(zlim=Range,legend.only=TRUE, col=coltab, legend.shrink = 1 ,
             legend.width=1, cex.lab=2, axis.args=list(cex.axis = 1, col.lab = Border, col.axis = Border))
  box()
  if (export_map) dev.off()
  
##########################################
  if(export_data) return(threatened_rec_cut)
  
}


