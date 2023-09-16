#' @title Internal function
#'
#' @description Count number of occupied cells given resolution, projection
#' 
#' @param nbe_rep integer
#' @param size integer
#' @param coord data.frame
#' @param export_shp logical
#' @param proj_type character string
#' 
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#'
#' @importFrom terra ext rast vect rasterize values as.polygons project
#' @importFrom stats runif
#' @import sf
cell.occupied <-
  function(nbe_rep = 0,
           size = 4,
           coord,
           export_shp = TRUE,
           proj_type = NULL)  {

    crs_proj <-
      proj_type

    Corners <- rbind(c(min(coord[, 1]),
                       max(coord[, 1])),
                     c(min(coord[, 2]),
                       max(coord[, 2])))


    if (nbe_rep == 0) {

      Occupied_cells <- vector(mode = "numeric", length = 4)
      rasts <- vector(mode = "list", length = 4)
      decal <- c(0, 1, 2, 3)

      
      for (h in decal) {

        ext <-
          terra::ext(
            floor(Corners[1, 1]) - h * (size * 1000 / 4) - 2 * size * 1000,
            floor(Corners[1, 2]) + h * (size * 1000 / 4) + 2 *
              size * 1000,
            floor(Corners[2, 1]) - h * (size * 1000 / 4) - 2 *
              size * 1000,
            floor(Corners[2, 2]) + h * (size * 1000 / 4) + 2 *
              size * 1000
          )
        
        r <-
          terra::rast(ext,
                      resolution = size * 1000,
                      crs = crs_proj$wkt)
        
        r_proj <- project(rast(r), "epsg:4326")
        
        # ext <-
        #   raster::extent(
        #     floor(Corners[1, 1]) - h * (size * 1000 / 4) - 2 * size * 1000,
        #     floor(Corners[1, 2]) + h * (size * 1000 / 4) + 2 *
        #       size * 1000,
        #     floor(Corners[2, 1]) - h * (size * 1000 / 4) - 2 *
        #       size * 1000,
        #     floor(Corners[2, 2]) + h * (size * 1000 / 4) + 2 *
        #       size * 1000
        #   )
        # 
        # r <-
        #   raster::raster(ext,
        #                  resolution = size * 1000,
        #                  crs = crs_proj)
        
        coord_vec <- coord[, 1:2]
        colnames(coord_vec) <- c("x", "y")
        coord_vec <- terra::vect(coord_vec, 
                                 crs = crs_proj$wkt, geom=c("x", "y"))
        coord_vec_proj <- 
          project(coord_vec, "epsg:4326")
        
        r2_ <-
          terra::rasterize(x = coord_vec_proj, y = r_proj)
        
        OCC <-
          length(which(!is.na(terra::values(r2_))))
        
        rasts[[h + 1]] <- 
          r2_
        
        Occupied_cells[h + 1] <-
          OCC
        
        ### If only one occupied cell, stop the production of raster
        if (OCC == 1)
          break

        # r2_ <-
        #   raster::rasterize(coord[, 1:2], r)
        # 
        # rasts[[h + 1]] <- 
        #   r2_
        # 
        # OCC <-
        #   length(which(!is.na(raster::values(r2_))))
        # 
        # Occupied_cells[h + 1] <- OCC
        # 
        # ### If only one occupied cell, stop the production of raster
        # if (OCC == 1)
        #   break
      }
      # h <- decal[which.min(Occupied_cells)]
      # Occupied_cells <- min(Occupied_cells)
    }

    if (nbe_rep > 0) {
      Occupied_cells <- vector(mode = "numeric", length = nbe_rep)
      rasts <- vector(mode = "list", length = 4)
      
      for (h in 1:nbe_rep) {
        rd.1 <- runif(1) * size * 1000
        rd.2 <- runif(1) * size * 1000

        ext <-
          terra::ext(
            floor(Corners[1, 1]) - rd.1 - 2 * size * 1000,
            floor(Corners[1, 2]) + rd.1 + 2 * size * 1000,
            floor(Corners[2, 1]) - rd.2 - 2 * size * 1000,
            floor(Corners[2, 2]) + rd.2 + 2 * size * 1000
          )
        
        r <-
          terra::rast(ext,
                      resolution = size * 1000,
                      crs = crs_proj$wkt)
        
        r_proj <- project(rast(r), "epsg:4326")
        
        # ext = raster::extent(
        #   floor(Corners[1, 1]) - rd.1 - 2 * size * 1000,
        #   floor(Corners[1, 2]) + rd.1 + 2 * size * 1000,
        #   floor(Corners[2, 1]) - rd.2 - 2 * size * 1000,
        #   floor(Corners[2, 2]) + rd.2 + 2 * size * 1000
        # )
        # r = raster::raster(ext, resolution = size * 1000, crs = crs_proj)
        # r
        
        coord_vec <- coord[, 1:2]
        colnames(coord_vec) <- c("x", "y")
        coord_vec <- terra::vect(coord_vec, 
                                 crs = crs_proj$wkt, geom=c("x", "y"))
        coord_vec_proj <- 
          project(coord_vec, "epsg:4326")
        
        r2_ <-
          terra::rasterize(x = coord_vec_proj, y = r_proj)
        
        # coord_vec <- coord[, 1:2]
        # colnames(coord_vec) <- c("x", "y")
        # coord_vec <- terra::vect(coord_vec, 
        #                          crs = as.character(crs_proj), geom=c("x", "y"))
        # r2_ <-
        #   terra::rasterize(x = coord_vec, y = r)
        
        OCC <-
          length(which(!is.na(terra::values(r2_))))
        
        rasts[[h]] <- 
          r2_
        
        Occupied_cells[h + 1] <-
          OCC
        
        ### If only one occupied cell, stop the production of raster
        if (OCC == 1)
          break
        
        
        # r2_ <- raster::rasterize(coord[, 1:2], r)
        # OCC <- length(which(!is.na(raster::values(r2_))))
        # Occupied_cells[h] <- OCC
        # 
        # rasts[[h]] <- 
        #   r2_
        # 
        # # rd.1.vec <- c(rd.1.vec, rd.1)
        # # rd.2.vec <- c(rd.2.vec, rd.2)
        # if (OCC == 1)
        #   break
      }

    }

    which_raster <- which.min(Occupied_cells[Occupied_cells>0])
    Occupied_cells <- Occupied_cells[Occupied_cells > 0]
    Occupied_cells <- min(Occupied_cells)
    

    if (export_shp) {
      
      # r2_ <-
      #   suppressWarnings(terra::project(x = rasts[[which_raster]],
      #                                   y = "epsg:4326", 
      #                                   gdal = FALSE))
      
      r2_ <- 
        rasts[[which_raster]]
      
      # 
      
      r2_pol <-
        terra::as.polygons(
          r2_)
      
      r2_pol_sf <- sf::st_as_sf(r2_pol)
      
      r2_pol_sf <- suppressWarnings(sf::st_cast(r2_pol_sf, "POLYGON"))
      row.names(r2_pol_sf) <- 1:nrow(r2_pol_sf)
      
      # st_as_sf(stars::st_as_stars(r2_), merge = TRUE)
      # 
      # microbenchmark::microbenchmark({r2_pol_sf <- sf::st_as_sf(r2_pol); sf::st_as_sf(r2_pol);  r2_pol_sf <- sf::st_cast(r2_pol_sf, "POLYGON")},
      #                                st_as_sf(stars::st_as_stars(r2_), merge = TRUE), times = 10)
      
    } else {
      r2_pol_sf <- NA
    }
    
    # if (export_shp) {
    #   r2_pol <-
    #     raster::rasterToPolygons(
    #       rasts[[which_raster]],
    #       fun = NULL,
    #       n = 4,
    #       na.rm = TRUE,
    #       digits = 6,
    #       dissolve = FALSE
    #     )
    #   
    #   r2_pol <-
    #     as(r2_pol, "sf")
    # }

    return(list(r2_pol_sf, Occupied_cells))
    
  }


# tm <- microbenchmark::microbenchmark(cell.occupied.old(nbe_rep = 10, size = 4, coord = coord),
#                      cell.occupied(nbe_rep = 10, size = 4, coord = coord, proj_type = proj_type), 
#                      times=10L)


### Not working for points and small polygons

# cell.occupied.fasterize <-
#   function(nbe_rep = 0,
#            size = 4,
#            coord,
#            export_shp = TRUE,
#            proj_type = NULL)  {
#     
#     crs_proj <-
#       proj_type
#     
#     Corners <- rbind(c(min(coord[, 1]),
#                        max(coord[, 1])),
#                      c(min(coord[, 2]),
#                        max(coord[, 2])))
#     
#     xy_sf <- st_as_sf(coord[, 1:2], coords = c(1, 2), crs = crs_proj)
#     
#     xy_sf <- st_buffer(xy_sf, 1)
#     
#     xy_sf <- st_sf(value = rep(1,nrow(xy_sf)),
#           geometry = xy_sf$geometry)
#     
#     
#     if (nbe_rep == 0) {
#       
#       Occupied_cells <- vector(mode = "numeric", length = 4)
#       decal <- c(0, 1, 2, 3)
#       
#       # if (abs(Corners[1, 1] - Corners[1, 2]) > abs(Corners[2, 1] - Corners[2, 2])) {
#       #   longer <-
#       #     abs(Corners[1, 1] - Corners[1, 2])
#       # } else{
#       #   longer <-
#       #     abs(Corners[2, 1] - Corners[2, 2])
#       # }
#       
#       
#       for (h in decal) {
#         
#         # bbox_sp <-
#         #   st_bbox(c(xmin = min(coord[, 1]), xmax =
#         #               max(coord[, 1]), ymax = max(coord[, 2]), ymin = min(coord[, 2])),
#         #           crs = proj_type)
#         #
#         # grid_sp <-
#         #   st_make_grid(st_as_sfc(bbox_sp), cellsize = cell_size*1000)
#         
#         
#         # test <- st_sfc(st_linestring(rbind(c(0, 2), c(0 , 0))), crs = 4326)
#         #
#         # st_length(test)
#         
#         ext <-
#           raster::extent(
#             floor(Corners[1, 1]) - h * (size * 1000 / 4) - 2 * size * 1000,
#             floor(Corners[1, 2]) + h * (size * 1000 / 4) + 2 *
#               size * 1000,
#             floor(Corners[2, 1]) - h * (size * 1000 / 4) - 2 *
#               size * 1000,
#             floor(Corners[2, 2]) + h * (size * 1000 / 4) + 2 *
#               size * 1000
#           )
#         
#         
#         # ext <-
#         #   terra::ext(
#         #     c(
#         #       floor(Corners[1, 1]) - h * (size * 1000 / 4) - 2 *
#         #         size * 1000,
#         #       floor(Corners[1, 2]) + h * (size * 1000 / 4) + 2 *
#         #         size * 1000,
#         #       floor(Corners[2, 1]) - h * (size * 1000 / 4) - 2 *
#         #         size * 1000,
#         #       floor(Corners[2, 2]) + h * (size * 1000 / 4) + 2 *
#         #         size * 1000
#         #     )
#         #   )
#         #
#         # r <-
#         #   terra::rast(extent = ext,
#         #                  resolution = size * 1000,
#         #                  crs = crs_proj)
#         
#         
#         r <-
#           raster::raster(ext,
#                          resolution = size * 1000,
#                          crs = crs_proj)
#         
#         r2_ <-
#           raster::rasterize(coord[, 1:2], r)
#         
#         r2_ <-
#           fasterize::fasterize(xy_sf, r)
#         
#         fasterize(xy_sf, r, field = "value", fun="count")
#         
#         values(r2_)[!is.na(values(r2_))]
#         
#         OCC <-
#           length(which(!is.na(raster::values(r2_))))
#         
#         Occupied_cells[h + 1] <- OCC
#         
#         ### If only one occupied cell, stop the production of raster
#         if (OCC == 1)
#           break
#       }
#       # h <- decal[which.min(Occupied_cells)]
#       # Occupied_cells <- min(Occupied_cells)
#     }
#     
#     if (nbe_rep > 0) {
#       Occupied_cells <- vector(mode = "numeric", length = nbe_rep)
#       
#       for (h in 1:nbe_rep) {
#         rd.1 <- runif(1) * size * 1000
#         rd.2 <- runif(1) * size * 1000
#         
#         ext = raster::extent(
#           floor(Corners[1, 1]) - rd.1 - 2 * size * 1000,
#           floor(Corners[1, 2]) + rd.1 + 2 * size * 1000,
#           floor(Corners[2, 1]) - rd.2 - 2 * size * 1000,
#           floor(Corners[2, 2]) + rd.2 + 2 * size * 1000
#         )
#         r = raster::raster(ext, resolution = size * 1000, crs = crs_proj)
#         # r
#         r2_ <- raster::rasterize(coord[, 1:2], r)
#         OCC <- length(which(!is.na(raster::values(r2_))))
#         Occupied_cells[h] <- OCC
#         # rd.1.vec <- c(rd.1.vec, rd.1)
#         # rd.2.vec <- c(rd.2.vec, rd.2)
#         if (OCC == 1)
#           break
#       }
#       
#     }
#     
#     Occupied_cells <- Occupied_cells[Occupied_cells > 0]
#     Occupied_cells <- min(Occupied_cells)
#     
#     ## CRS object has comment, which is lost in output warning
#     # if (export_shp)
#     #   r2_ <-
#     #   suppressWarnings(raster::projectRaster(from = r2_,
#     #                                          crs = "+proj=longlat +datum=WGS84 +no_defs"))
#     
#     if (export_shp) {
#       r2_pol <-
#         raster::rasterToPolygons(
#           r2_,
#           fun = NULL,
#           n = 4,
#           na.rm = TRUE,
#           digits = 6,
#           dissolve = FALSE
#         )
#       
#       r2_pol <-
#         as(r2_pol, "sf")      
#     }
#     
#     
#     if (export_shp)
#       return(list(r2_pol, Occupied_cells))
#     if (!export_shp)
#       return(list(NA, Occupied_cells))
#     
#   }


# cell.occupied <-
#   function(nbe_rep = 0,
#            size = 4,
#            coord,
#            export_shp = TRUE,
#            proj_type = NULL) {
#     
#     crs_proj <-
#       proj_type
#     
#     Corners <- rbind(c(min(coord[, 1]),
#                        max(coord[, 1])),
#                      c(min(coord[, 2]),
#                        max(coord[, 2])))
#     
#     if (nbe_rep == 0) {
#       
#       Occupied_cells <- vector(mode = "numeric", length = 4)
#       decal <- c(0, 1, 2, 3)
#       
#       for (h in decal) {
#         
#         bbox_sp <-
#           st_bbox(c(xmin = floor(Corners[1, 1]) - h * (size * 1000 / 4) - 2 * size * 1000, 
#                     xmax = floor(Corners[1, 2]) + h * (size * 1000 / 4) + 2 * size * 1000, 
#                     ymax = floor(Corners[2, 2]) + h * (size * 1000 / 4) + 2 * size * 1000, 
#                     ymin = floor(Corners[2, 1]) - h * (size * 1000 / 4) - 2 * size * 1000),
#                   crs = proj_type)
#         
#         # ext <-
#         #   raster::extent(
#         #     floor(Corners[1, 1]) - h * (size * 1000 / 4) - 2 * size * 1000,
#         #     floor(Corners[1, 2]) + h * (size * 1000 / 4) + 2 *
#         #       size * 1000,
#         #     floor(Corners[2, 1]) - h * (size * 1000 / 4) - 2 *
#         #       size * 1000,
#         #     floor(Corners[2, 2]) + h * (size * 1000 / 4) + 2 *
#         #       size * 1000
#         #   )
#         
#         
#         r <-
#           stars::st_as_stars(
#             bbox_sp,
#             dx = size * 1000,
#             dy = size * 1000,
#             values = NA,
#             name = "aoo",
#             ignore_file = TRUE
#           )
#         
#         xy_sf <- 
#           st_as_sf(coord[, 1:2], coords = c(1, 2), crs = crs_proj)
#         
#         r <- stars::st_rasterize(xy_sf, r)
#         
#         # r_rast <-
#         #   raster::raster(ext,
#         #                  resolution = size * 1000,
#         #                  crs = crs_proj)
#         
#         # r2_ <-
#         #   raster::rasterize(coord[, 1:2], r_rast)
#         
#         # OCC <-
#         #   length(which(!is.na(raster::values(r2_))))
#         
#         OCC <-
#           length(r$ID[!is.na(r$ID)])
#         
#         rm(r)
#         
#         Occupied_cells[h + 1] <- OCC
#         
#         ### If only one occupied cell, stop the production of raster
#         if (OCC == 1)
#           break
#       }
#       h <- decal[which.min(Occupied_cells)]
#       
#       bbox_sp <-
#         st_bbox(c(xmin = floor(Corners[1, 1]) - h * (size * 1000 / 4) - 2 * size * 1000, 
#                   xmax = floor(Corners[1, 2]) + h * (size * 1000 / 4) + 2 * size * 1000, 
#                   ymax = floor(Corners[2, 2]) + h * (size * 1000 / 4) + 2 * size * 1000, 
#                   ymin = floor(Corners[2, 1]) - h * (size * 1000 / 4) - 2 * size * 1000),
#                 crs = proj_type)
#       
#       r <-
#         stars::st_as_stars(
#           bbox_sp,
#           dx = size * 1000,
#           dy = size * 1000,
#           values = NA,
#           name = "aoo",
#           ignore_file = TRUE
#         )
#       
#       xy_sf <- 
#         st_as_sf(coord[, 1:2], coords = c(1, 2), crs = crs_proj)
#       
#       r <- stars::st_rasterize(xy_sf, r, file = "temp_rast")
#       
#       # Occupied_cells <- min(Occupied_cells)
#     }
#     
#     if (nbe_rep > 0) {
#       
#       Occupied_cells <- vector(mode = "numeric", length = nbe_rep)
#       rd.1.v <- vector(mode = "numeric", length = nbe_rep)
#       rd.2.v <- vector(mode = "numeric", length = nbe_rep)
#       
#       for (h in 1:nbe_rep) {
#         rd.1 <- runif(1) * size * 1000
#         rd.2 <- runif(1) * size * 1000
#         
#         bbox_sp <-
#           st_bbox(c(xmin = floor(Corners[1, 1]) - rd.1 - 2 * (size * 1000 / 4) - 2 * size * 1000, 
#                     xmax = floor(Corners[1, 2]) + rd.1 + 2 * (size * 1000 / 4) + 2 * size * 1000, 
#                     ymax = floor(Corners[2, 2]) + rd.2 + 2 * (size * 1000 / 4) + 2 * size * 1000, 
#                     ymin = floor(Corners[2, 1]) - rd.2 - 2 * (size * 1000 / 4) - 2 * size * 1000),
#                   crs = proj_type)
#         
#         # ext = raster::extent(
#         #   floor(Corners[1, 1]) - rd.1 - 2 * size * 1000,
#         #   floor(Corners[1, 2]) + rd.1 + 2 * size * 1000,
#         #   floor(Corners[2, 1]) - rd.2 - 2 * size * 1000,
#         #   floor(Corners[2, 2]) + rd.2 + 2 * size * 1000
#         # )
#         # r = raster::raster(ext, resolution = size * 1000, crs = crs_proj)
#         # # r
#         # r2_ <- raster::rasterize(coord[, 1:2], r)
#         # OCC <- length(which(!is.na(raster::values(r2_))))
#         # Occupied_cells[h] <- OCC
#         
#         r <-
#           stars::st_as_stars(
#             bbox_sp,
#             dx = size * 1000,
#             dy = size * 1000,
#             values = NA,
#             name = "aoo",
#             ignore_file = TRUE
#           )
#         
#         xy_sf <- 
#           st_as_sf(coord[, 1:2], coords = c(1, 2), crs = crs_proj)
#         
#         r <- stars::st_rasterize(xy_sf, r, file = "temp_rast2")
#         
#         OCC <-
#           length(r$ID[!is.na(r$ID)])
#         
#         rm(r)
#         
#         Occupied_cells[h] <- OCC
#         rd.1.v[h] <- rd.1
#         rd.2.v[h] <- rd.2
#         
#         # rd.1.vec <- c(rd.1.vec, rd.1)
#         # rd.2.vec <- c(rd.2.vec, rd.2)
#         if (OCC == 1)
#           break
#       }
#       
#       rd.1 <- rd.1.v[which.min(Occupied_cells)]
#       rd.2 <- rd.2.v[which.min(Occupied_cells)]
#       
#       bbox_sp <-
#         st_bbox(c(xmin = floor(Corners[1, 1]) - rd.1 - 2 * (size * 1000 / 4) - 2 * size * 1000, 
#                   xmax = floor(Corners[1, 2]) + rd.1 + 2 * (size * 1000 / 4) + 2 * size * 1000, 
#                   ymax = floor(Corners[2, 2]) + rd.2 + 2 * (size * 1000 / 4) + 2 * size * 1000, 
#                   ymin = floor(Corners[2, 1]) - rd.2 - 2 * (size * 1000 / 4) - 2 * size * 1000),
#                 crs = proj_type)
#       
#       r <-
#         stars::st_as_stars(
#           bbox_sp,
#           dx = size * 1000,
#           dy = size * 1000,
#           values = NA,
#           name = "aoo",
#           ignore_file = TRUE
#         )
#       
#       xy_sf <- 
#         st_as_sf(coord[, 1:2], coords = c(1, 2), crs = crs_proj)
#       
#       r <- stars::st_rasterize(xy_sf, r)
#       
#     }
#     
#     Occupied_cells <- Occupied_cells[Occupied_cells > 0]
#     Occupied_cells <- min(Occupied_cells)
#     
#     ## CRS object has comment, which is lost in output warning
#     # if (export_shp)
#     #   r2_ <-
#     #   suppressWarnings(raster::projectRaster(from = r2_,
#     #                                          crs = "+proj=longlat +datum=WGS84 +no_defs"))
#     
#     
#     r$ID[!is.na(r$ID)] <- 1
#     
#     if (export_shp)
#       r2_pol <- 
#       suppressWarnings(st_as_sf(r, na.rm = TRUE, merge = FALSE))
#     
#     # rm(r)
#     # 
#     # removeTmpFiles()
#     
#     if (export_shp)
#       return(list(r2_pol, Occupied_cells))
#     if (!export_shp)
#       return(list(NA, Occupied_cells))
#   }


# .cell.occupied.terra <-
#   function(nbe_rep = 0,
#            size = 4,
#            coord,
#            export_shp = TRUE,
#            proj_type = NULL) {
# 
#     crs_proj <-
#       proj_crs()
# 
#     Corners <- rbind(c(min(coord[, 1]),
#                        max(coord[, 1])),
#                      c(min(coord[, 2]),
#                        max(coord[, 2])))
# 
#     if (nbe_rep == 0) {
# 
#       Occupied_cells <- vector(mode = "numeric", length = 4)
#       decal <- c(0, 1, 2, 3)
#       rasts <- vector(mode = "list", length = 4)
# 
#       if (abs(Corners[1, 1] - Corners[1, 2]) > abs(Corners[2, 1] - Corners[2, 2])) {
#         longer <-
#           abs(Corners[1, 1] - Corners[1, 2])
#       } else{
#         longer <-
#           abs(Corners[2, 1] - Corners[2, 2])
#       }
# 
# 
#       for (h in decal) {
# 
#         ext <-
#           terra::ext(
#             floor(Corners[1, 1]) - h * (size * 1000 / 4) - 2 * size * 1000,
#             floor(Corners[1, 2]) + h * (size * 1000 / 4) + 2 *
#               size * 1000,
#             floor(Corners[2, 1]) - h * (size * 1000 / 4) - 2 *
#               size * 1000,
#             floor(Corners[2, 2]) + h * (size * 1000 / 4) + 2 *
#               size * 1000
#           )
# 
#         r <-
#           terra::rast(ext,
#                       resolution = size * 1000,
#                       crs = as.character(crs_proj))
# 
# 
#         # ext_rast <-
#         #   raster::extent(
#         #     floor(Corners[1, 1]) - h * (size * 1000 / 4) - 2 * size * 1000,
#         #     floor(Corners[1, 2]) + h * (size * 1000 / 4) + 2 *
#         #       size * 1000,
#         #     floor(Corners[2, 1]) - h * (size * 1000 / 4) - 2 *
#         #       size * 1000,
#         #     floor(Corners[2, 2]) + h * (size * 1000 / 4) + 2 *
#         #       size * 1000
#         #   )
#         # 
#         # r <-
#         #   raster::raster(ext_rast,
#         #                  resolution = size * 1000,
#         #                  crs = crs_proj)
# 
# 
#         # terra::values(r) <- NA
# 
#         coord_vec <- coord[, 1:2]
#         colnames(coord_vec) <- c("x", "y")
#         coord_vec <- terra::vect(coord_vec, 
#                                  crs = as.character(crs_proj), geom=c("x", "y"))
#         r2_ <-
#           terra::rasterize(x = coord_vec, y = r)
# 
#         OCC <-
#           length(which(!is.na(terra::values(r2_))))
#         
#         rasts[[h + 1]] <- 
#           r2_
# 
#         Occupied_cells[h + 1] <-
#           OCC
# 
#         ### If only one occupied cell, stop the production of raster
#         if (OCC == 1)
#           break
#       }
#       # h <- decal[which.min(Occupied_cells)]
#       # Occupied_cells <- min(Occupied_cells)
#     }
# 
#     if (nbe_rep > 0) {
# 
#       Occupied_cells <- vector(mode = "numeric", length = nbe_rep)
# 
#       for (h in 1:nbe_rep) {
#         rd.1 <- runif(1) * size * 1000
#         rd.2 <- runif(1) * size * 1000
# 
#         ext <-
#           terra::ext(
#             floor(Corners[1, 1]) - rd.1 - 2 * size * 1000,
#             floor(Corners[1, 2]) + rd.1 + 2 * size * 1000,
#             floor(Corners[2, 1]) - rd.2 - 2 * size * 1000,
#             floor(Corners[2, 2]) + rd.2 + 2 * size * 1000
#           )
# 
#         # ext = raster::extent(
#         #   floor(Corners[1, 1]) - rd.1 - 2 * size * 1000,
#         #   floor(Corners[1, 2]) + rd.1 + 2 * size * 1000,
#         #   floor(Corners[2, 1]) - rd.2 - 2 * size * 1000,
#         #   floor(Corners[2, 2]) + rd.2 + 2 * size * 1000
#         # )
# 
#         coord_vec <- coord[, 1:2]
#         colnames(coord_vec) <- c("x", "y")
#         coord_vec <- terra::vect(coord_vec, crs = as.character(crs_proj))
#         r2_ <-
#           terra::rasterize(x = coord_vec, y = r)
# 
#         # r <- terra::rast(ext, resolution = size * 1000, crs = crs_proj)
#         # # r
#         # r2_ <- terra::rasterize(coord[, 1:2], r)
#         OCC <- length(which(!is.na(terra::values(r2_))))
#         Occupied_cells[h] <- OCC
#         # rd.1.vec <- c(rd.1.vec, rd.1)
#         # rd.2.vec <- c(rd.2.vec, rd.2)
#         if (OCC == 1)
#           break
#       }
# 
#     }
# 
#     Occupied_cells <- Occupied_cells[Occupied_cells > 0]
#     Occupied_cells <- min(Occupied_cells)
#     which_raster <- which.min(Occupied_cells[Occupied_cells>0])
#     
#     if (export_shp) {
#       
#       r2_ <-
#         suppressWarnings(terra::project(x = rasts[[which_raster]],
#                                         y = "epsg:4326"))
#       
#       r2_pol <-
#         terra::as.polygons(
#           r2_)
#       
#       r2_pol_sf <- sf::st_as_sf(r2_pol)
#       
#     }
# 
#     ## CRS object has comment, which is lost in output warning
#     # if (export_shp)
#     #   r2_ <-
#     #   suppressWarnings(terra::project(y = r2_,
#     #                                   crs = "+proj=longlat +datum=WGS84 +no_defs"))
#     # 
#     # 
#     # r2_proj <-
#     #   terra::rast(ext, crs = as.character("+proj=longlat +datum=WGS84 +no_defs"))
# 
#     # if (export_shp)
#     #   r2_pol <-
#     #   terra::as.polygons(
#     #     r2_
#     #   )
# 
#     if (export_shp)
#       return(list(r2_pol_sf, Occupied_cells))
#     if (!export_shp)
#       return(list(NA, Occupied_cells))
# 
#   }