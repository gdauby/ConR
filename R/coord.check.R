#' Internal function
#'
#' Coordinates check
#'
#' @param XY data.frame, of at least two columns (coordinates), third is taxa
#' @param listing logical, whether the dataset should be splitted in a list by taxa
#' @param proj_type crs
#' @param listing_by_valid logical
#' @param check_eoo logical
#' 
#' 
#' @return a list
#' 
#' @keywords internal
#'
#'
coord.check <-
  function(XY,
           listing = TRUE,
           proj_type = NULL,
           listing_by_valid = FALSE,
           cell_size = NULL,
           check_eoo = TRUE) {
    
    XY <- as.data.frame(XY)
    
    if (ncol(XY) < 3) stop("At least three columns are expected in the following order : latitude, longitude and species names")
    
    
    # if (length(grep("[?]", XY[, 3])) > 0)
      XY[, 3] <- gsub("[?]", "_", XY[, 3])
    # if (length(grep("[/]", XY[, 3])) > 0)
      XY[, 3] <- gsub("[/]", "_", XY[, 3])
    
    # if (length(grep("\\(", XY[, 3])) > 0)
      XY[, 3] <- gsub("\\(", "", XY[, 3])
    # if (length(grep("\\)", XY[, 3])) > 0)
      XY[, 3] <- gsub("\\)", "", XY[, 3])
      
      
    if (!is.numeric(XY[,1])) XY[,1] <- as.numeric(XY[,1])
    
     
    if (any(is.na(XY[, c(1:2)]))) {
      print(
        paste(
          "Skipping",
          length(which(rowMeans(is.na(
            XY[, 1:2]
          )) > 0)) ,
          "occurrences because of missing coordinates for",
          # if(verbose)
          length(unique(XY[which(rowMeans(is.na(XY[, 1:2])) >
                                   0), 3])),
          "taxa"
        )
      )
      XY <- XY[which(!is.na(XY[, 1])),]
      XY <- XY[which(!is.na(XY[, 2])),]
    }
    
    if (any(XY[, 2] > 180) ||
        any(XY[, 2] < -180) ||
        any(XY[, 1] < -180) ||
        any(XY[, 1] > 180))
      stop("coordinates are outside of expected range")
    
    colnames(XY)[1:3] <- c("ddlat", "ddlon", "tax")
    XY$tax <- as.character(XY$tax)
    list_data <- split(XY, f = XY$tax)
    
    ### check distances to antimeredian
    if (!is.null(cell_size)) {
      
      check_dist_antimeridian <-  function(x, siz) {
        
        kk <- x
        
        if (any(x > 0))
          kk[kk > 0] <- 180 - x[x > 0]
        
        if (any(x <= 0))
          kk[kk <= 0] <- 180 - (-x[x <= 0])
        
        kk <- kk*40000/360
        
        return(any(kk <= siz*2))
      }
      
      issue_close_to_anti <- lapply(list_data, FUN = function(x) check_dist_antimeridian(x = x[,2], siz = cell_size))
      issue_close_to_anti <- unlist(issue_close_to_anti)
      issue_close_to_anti <- which(issue_close_to_anti)
      
    } else {
      
      issue_close_to_anti <- NA
      
    }
    
    if (check_eoo) {
      
      check_max_long_dist <- function(x) {
        
        if (nrow(x) > 1) {
          if (max(dist(x[, 2]), na.rm = T) >= 180) {
            ver <- TRUE
          } else {
            ver <- FALSE
          }        
        } else {
          ver <- FALSE
        }
        
        return(ver)      
      }
      
      check_dist <- lapply(list_data, FUN = function(x) check_max_long_dist(x = x))
      check_dist <- unlist(check_dist)
      issue_long_span <- which(check_dist)
      
      check_nrow <- lapply(list_data, FUN = function(x) {nrow(unique(x)) < 3})
      check_nrow <- unlist(check_nrow)
      issue_nrow <- which(check_nrow)
      
      
    } else {
      issue_long_span <- issue_nrow <- NA
    }
    
    if (!is.null(proj_type)) {
      XY_proj <-
        sf::sf_project(
          from = sf::st_crs(4326),
          to =
            sf::st_crs(proj_type),
          pts = XY[, c(2, 1)]
        )
      
      XY_proj <- as.data.frame(XY_proj)
      
      XY_proj <-  XY_proj[, c(2, 1)]
      
      XY[, c(1, 2)] <-
        XY_proj[, c(1, 2)]
      
    }
    
    if (listing) {
      

        
        if (listing_by_valid) {
          list_data <-
            vector('list', length(unique(XY$classes)) * length(unique(XY$tax)))
          
          for (i in sort(unique(XY$classes))) {
            tax_classes <- vector(mode = "character", length = nrow(XY))
            
            if (i == 1) {
              XY_subset <-
                data.frame(
                  XY[XY$classes == i, c(1, 2)],
                  tax = XY[XY$classes == i, which(colnames(XY) == "tax")],
                  valid = XY[XY$classes == i, which(colnames(XY) == "valid")],
                  recordID = XY[XY$recordID == i, which(colnames(XY) == "recordID")],
                  classes = XY[XY$classes == i, which(colnames(XY) == "classes")],
                  tax_class = paste(XY$tax[XY$classes == i], i, sep = "___")
                )
              
              id_list <-
                which(unlist(lapply(list_data, is.null)))[1:length(unique(XY_subset$tax_class))]
              
              splited_data <-
                split(XY_subset, f = XY_subset$tax_class)
              
              list_data[id_list] <-
                splited_data
              
              names(list_data)[id_list] <-
                names(splited_data)
              
            }
            
            if (i > 1) {
              XY_subset <-
                data.frame(
                  XY[XY$classes %in% seq(1, i, 1), c(1, 2)],
                  tax = XY[XY$classes %in% seq(1, i, 1), which(colnames(XY) == "tax")],
                  valid = XY[XY$classes %in% seq(1, i, 1), which(colnames(XY) == "valid")],
                  recordID = XY[XY$classes %in% seq(1, i, 1), which(colnames(XY) == "recordID")],
                  classes = XY[XY$classes %in% seq(1, i, 1), which(colnames(XY) == "classes")],
                  tax_class = paste(XY$tax[XY$classes %in% seq(1, i, 1)], paste(seq(1, i, 1), collapse = "_"), sep = "___")
                )
              
              id_list <-
                which(unlist(lapply(list_data, is.null)))[1:length(unique(XY_subset$tax_class))]
              
              splited_data <-
                split(XY_subset, f = XY_subset$tax_class)
              
              list_data[id_list] <-
                splited_data
              
              names(list_data)[id_list] <-
                names(splited_data)
              
              
            }
            
          }
          
          # list_data <-
          #   do.call("rbind", lapply(list_data, function(x) do.call("rbind", x)))
          
          
          # XY$tax_class <-
          #   paste(XY$tax, XY$classes, sep = "___")
          # list_data <-
          #   split(XY, f = XY$tax_class)
          
          list_data <- list_data[!unlist(lapply(list_data, is.null))]
          
        } else {
          
          list_data <- split(unique(XY), f = unique(XY)$tax)
          
        }
        
    } else{
      
      list_data <-
        XY
      
    }
    
    return(list(list_data = list_data, 
                issue_close_to_anti = issue_close_to_anti,
                issue_long_span = issue_long_span,
                issue_nrow = issue_nrow))
  }