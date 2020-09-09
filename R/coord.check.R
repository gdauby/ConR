#' Internal function
#'
#' Coordinates check
#'
#' @param XY data.frame, of at least two columns (coordinates), third is taxa
#' @param listing logical, whether the dataset should be splitted in a list by taxa
#' @param proj_type ...
#' @param listing_by_valid ...
#'
#'
coord.check <-
  function(XY,
           listing = TRUE,
           proj_type = NULL,
           listing_by_valid = FALSE) {
    
    XY <- as.data.frame(XY)
    
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
    
    if (!is.null(proj_type)) {
      XY_proj <-
        sf::sf_project(
          from = sf::st_crs(4326),
          to =
            sf::st_crs(proj_type),
          pts = XY[, c(2, 1)]
        )[, c(2, 1)]
      
      XY[, c(1, 2)] <-
        XY_proj[, c(1, 2)]
      
    }
    
    
    if (listing) {
      if (ncol(XY) > 2) {
        colnames(XY)[1:3] <- c("ddlat", "ddlon", "tax")
        XY$tax <- as.character(XY$tax)
        
        if (length(grep("[?]", XY[, 3])) > 0)
          XY[, 3] <- gsub("[?]", "_", XY[, 3])
        if (length(grep("[/]", XY[, 3])) > 0)
          XY[, 3] <- gsub("[/]", "_", XY[, 3])
        
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
          list_data <- split(XY, f = XY$tax)
          
        }
        
      } else{
        colnames(XY)[1:2] <- c("ddlat", "ddlon")
        list_data <- list(XY)
        
      }
      
    } else{
      list_data <-
        XY
      
    }
    
    return(list_data)
  }