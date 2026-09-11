
library(ConR)

context("Test AOH.estimation with a multi-layer SpatRaster")

test_that("AOH.estimation uses every layer and keeps pixels equal to hab.class", {

  skip_if_not_installed("lwgeom")

  set.seed(1)
  # sp A in the west half of the map, sp B in the east half
  XY <- data.frame(ddlat = c(runif(10, -5, -1), runif(10, -5, -1)),
                   ddlon = c(runif(10, 20, 23), runif(10, 27, 30)),
                   tax = rep(c("sp A", "sp B"), each = 10))

  hab.map <- terra::rast(nrows = 140, ncols = 200,
                         xmin = 15, xmax = 35, ymin = -12, ymax = 2,
                         crs = "EPSG:4326", nlyrs = 2)
  # layer 1: class 18 in the west, class 10 in the east
  terra::values(hab.map[[1]]) <-
    ifelse(terra::xFromCell(hab.map, 1:terra::ncell(hab.map)) < 25, 18, 10)
  # layer 2: class 18 or NA on a random half of the pixels
  terra::values(hab.map[[2]]) <-
    sample(c(NA, 18), terra::ncell(hab.map), replace = TRUE)

  res <- suppressWarnings(suppressMessages(
    AOH.estimation(XY = XY,
                   hab.map = hab.map,
                   hab.class = 18,
                   years = c(1992, 2022),
                   hab.map.type = TRUE,
                   exclude.area = FALSE,
                   show_progress = FALSE,
                   simplifiy_poly = FALSE)
  ))

  AOH <- res$AOH[order(res$AOH$species), ]

  expect_true(all(c("hab.map_18_1992", "hab.map_18_2022") %in% names(AOH)))

  # layer 1: the whole EOO of sp A is habitat, none of sp B
  expect_equal(AOH$hab.map_18_1992[1] / AOH$eoo[1], 1, tolerance = 0.02)
  expect_equal(AOH$hab.map_18_1992[2], 0)

  # layer 2: about half of each EOO is habitat, classes are not blended
  # when the raster is projected
  prop_2022 <- AOH$hab.map_18_2022 / AOH$eoo
  expect_true(all(prop_2022 > 0.4 & prop_2022 < 0.6))
})
