library(s2)
library(sf)
library(readr)

#-------------------------------------------
# Read data, DateTime to R datetime object
#-------------------------------------------
mw_data = read_delim("RAINLINK2/mw_data_subset.txt", delim = " ")
mw_data$DateTime = as.POSIXct(as.character(mw_data$DateTime),
                              format = "%Y%m%d%H%M")

#-------------------------------------------
# functions
#-------------------------------------------
CreateTowersSpatial <-  function(mw_data) {
  # Helper function to prepare spatial layer of MW towers
  # Create sf object of endpoints of links, using s2 
  IDs = unique(mw_data$ID)
  sf_use_s2(TRUE)
  towers_list = lapply(IDs, function(i) {
    # Create point sf object with two features, start and end of link
    link = mw_data[mw_data$ID == i,]
    # All time slots have the same coords. Get just the first:
    x1 = link$XStart[1]
    y1 = link$YStart[1]
    x2 = link$XEnd[1]
    y2 = link$YEnd[1]
    towers_df = data.frame("x" = c(x1, x2), "y" = c(y1, y2),
                           "ID"=i, "Side" = c("Start", "End"))
    
    # Create sf object of point features
    if (is.null(CoorSystemInputData)) {
      towers_sf = st_as_sf(towers_df, coords = c("x", "y"),
                           crs = 4326, agr = "constant")
    } else {
      # User supplied a CRS
      # TODO: Add to documentation that user CRS must be an EPSG code
      towers_sf = st_as_sf(towers_df, coords = c("x", "y"),
                           crs = CoorSystemInputData, agr = "constant")
    }
    return(towers_sf)
  })
  towers = do.call(rbind, towers_list)
  return(towers)
}

ApplyStep8 = function(ID, DateTime, interval) {
  # Helper Function to apply "Step 8"
  # Sets adjacent time slots to "wet"
  prev_1 = which(mw_data$ID == ID &
                   mw_data$DateTime == DateTime-interval)
  prev_2 = which(mw_data$ID == ID &
                   mw_data$DateTime == DateTime-(2*interval))
  next_1 = which(mw_data$ID == ID &
                   mw_data$DateTime == DateTime+interval)
  Dry_F$Dry[prev_1] <- 1
  Dry_F$Dry[prev_2] <- 1
  Dry_F$Dry[next_1] <- 1
}

# Plot function for testing/visualizing Nearby Links
PlotTowers <- function(tfull, tsubset, tlink) {
  library(leaflet)
  library(htmltools)
  m = leaflet() %>% addTiles()  %>%
    addCircleMarkers(data=tsubset, radius = 5,
                     color = "red", fillColor = "red",
                     fillOpacity = 0.4) %>%
    addCircleMarkers(data=tfull, radius = 10,
                     color="grey", fill = FALSE,
                     popup = ~htmlEscape(ID)) %>%
    addMarkers(data = tlink)
  m
}

#-------------------------------------------
# WetDryNearbyLinks function
#-------------------------------------------
WetDryNearbyLinksNew <- function(mw_data,
                                CoorSystemInputData=NULL,
                                MinHoursPmin=6,
                                PeriodHoursPmin=24,
                                Radius=10000,
                                Step8=TRUE,
                                ThresholdMedian=-1.4,
                                ThresholdMedianL=-0.7,
                                ThresholdNumberLinks=3,
                                ThresholdWetDry=2)  {

  towers = CreateTowersSpatial(mw_data)
  # Prepare distance matrix (one time operation)
  if (is.null(CoorSystemInputData)) {
    dist_matrix = s2_distance_matrix(towers, towers)
  } else {
    sf_use_s2(FALSE)
    dist_matrix = st_distance(towers, towers)
  }
  
  # Loop thru all rows of data,
  # For each row, get Pmin for nearby towers (from dist_matrix) 
  # and time slots during last PeriodHoursPmin hours
  n_rows = nrow(mw_data)
  Dry_F <- data.frame("Dry" = rep(NA, n_rows),
                      "F" = rep(NA, n_rows))
  # Determine time interval length (in seconds)
  uniq_times = sort(unique(strftime(mw_data$DateTime, "%s")))
  this_interval <- min(diff(as.numeric(uniq_times)))
  
  for (i in 1:n_rows) {
    this_row = mw_data[i,]
    this_ID = this_row$ID
    this_Pmin = this_row$Pmin
    this_Pmax = this_row$Pmax
    this_datetime = this_row$DateTime
    this_Length = this_row$PathLength  # PathLength is in km

    # Which tower ID matches this data row
    # Should return 2 rows, Start and End
    idx = as.numeric(row.names(towers[towers$ID == this_ID,]))
    # Build list of ID's that are within Radius, using the dist_matrix
    idx_start = which(dist_matrix[idx[1],] < Radius)
    idx_end = which(dist_matrix[idx[2],] < Radius)
    all_idx = c(idx_start, idx_end)
    tower_IDs_near = unique(towers[all_idx,]$ID)
    
    # we want data from PeriodHoursPmin ago: 
    prev_datetime = this_datetime - PeriodHoursPmin * 3600
    
    # Now subset original dataframe both by time (from MinHoursPmin ago)
    # and by spatial location, within Radius meters from this link
    mw_subset = mw_data[mw_data$DateTime > prev_datetime &
                          mw_data$DateTime <= this_datetime,]
    mw_subset = mw_subset[mw_subset$ID %in% tower_IDs_near,]
    #-------------------------------------------
    # Test/plot location of nearby links
    # Create spatial objects for plotting
    #towers_subset = CreateTowersSpatial(mw_subset)
    #tower_pair = CreateTowersSpatial(mw_subset[mw_subset$ID == this_ID,])
    #PlotTowers(towers, towers_subset, tower_pair)
    #-------------------------------------------

    # Do we have enough rows in the subset
    uniq_IDs = unique(mw_subset[mw_subset$DateTime == this_datetime,]$ID)
    if (length(uniq_IDs) > ThresholdNumberLinks) {
      # Calculate median of delta between Pmin and maximum Pmin value,
      # and determine Dry value [0,1]
      Pmin_max = max(mw_subset$Pmin, na.rm = TRUE)
      DeltaP = mw_subset$Pmin - Pmin_max
      DeltaP_L = DeltaP / this_Length
      medianDeltaP = median(DeltaP, na.rm = TRUE)
      medianDeltaP_L = median(DeltaP_L, na.rm = TRUE)
      if (medianDeltaP >= ThresholdMedian | 
          medianDeltaP_L >= ThresholdMedianL) {
        dry_val = 1
        # Should we apply Step8 to set previous and next time slots to wet?
        if (Step8==TRUE & (Pmin_max - this_Pmin > ThresholdWetDry)) {
          ApplyStep8(this_ID, this_datetime, this_interval)
        }
      } else if (is.na(medianDeltaP) | is.na(medianDeltaP_L)) {
        dry_val = NA
      } else {
        dry_val = 0
      }
      F_val <- (sum(DeltaP_L) - medianDeltaP_L) * this_interval / 3600
      Dry_F$Dry[i] = dry_val
      Dry_F$F[i] = F_val
    }
  }
  return(Dry_F)
}

mw_data_wd = WetDryNearbyLinksNew(mw_data)
