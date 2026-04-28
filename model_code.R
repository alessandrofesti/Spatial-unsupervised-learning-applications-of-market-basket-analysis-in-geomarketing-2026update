############################################################
# CHAPTER 11 - SPATIAL MARKET BASKET ANALYSIS
# Reproducible clean script aligned with the chapter
# Author: Alessandro Festi
############################################################

rm(list = ls())
set.seed(123)

############################################################
# 0. PACKAGES
############################################################

# Load all packages required by the workflow.
required_packages <- c(
  "MASS",
  "dplyr",
  "tidyr",
  "sf",
  "sp",
  "leaflet",
  "leaflet.extras2",
  "arules",
  "arulesViz",
  "arulesSequences",
  "osrm",
  "htmlwidgets",
  "webshot2",
  "knitr",
  "rmarkdown"
)

to_install <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)

invisible(lapply(required_packages, library, character.only = TRUE))

############################################################
# 1. PARAMETERS AND PROJECT FOLDERS
############################################################

# Simulation parameters.
n_individuals <- 50
n_paths <- 20
n_locations <- 100

# Reference point: Palazzo d'Accursio, Bologna.
palazzo_accursio_latitude <- 44.493674
palazzo_accursio_longitude <- 11.342220

# Path generation and spatial matching parameters.
variability_constant <- 200
correlation <- 0.7
step_scale <- 0.0003
max_match_distance_km <- 0.025   # 25 metres

# Project folders.
input_file <- file.path("data", "commercials.csv")
output_dir <- "outputs"
figure_dir <- "Images"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# Save a leaflet map both as HTML and PNG.
save_leaflet_png <- function(map_object, file_name, width = 1200, height = 800) {
  html_file <- normalizePath(
    file.path(figure_dir, paste0(file_name, ".html")),
    mustWork = FALSE
  )
  
  png_file <- normalizePath(
    file.path(figure_dir, paste0(file_name, ".png")),
    mustWork = FALSE
  )
  
  htmlwidgets::saveWidget(
    map_object,
    html_file,
    selfcontained = TRUE
  )
  
  webshot2::webshot(
    url = html_file,
    file = png_file,
    vwidth = width,
    vheight = height,
    delay = 1
  )
}

# Save the first five rows of an object as a CSV preview.
write_head5 <- function(x, file_name) {
  write.csv(
    head(x, 5),
    file.path(output_dir, file_name),
    row.names = FALSE
  )
}

# Stop early if the input file is missing.
if (!file.exists(input_file)) {
  stop(
    "Input file not found. Put commercials.csv in the data/ folder: ",
    normalizePath(input_file, mustWork = FALSE)
  )
}

############################################################
# 1B. STYLISED EXAMPLES OF POSSIBLE PATHS
############################################################

# Create simple stylised examples of spatial paths.
x <- runif(50, -1, 1)
y1 <- 3 * x + 5
y2 <- 12 * x^4 - 10 * x^2 + x - 4
y3 <- 2 * x^7 + 7

png(file.path(figure_dir, "figure_11_3_stylised_paths.png"), width = 1200, height = 800)
par(mfrow = c(1, 3))

plot(x, y1, xlab = "x", ylab = "y", pch = 16, main = "Linear path")
lines(x[order(x)], y1[order(x)])

plot(x, y2, xlab = "x", ylab = "y", pch = 16, main = "Curved path")
lines(x[order(x)], y2[order(x)])

plot(x, y3, xlab = "x", ylab = "y", pch = 16, main = "Non-linear path")
lines(x[order(x)], y3[order(x)])

par(mfrow = c(1, 1))
dev.off()

############################################################
# 2. SIMULATE INDIVIDUAL CENTRES AROUND PALAZZO D'ACCURSIO
############################################################

# Generate individual centres around the city centre.
individual_centres <- data.frame(
  lon = rnorm(n_individuals, mean = 0, sd = 1) / variability_constant + palazzo_accursio_longitude,
  lat = rnorm(n_individuals, mean = 0, sd = 1) / variability_constant + palazzo_accursio_latitude
)

write_head5(
  individual_centres,
  "output_11_0_individual_centres_preview.csv"
)

# Map the simulated individual centres.
centres_map <- leaflet(individual_centres) %>%
  addProviderTiles("CartoDB.Positron") %>%
  setView(lng = palazzo_accursio_longitude, lat = palazzo_accursio_latitude, zoom = 15) %>%
  addCircleMarkers(
    lng = ~lon,
    lat = ~lat,
    radius = 4,
    color = "red",
    fillColor = "red",
    fillOpacity = 0.8,
    stroke = FALSE,
    popup = ~paste0("Individual centre")
  ) %>%
  addCircleMarkers(
    lng = palazzo_accursio_longitude,
    lat = palazzo_accursio_latitude,
    radius = 7,
    color = "black",
    fillColor = "black",
    fillOpacity = 1,
    popup = "Palazzo d'Accursio",
    label = "Palazzo d'Accursio"
  )

save_leaflet_png(centres_map, "figure_11_4_individual_centres")

############################################################
# 3. SIMULATE CORRELATED PATHS
############################################################

# Generate serially correlated spatial paths.
generate_paths <- function(individual_centres,
                           n_paths,
                           n_locations,
                           correlation = 0.7,
                           step_scale = 0.0003) {
  time_index <- seq_len(n_locations)
  
  sigma_time <- outer(
    time_index,
    time_index,
    function(i, j) correlation^abs(i - j)
  )
  
  all_paths <- vector("list", length = nrow(individual_centres))
  
  for (j in seq_len(nrow(individual_centres))) {
    individual_paths <- vector("list", length = n_paths)
    
    for (i in seq_len(n_paths)) {
      start_lon <- individual_centres$lon[j] + rnorm(1, 0, step_scale * 1.5)
      start_lat <- individual_centres$lat[j] + rnorm(1, 0, step_scale * 1.5)
      
      delta_lon <- as.numeric(
        MASS::mvrnorm(n = 1, mu = rep(0, n_locations), Sigma = sigma_time)
      ) * step_scale
      
      delta_lat <- as.numeric(
        MASS::mvrnorm(n = 1, mu = rep(0, n_locations), Sigma = sigma_time)
      ) * step_scale
      
      individual_paths[[i]] <- data.frame(
        Longitude = start_lon + cumsum(delta_lon),
        Latitude = start_lat + cumsum(delta_lat),
        Individual_ID = j,
        Path_ID = i,
        Location_ID = seq_len(n_locations)
      )
    }
    
    all_paths[[j]] <- bind_rows(individual_paths)
  }
  
  bind_rows(all_paths)
}

paths <- generate_paths(
  individual_centres = individual_centres,
  n_paths = n_paths,
  n_locations = n_locations,
  correlation = correlation,
  step_scale = step_scale
)

write_head5(
  paths,
  "output_11_1_paths_head.csv"
)

############################################################
# 4. MAP OF TWO SAMPLE SIMULATED PATHS
############################################################

# Select two paths for visual inspection.
sample_paths <- paths %>%
  filter((Individual_ID == 1 & Path_ID == 1) | (Individual_ID == 2 & Path_ID == 1)) %>%
  arrange(Individual_ID, Path_ID, Location_ID)

sample_path_map <- leaflet() %>%
  addProviderTiles("CartoDB.Positron")

for (id in c(1, 2)) {
  path_subset <- sample_paths %>%
    filter(Individual_ID == id) %>%
    arrange(Location_ID)
  
  line_color <- if (id == 1) "blue" else "red"
  
  sample_path_map <- sample_path_map %>%
    addPolylines(
      lng = path_subset$Longitude,
      lat = path_subset$Latitude,
      color = line_color,
      weight = 4,
      opacity = 0.85,
      popup = paste("Individual", id)
    ) %>%
    addCircleMarkers(
      lng = path_subset$Longitude,
      lat = path_subset$Latitude,
      radius = 2,
      color = line_color,
      fillColor = line_color,
      fillOpacity = 0.7,
      stroke = FALSE
    )
}

start_points <- sample_paths %>%
  group_by(Individual_ID) %>%
  slice_min(Location_ID, n = 1) %>%
  ungroup()

end_points <- sample_paths %>%
  group_by(Individual_ID) %>%
  slice_max(Location_ID, n = 1) %>%
  ungroup()

sample_path_map <- sample_path_map %>%
  addCircleMarkers(
    data = start_points,
    lng = ~Longitude,
    lat = ~Latitude,
    radius = 6,
    color = "darkgreen",
    fillColor = "limegreen",
    fillOpacity = 1,
    weight = 2,
    popup = ~paste("Start - Individual", Individual_ID)
  ) %>%
  addCircleMarkers(
    data = end_points,
    lng = ~Longitude,
    lat = ~Latitude,
    radius = 6,
    color = "black",
    fillColor = "black",
    fillOpacity = 1,
    weight = 2,
    popup = ~paste("End - Individual", Individual_ID)
  ) %>%
  addLegend(
    "bottomright",
    colors = c("blue", "red", "limegreen", "black"),
    labels = c("Individual 1", "Individual 2", "Start", "End"),
    title = "Legend"
  )

save_leaflet_png(sample_path_map, "figure_11_5_sample_paths")

############################################################
# 5. LOAD COMMERCIAL ACTIVITIES DATASET
############################################################

# Read and standardise the commercial activities dataset.
commercial_activities <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)

names(commercial_activities) <- trimws(names(commercial_activities))
names(commercial_activities) <- gsub("\\.+", "_", names(commercial_activities))

if ("Latitudine" %in% names(commercial_activities)) {
  names(commercial_activities)[names(commercial_activities) == "Latitudine"] <- "lat"
}
if ("Longitudine" %in% names(commercial_activities)) {
  names(commercial_activities)[names(commercial_activities) == "Longitudine"] <- "lon"
}
if ("LATITUDINE" %in% names(commercial_activities)) {
  names(commercial_activities)[names(commercial_activities) == "LATITUDINE"] <- "lat"
}
if ("LONGITUDINE" %in% names(commercial_activities)) {
  names(commercial_activities)[names(commercial_activities) == "LONGITUDINE"] <- "lon"
}

required_cols <- c("lat", "lon")
missing_cols <- setdiff(required_cols, names(commercial_activities))

if (length(missing_cols) > 0) {
  stop("Missing columns in commercials.csv: ", paste(missing_cols, collapse = ", "))
}

commercial_activities <- commercial_activities %>%
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  ) %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  mutate(indexed_company = row_number())

commercial_activities_clean <- commercial_activities %>%
  mutate(matched_idx = row_number())

write_head5(
  commercial_activities_clean,
  "output_11_0_commercial_activities_preview.csv"
)

# Map all commercial activities.
commercial_activities_sf <- st_as_sf(
  commercial_activities,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

commercial_map <- leaflet(commercial_activities_sf) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(
    radius = 3,
    stroke = FALSE,
    fillOpacity = 0.6
  )

save_leaflet_png(commercial_map, "figure_11_6_commercial_activities")

############################################################
# 6. MATCH SIMULATED POINTS TO NEAREST COMMERCIAL ACTIVITY
############################################################

# Convert path and activity coordinates into numeric matrices.
path_matrix <- as.matrix(data.frame(
  lon = as.numeric(paths$Longitude),
  lat = as.numeric(paths$Latitude)
))

activity_matrix <- as.matrix(data.frame(
  lon = as.numeric(commercial_activities$lon),
  lat = as.numeric(commercial_activities$lat)
))

activity_matrix <- activity_matrix[complete.cases(activity_matrix), , drop = FALSE]

# Assign the nearest activity to each simulated path point.
assign_nearest_activity <- function(path_matrix,
                                    activity_matrix,
                                    max_match_distance_km = 0.025) {
  stopifnot(is.matrix(path_matrix), ncol(path_matrix) == 2)
  stopifnot(is.matrix(activity_matrix), ncol(activity_matrix) == 2)
  
  matched_points <- data.frame(
    distance_km = numeric(nrow(path_matrix)),
    matched_idx = rep(NA_integer_, nrow(path_matrix)),
    matched_lon = rep(NA_real_, nrow(path_matrix)),
    matched_lat = rep(NA_real_, nrow(path_matrix))
  )
  
  for (j in seq_len(nrow(path_matrix))) {
    distance_vector <- spDistsN1(
      pts = activity_matrix,
      pt = path_matrix[j, ],
      longlat = TRUE
    )
    
    nearest_idx <- which.min(distance_vector)
    nearest_distance <- distance_vector[nearest_idx]
    
    matched_points$distance_km[j] <- nearest_distance
    
    if (nearest_distance < max_match_distance_km) {
      matched_points$matched_idx[j] <- nearest_idx
      matched_points$matched_lon[j] <- activity_matrix[nearest_idx, 1]
      matched_points$matched_lat[j] <- activity_matrix[nearest_idx, 2]
    }
  }
  
  matched_points
}

matched_points <- assign_nearest_activity(
  path_matrix = path_matrix,
  activity_matrix = activity_matrix,
  max_match_distance_km = max_match_distance_km
)

matched_data <- data.frame(
  matched_idx = matched_points$matched_idx,
  matched_lon = matched_points$matched_lon,
  matched_lat = matched_points$matched_lat,
  distance_km = matched_points$distance_km,
  Individual_ID = paths$Individual_ID,
  Path_ID = paths$Path_ID,
  Location_ID = paths$Location_ID
) %>%
  tidyr::drop_na(matched_idx)

matched_data <- matched_data %>%
  left_join(commercial_activities_clean, by = "matched_idx") %>%
  rename(
    activity_lon = lon,
    activity_lat = lat
  ) %>%
  mutate(itemset_id = paste(Individual_ID, Path_ID, sep = "_"))

matched_data_preview <- matched_data %>%
  select(
    Individual_ID,
    Path_ID,
    Location_ID,
    itemset_id,
    indexed_company,
    activity_lon,
    activity_lat,
    distance_km
  ) %>%
  arrange(Individual_ID, Path_ID, Location_ID)

write_head5(
  matched_data_preview,
  "output_11_2_matched_data_preview.csv"
)

############################################################
# 7. SPATIAL ASSOCIATION RULES
############################################################

# Build path-level baskets from matched activities.
basket_data <- matched_data %>%
  select(
    itemset_id,
    Individual_ID,
    Path_ID,
    indexed_company,
    lon = activity_lon,
    lat = activity_lat,
    distance_km,
    Location_ID
  ) %>%
  distinct(itemset_id, indexed_company, .keep_all = TRUE)

write_head5(
  basket_data,
  "output_11_3_basket_data_preview.csv"
)

# Collapse each path into a comma-separated item list.
transaction_data <- basket_data %>%
  group_by(itemset_id) %>%
  summarise(
    items = paste(indexed_company, collapse = ","),
    .groups = "drop"
  )

write_head5(
  transaction_data,
  "output_11_4_transaction_data_preview.csv"
)

transaction_list <- strsplit(transaction_data$items, ",")
transactions <- as(transaction_list, "transactions")

# Extract spatial association rules with Apriori.
association_rules <- apriori(
  transactions,
  parameter = list(
    supp = 0.01,
    conf = 0.30,
    minlen = 2,
    maxlen = 4
  )
)

association_rules_summary <- capture.output(
  summary(association_rules)
)

writeLines(
  association_rules_summary,
  con = file.path(output_dir, "output_11_5a_association_rules_summary.txt")
)

association_rules_lift <- sort(association_rules, by = "lift")

association_rules_table <- as(association_rules, "data.frame") %>%
  arrange(desc(support), desc(confidence), desc(lift)) %>%
  head(5)

write.csv(
  association_rules_table,
  file.path(output_dir, "output_11_5_top_spatial_rules.csv"),
  row.names = FALSE
)

# Save graph of the top 20 spatial association rules.
png(
  filename = file.path(figure_dir, "figure_11_9_association_rules_graph_top20.png"),
  width = 1400,
  height = 1000,
  res = 150
)

plot(
  head(association_rules_lift, 20),
  method = "graph",
  engine = "igraph"
)

dev.off()

# Save graph of the top 100 spatial association rules.
png(
  filename = file.path(figure_dir, "figure_11_9_association_rules_graph_top100.png"),
  width = 1800,
  height = 1300,
  res = 150
)

plot(
  head(association_rules_lift, 100),
  method = "graph",
  engine = "igraph"
)

dev.off()

# Backward-compatible figure name used in older README versions.
file.copy(
  from = file.path(figure_dir, "figure_11_9_association_rules_graph_top20.png"),
  to = file.path(figure_dir, "figure_11_9_association_rules_graph.png"),
  overwrite = TRUE
)

############################################################
# 8. LOOKUP TABLE FOR COMMERCIAL ACTIVITIES
############################################################

# Export a lookup table to interpret activity IDs.
lookup_cols <- intersect(
  c("indexed_company", "Ubicazione", "Quartiere", "Zona", "Sottoarea"),
  names(commercial_activities_clean)
)

activity_lookup <- commercial_activities_clean %>%
  select(all_of(lookup_cols))

write.csv(
  activity_lookup,
  file.path(output_dir, "table_11_3_activity_lookup.csv"),
  row.names = FALSE
)

############################################################
# 9. BUILD ORDERED SEQUENCES
############################################################

# Keep the first occurrence of each activity inside each path.
sequence_data <- matched_data %>%
  arrange(Individual_ID, Path_ID, Location_ID) %>%
  group_by(itemset_id, indexed_company) %>%
  summarise(
    first_location_id = min(Location_ID),
    lon = first(activity_lon),
    lat = first(activity_lat),
    .groups = "drop"
  ) %>%
  separate(
    itemset_id,
    into = c("individual_chr", "path_chr"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    individual_num = as.integer(individual_chr),
    path_num = as.integer(path_chr)
  ) %>%
  arrange(individual_num, path_num, first_location_id) %>%
  group_by(itemset_id, individual_num, path_num) %>%
  mutate(
    order_in_path = row_number(),
    sequenceID = cur_group_id()
  ) %>%
  ungroup()

sequence_preview <- sequence_data %>%
  select(
    itemset_id,
    indexed_company,
    first_location_id,
    lon,
    lat,
    order_in_path,
    sequenceID
  ) %>%
  mutate(
    lon = round(lon, 6),
    lat = round(lat, 6)
  )

write_head5(
  sequence_preview,
  "output_11_6_sequence_preview.csv"
)

############################################################
# 10. SEQUENTIAL SPATIAL RULES
############################################################

# Convert ordered paths into the format required by arulesSequences.
sequence_input <- sequence_data %>%
  arrange(sequenceID, order_in_path) %>%
  transmute(
    sequenceID = sequenceID,
    eventID = order_in_path,
    SIZE = 1,
    items = as.character(indexed_company)
  )

write_head5(
  sequence_input,
  "output_11_7_sequence_input_preview.csv"
)

tmp_sequence_file <- tempfile(fileext = ".txt")

write.table(
  sequence_input,
  file = tmp_sequence_file,
  sep = " ",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

sequence_transactions <- read_baskets(
  con = tmp_sequence_file,
  info = c("sequenceID", "eventID", "SIZE")
)

# Extract sequential patterns and derive sequential rules.
sequential_patterns <- cspade(
  sequence_transactions,
  parameter = list(support = 0.005),
  control = list(verbose = TRUE)
)

sequential_rules <- ruleInduction(
  sequential_patterns,
  confidence = 0.15
)

# Inspect the strongest sequential rules by confidence, including rules
# that may contain more than one item or more than one event.
top_sequential_rules_by_confidence <- sort(
  sequential_rules,
  by = "confidence",
  decreasing = TRUE
) %>%
  head(5)

top_sequential_rules_by_confidence_df <- as(
  top_sequential_rules_by_confidence,
  "data.frame"
)

write.csv(
  top_sequential_rules_by_confidence_df,
  file.path(output_dir, "output_11_8a_top_sequential_rules_by_confidence.csv"),
  row.names = FALSE
)

# Extract rules with one item on the left and one item on the right.
extract_single_rule_id <- function(label_value) {
  ids <- regmatches(label_value, gregexpr("[0-9]+", label_value))[[1]]
  if (length(ids) == 1) as.integer(ids) else NA_integer_
}

sequential_rules_df <- as(sequential_rules, "data.frame")

sequential_rules_labels <- data.frame(
  rule = sequential_rules_df$rule,
  support = sequential_rules_df$support,
  confidence = sequential_rules_df$confidence,
  lift = sequential_rules_df$lift,
  lhs_id = vapply(
    as.character(labels(lhs(sequential_rules))),
    extract_single_rule_id,
    integer(1)
  ),
  rhs_id = vapply(
    as.character(labels(rhs(sequential_rules))),
    extract_single_rule_id,
    integer(1)
  )
)

top_simple_sequential_rules <- sequential_rules_labels %>%
  filter(!is.na(lhs_id), !is.na(rhs_id)) %>%
  arrange(desc(confidence), desc(support), desc(lift)) %>%
  head(5)

write.csv(
  top_simple_sequential_rules,
  file.path(output_dir, "output_11_8_top_sequential_rules.csv"),
  row.names = FALSE
)

############################################################
# 11. REPRESENTATIVE ORDERED PATH VISUALISATION
############################################################

# Select a path with a readable number of matched activities.
path_sizes <- sequence_data %>%
  group_by(itemset_id) %>%
  summarise(
    n_companies = n_distinct(indexed_company),
    .groups = "drop"
  ) %>%
  arrange(n_companies)

candidate_paths <- path_sizes %>%
  filter(n_companies >= 5, n_companies <= 8)

selected_itemset <- if (nrow(candidate_paths) > 0) {
  candidate_paths$itemset_id[1]
} else {
  path_sizes$itemset_id[1]
}

ordered_path <- sequence_data %>%
  filter(itemset_id == selected_itemset) %>%
  arrange(order_in_path)

start_point <- ordered_path %>% slice(1)
end_point <- ordered_path %>% slice(n())

middle_points <- if (nrow(ordered_path) > 2) {
  ordered_path %>% slice(2:(n() - 1))
} else {
  ordered_path[0, ]
}

# Map the ordered path with start, middle and end points.
ordered_path_map <- leaflet() %>%
  addProviderTiles("CartoDB.Positron") %>%
  addPolylines(
    data = ordered_path,
    lng = ~lon,
    lat = ~lat,
    color = "blue",
    weight = 4,
    opacity = 0.85
  )

if (nrow(middle_points) > 0) {
  ordered_path_map <- ordered_path_map %>%
    addCircleMarkers(
      data = middle_points,
      lng = ~lon,
      lat = ~lat,
      radius = 5,
      color = "red",
      fillColor = "red",
      fillOpacity = 0.95,
      popup = ~paste0(
        "<b>Company:</b> ", indexed_company,
        "<br><b>Order:</b> ", order_in_path
      )
    )
}

ordered_path_map <- ordered_path_map %>%
  addCircleMarkers(
    data = start_point,
    lng = ~lon,
    lat = ~lat,
    radius = 7,
    color = "darkgreen",
    fillColor = "limegreen",
    fillOpacity = 1,
    weight = 2,
    popup = ~paste0(
      "<b>START</b><br><b>Company:</b> ", indexed_company,
      "<br><b>Order:</b> ", order_in_path
    )
  ) %>%
  addCircleMarkers(
    data = end_point,
    lng = ~lon,
    lat = ~lat,
    radius = 7,
    color = "black",
    fillColor = "black",
    fillOpacity = 1,
    weight = 2,
    popup = ~paste0(
      "<b>END</b><br><b>Company:</b> ", indexed_company,
      "<br><b>Order:</b> ", order_in_path
    )
  ) %>%
  addLabelOnlyMarkers(
    data = ordered_path,
    lng = ~lon,
    lat = ~lat,
    label = ~as.character(order_in_path),
    labelOptions = labelOptions(
      noHide = TRUE,
      textOnly = TRUE,
      direction = "top",
      offset = c(0, -10)
    )
  ) %>%
  addLegend(
    "bottomright",
    colors = c("blue", "red", "limegreen", "black"),
    labels = c("Ordered path", "Visited companies", "Start", "End"),
    title = "Legend"
  ) %>%
  fitBounds(
    min(ordered_path$lon, na.rm = TRUE) - 0.0001,
    min(ordered_path$lat, na.rm = TRUE) - 0.0001,
    max(ordered_path$lon, na.rm = TRUE) + 0.0001,
    max(ordered_path$lat, na.rm = TRUE) + 0.0001
  )

save_leaflet_png(ordered_path_map, "figure_11_7_ordered_path")

############################################################
# 12. SAME ORDERED PATH ON THE PEDESTRIAN ROAD NETWORK
############################################################

# Use OSRM only to draw a road-network route between the already ordered points.
options(osrm.server = "https://routing.openstreetmap.de/")
options(osrm.profile = "foot")

route_points <- ordered_path %>%
  select(lon, lat)

road_route <- tryCatch(
  osrmRoute(
    loc = route_points,
    overview = "full"
  ),
  error = function(e) {
    message("osrmRoute failed: ", conditionMessage(e))
    NULL
  }
)

ordered_path_road_map <- leaflet() %>%
  addProviderTiles("CartoDB.Positron")

if (!is.null(road_route)) {
  ordered_path_road_map <- ordered_path_road_map %>%
    addPolylines(
      data = road_route,
      color = "blue",
      weight = 4,
      opacity = 0.9
    )
}

ordered_path_road_map <- ordered_path_road_map %>%
  addCircleMarkers(
    data = ordered_path,
    lng = ~lon,
    lat = ~lat,
    radius = 5,
    color = "red",
    fillColor = "red",
    fillOpacity = 0.95,
    popup = ~paste0(
      "<b>Company:</b> ", indexed_company,
      "<br><b>Order:</b> ", order_in_path
    )
  ) %>%
  addLabelOnlyMarkers(
    data = ordered_path,
    lng = ~lon,
    lat = ~lat,
    label = ~as.character(order_in_path),
    labelOptions = labelOptions(
      noHide = TRUE,
      textOnly = TRUE,
      direction = "top",
      offset = c(0, -10)
    )
  ) %>%
  fitBounds(
    min(ordered_path$lon, na.rm = TRUE) - 0.0001,
    min(ordered_path$lat, na.rm = TRUE) - 0.0001,
    max(ordered_path$lon, na.rm = TRUE) + 0.0001,
    max(ordered_path$lat, na.rm = TRUE) + 0.0001
  )

save_leaflet_png(ordered_path_road_map, "figure_11_8_ordered_path_osrm")

############################################################
# 13. MAP OF TOP SIMPLE SEQUENTIAL RULES
############################################################

seq_rules_df <- as(sequential_rules, "data.frame")
seq_rules_df$rule_index <- seq_len(nrow(seq_rules_df))

simple_rule_indices <- seq_rules_df %>%
  mutate(
    lhs_id = vapply(
      as.character(labels(lhs(sequential_rules))),
      extract_single_rule_id,
      integer(1)
    ),
    rhs_id = vapply(
      as.character(labels(rhs(sequential_rules))),
      extract_single_rule_id,
      integer(1)
    )
  ) %>%
  filter(!is.na(lhs_id), !is.na(rhs_id))


rules_to_map <- simple_rule_indices %>%
  left_join(
    commercial_activities_clean %>%
      select(lhs_id = indexed_company, lhs_lon = lon, lhs_lat = lat),
    by = "lhs_id"
  ) %>%
  left_join(
    commercial_activities_clean %>%
      select(rhs_id = indexed_company, rhs_lon = lon, rhs_lat = lat),
    by = "rhs_id"
  ) %>%
  filter(!is.na(lhs_lon), !is.na(rhs_lon)) %>%
  rowwise() %>%
  mutate(
    rule_distance_km = sp::spDistsN1(
      pts = matrix(c(rhs_lon, rhs_lat), ncol = 2),
      pt = c(lhs_lon, lhs_lat),
      longlat = TRUE
    )
  ) %>%
  ungroup() %>%
  arrange(desc(rule_distance_km), desc(confidence), desc(support), desc(lift)) %>%
  slice(1)

sequential_rule_map <- leaflet() %>%
  addProviderTiles("CartoDB.Positron")

if (nrow(rules_to_map) > 0) {
  for (i in seq_len(nrow(rules_to_map))) {
    sequential_rule_map <- sequential_rule_map %>%
      addArrowhead(
        lng = c(rules_to_map$lhs_lon[i], rules_to_map$rhs_lon[i]),
        lat = c(rules_to_map$lhs_lat[i], rules_to_map$rhs_lat[i]),
        color = "black",
        weight = 3,
        opacity = 0.9,
        options = arrowheadOptions(frequency = "endonly", size = "15px"),
        popup = paste0(
          rules_to_map$lhs_id[i], " → ", rules_to_map$rhs_id[i],
          "<br>Support: ", round(rules_to_map$support[i], 4),
          "<br>Confidence: ", round(rules_to_map$confidence[i], 4),
          "<br>Lift: ", round(rules_to_map$lift[i], 2),
          "<br>Distance: ", round(rules_to_map$rule_distance_km[i], 3), " km"
        )
      )
  }
  
  all_points <- data.frame(
    indexed_company = c(rules_to_map$lhs_id, rules_to_map$rhs_id),
    lon = c(rules_to_map$lhs_lon, rules_to_map$rhs_lon),
    lat = c(rules_to_map$lhs_lat, rules_to_map$rhs_lat)
  ) %>%
    distinct(indexed_company, .keep_all = TRUE)
  
  sequential_rule_map <- sequential_rule_map %>%
    addCircleMarkers(
      data = all_points,
      lng = ~lon,
      lat = ~lat,
      radius = 5,
      color = "black",
      fillColor = "white",
      fillOpacity = 1,
      weight = 2,
      label = ~as.character(indexed_company)
    ) %>%
    fitBounds(
      min(all_points$lon, na.rm = TRUE) - 0.001,
      min(all_points$lat, na.rm = TRUE) - 0.001,
      max(all_points$lon, na.rm = TRUE) + 0.001,
      max(all_points$lat, na.rm = TRUE) + 0.001
    )
}

save_leaflet_png(sequential_rule_map, "figure_11_10_top_sequential_rules_map")

############################################################
# 14. SESSION INFO
############################################################

# Save session information for reproducibility.
capture.output(
  sessionInfo(),
  file = file.path(output_dir, "session_info.txt")
)

message("Script completed. Outputs saved in: ", normalizePath(output_dir))
