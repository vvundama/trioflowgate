library(flowWorkspace) 
library(flowCore)
library(openCyto)
library(CytoML) 

input_base <- "/Users/vishnuvundamati/Desktop/REPLACE WITH FULL PATH TO .fcs_gs FOLDER "
output_dir <- "/Users/vishnuvundamati/Desktop/REPLACE WITH FULL PATH TO EMPTY FOLDER"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Parsing through input_base to find desired file type 
gs_files <- list.files(input_base, pattern = "\\.fcs_gs$", full.names = TRUE)

message("Found ", length(gs_files), " .fcs_gs files")

# Density-based threshold finder 
robust_mindensity <- function(values, channel_name, sample_name, gate_range = NULL) {
  values <- values[!is.na(values) & is.finite(values)]
  values <- values[values >= 0]
  
  if (!is.null(gate_range) && length(gate_range) == 2) {
    values <- values[values >= gate_range[1] & values <= gate_range[2]]
  }
  # Accommodating empty gate range 
  if (length(values) == 0) {
    warning("No values in gate range for ", channel_name, " in ", sample_name, " default to full range median.")
    values_full <- values[!is.na(values) & is.finite(values) & values >= 0]
    return(median(values_full))
  }
  # Accomodating too few unique values
  if (length(unique(values)) < 10) {
    warning("Too few unique values in ", channel_name, " for ", sample_name, " using median.")
    return(median(values))
  }
  # Fallback to first minimum in kernel density 
  thresh <- suppressWarnings(tryCatch({
    openCyto:::mindensity(values)
  }, error = function(e) NA))
  
  if (is.na(thresh)) {
    d <- density(values, bw = "nrd0", adjust = 1.5)
    local_min_idx <- which(diff(sign(diff(d$y))) == 2)
    if (length(local_min_idx) > 0) {
      thresh <- d$x[local_min_idx[1]]
    } else {
      warning("Fallback to median for ", channel_name, " in ", sample_name)
      thresh <- median(values)
    }
  }
  
  return(thresh)
}

# Creating NKT-like gate
.gate_CD3_CD56_nkt <- function(fr, channels) {
  sample_name <- identifier(fr)
  dat <- exprs(fr)[, channels]
  cd3_thresh <- robust_mindensity(dat[, 1], channels[1], sample_name)
  cd56_thresh <- robust_mindensity(dat[, 2], channels[2], sample_name, gate_range = c(75, 125))
  rectangleGate(filterId = "NKTlike",
                .gate = setNames(list(c(cd3_thresh, Inf), c(cd56_thresh, Inf)), channels))
}

# Creating Senescent T gate
.gate_CD3_CD57_senescent <- function(fr, channels) {
  sample_name <- identifier(fr)
  dat <- exprs(fr)[, channels]
  cd3_thresh <- robust_mindensity(dat[, 1], channels[1], sample_name)
  cd57_thresh <- robust_mindensity(dat[, 2], channels[2], sample_name, gate_range = c(113, 125))
  rectangleGate(filterId = "SenescentT",
                .gate = setNames(list(c(cd3_thresh, Inf), c(cd57_thresh, Inf)), channels))
}

# Process each .fcs_gs file
for (g in gs_files) {
  base_name <- sub("\\.fcs_gs$", "", basename(g))
  message("Processing: ", base_name)
  
  tryCatch({
    gs <- load_gs(g)
    
    # Get flowFrame from NKT-like & SenescentT parent "LiveCells"
    fr <- gh_pop_get_data(gs[[1]], "/nonDebris/Lymphocytes/SingletCells/LiveCells")
    
    # Define and add NKT-like & SenescentT gate 
    gate_nkt <- .gate_CD3_CD56_nkt(fr, channels = c("APC-A", "BV510-A"))
    gate_senescent <- .gate_CD3_CD57_senescent(fr, channels = c("APC-A", "PE-CF594-A"))
    
    if (!"NKTlike" %in% gs_get_pop_paths(gs)) {
      message("Adding NKTlike gate")
      gs_pop_add(gs, 
                 gate_nkt, 
                 parent = "/nonDebris/Lymphocytes/SingletCells/LiveCells")
    }
    
    if (!"SenescentT" %in% gs_get_pop_paths(gs)) {
      message("Adding SenescentT gate")
      gs_pop_add(gs, 
                 gate_senescent, 
                 parent = "/nonDebris/Lymphocytes/SingletCells/LiveCells")
    }
    
    # Define and add Double Positive T Cell gate from existing parent threshold
    gh <- gs[[1]]
    fr <- gh_pop_get_data(gh, "/nonDebris/Lymphocytes/SingletCells/LiveCells/Tcells")
    
    cd8_min<-gs_pop_get_gate(gh,"CytotoxicT")[[1]]@min[[2]]
    cd4_min<-gs_pop_get_gate(gh,"CytotoxicT")[[1]]@max[[1]]
    
    double_pos_gate <- rectangleGate(
      filterId = "DoublePositive",
      .gate = list("APC-Cy7-A" = c(cd4_min, Inf),"BUV395-A" = c(cd8_min, Inf))
    )
    
    message("Adding DoublePositive gate")
    
    gs_pop_add(
      gs, double_pos_gate,
      parent = "/nonDebris/Lymphocytes/SingletCells/LiveCells/Tcells",
      name = "DoublePositive")
    
    recompute(gs)
    
    # Saving the modified GatingSet with appended gates 
    save_gs(gs, file.path(output_dir, paste0(base_name, "_gs_appended")))
    message("All three gates appended: ", base_name)
    
  }, error = function(e) {
    message("Sample ", base_name, " failed: ", conditionMessage(e))
  })
  
}
