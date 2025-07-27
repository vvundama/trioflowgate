setwd("/Users/vishnuvundamati/Desk...YOUR OWN DIRECTORY")

.libPaths(new="/Library/Frameworks/R.framework/Versions/4.5-arm64...REPLACE WITH R LOCATION")

library(openCyto)
library(data.table)
library(flowWorkspaceData)
library(flowWorkspace)
library(ncdfFlow)
library(CytoML)
library(ggcyto)
library(cytoUtils)
library(tidyverse)
library(flowCore)
library(flowStats)

# Manually calling openCyto tailgate function 
openCyto:::register_plugins(fun = cytoUtils:::.gate_tail,methodName="tailgate")

.truncate_flowframe <- function(fr, channels, min = NULL, max = NULL) {return(fr)}
openCyto::register_plugins(fun = .truncate_flowframe, methodName = ".truncate_flowframe")

# Robust mindensity helper function
robust_mindensity <- function(values, channel_name, sample_name) {
  values <- values[!is.na(values) & is.finite(values)]
  values <- values[values >= 0]
  
  if (length(unique(values)) < 10) {
    warning("Too few unique values in ", channel_name, " for ", sample_name, "- using median as threshold.")
    return(median(values))
  }
  
  thresh <- suppressWarnings(tryCatch({
    openCyto:::mindensity(values)
  }, error = function(e) NA))
  
  if (is.na(thresh)) {
    d <- density(values, bw = "nrd0", adjust = 1.5)
    local_min_idx <- which(diff(sign(diff(d$y))) == 2)
    
    if (length(local_min_idx) > 0) {
      thresh <- d$x[local_min_idx[1]]
    } else {
      warning("Fallback to median quantile for ", channel_name, " in ", sample_name)
      thresh <- median(values)
    }
  }
  
  return(thresh)
}

# NKT-like gating function
.gate_CD3_CD56_nkt <- function(fr, pp_res, ...) {
  cd3_channel <- "APC-A"
  cd56_channel <- "BV510-A"
  sample_name <- identifier(fr)
  
  cd3_vals <- exprs(fr)[, cd3_channel]
  cd56_vals <- exprs(fr)[, cd56_channel]
  
  cd3_thresh <- robust_mindensity(cd3_vals, cd3_channel, sample_name)
  cd56_thresh <- robust_mindensity(cd56_vals, cd56_channel, sample_name)
  
  gate <- rectangleGate(
    filterId = "NKTlike",
    .gate = setNames(list(c(cd3_thresh, Inf), c(cd56_thresh, Inf)), c(cd3_channel, cd56_channel))
  )
  
  return(setNames(list(gate), sample_name))
}
openCyto::register_plugins(fun = .gate_CD3_CD56_nkt, methodName = "gate_nkt")

# Senescent T gating function
.gate_CD3_CD57_senescent <- function(fr, pp_res, ...) {
  cd3_channel <- "APC-A"
  cd57_channel <- "PE-CF594-A"
  sample_name <- identifier(fr)
  
  cd3_vals <- exprs(fr)[, cd3_channel]
  cd57_vals <- exprs(fr)[, cd57_channel]
  
  cd3_thresh <- robust_mindensity(cd3_vals, cd3_channel, sample_name)
  cd57_thresh <- robust_mindensity(cd57_vals, cd57_channel, sample_name)
  
  gate <- rectangleGate(
    filterId = "SenescentT", .gate = setNames(list(c(cd3_thresh, Inf), c(cd57_thresh, Inf)), c(cd3_channel, cd57_channel))
  )
  
  return(setNames(list(gate), sample_name))
}
openCyto::register_plugins(fun = .gate_CD3_CD57_senescent, methodName = "gate_senescent_t")

gtFile <- "/Users/vishnuvundamati/Desktop/ENTER YOU OWN.csv"
template_test <- gatingTemplate(gtFile)
template_test

flowDataPath <- "/Users/vishnuvundamati/Desktop/LLFSregate/LLFSflojoCOMPAR"
gs_outDir <- "/Users/vishnuvundamati/Desktop/LLFSregate/gs"
fcsFiles <- list.files(pattern = ".fcs$", flowDataPath, full = TRUE)
fcsFiles

nodesToHide <- c("CD19+","CD3+","CD19-","CD3-","CD8+","CD4+","CD4-CD8-","CD4+CD8-","CD4-CD8+",
                 "CD4+CD8+","CytotoxicT/CCR7+","CytotoxicT/CCR7-","CytotoxicT/CD45RA+","CytotoxicT/CD45RA-",
                 "HelperT/CCR7+","HelperT/CCR7-","HelperT/CD45RA+","HelperT/CD45RA-",
                 "CD27+","IgD+","IgD-","CD27-")

fcs_cdproblems <- vector()
fcs_gatingproblems <- vector()
C=0
D=0
for (f in fcsFiles) {
  
  print(basename(f))
  
  fs <- read.flowSet(f)
  compMat <- spillover(fs[[1]])[[1]]
  
  gs <- GatingSet(fs)
  gs
  print(paste0("gating set successful...",basename(f)))
  
  if(is_empty(markernames(gs))){
    C=C+1
    fcs_cdproblems[C] <- basename(f)
    print("CD names not found")
    next
  }
  
  # Compensation 
  gs <- compensate(gs,compMat)
  print(paste0("compensation complete...",basename(f)))
  
  # Transformation 
  chnls <- colnames(compMat)
  biexpTrans <-
    flowjo_biexp_trans(
      channelRange = 256,
      maxValue = 262144.0000291775 ,
      pos = 4.418539922,
      neg = 0,
      widthBasis = -100
    )
  tf <- transformerList(chnls, biexpTrans)
  gs <- flowCore::transform(gs,tf)
  print(paste0("biex transformation complete...",basename(f)))
  
  #applying the template to the gating set
  gs_call <- try(do.call("gt_gating",args=list(template_test,gs)))
  if (class(gs_call)=="try-error"){
    D=D+1
    fcs_gatingproblems[D] <- basename(f)
    print("gating unsuccessful")
    next
  }
  print(paste0("automated gating complete...",basename(f)))
  
  flowFrame_data <- gh_pop_get_data(gs[[1]], "/nonDebris/Lymphocytes/SingletCells/LiveCells")
  
  # NKT-like gate
  nkt_gate <- .gate_CD3_CD56_nkt(
    fr = flowFrame_data,
    pp_res = NULL,
    channels = c("APC-A", "BV510-A")
  )
  gs_pop_add(gs, nkt_gate, parent = "/nonDebris/Lymphocytes/SingletCells/LiveCells")
  recompute(gs)
  
  # Senescent T cell gate
  senescent_gate <- .gate_CD3_CD57_senescent(
    fr = flowFrame_data,
    pp_res = NULL,
    channels = c("APC-A", "PE-CF594-A")
  )
  gs_pop_add(gs, senescent_gate, parent = "/nonDebris/Lymphocytes/SingletCells/LiveCells")
  recompute(gs)
  
  
  
  lapply(nodesToHide, function(thisNode) gs_pop_set_visibility(gs, thisNode, FALSE))
  
  fcsname <- basename(f)
  outFile = paste0(gs_outDir,fcsname,"_gs")
  print(paste0("adding gated gating set file to..", outFile ))
  save_gs(gs,outFile)
  print(file.exists(outFile))
}

# 1. Lymphocytes: FSC-A vs SSC-A
p1 <- ggcyto(gs, subset = "/nonDebris", aes(x = "FSC-A", y = "SSC-A")) +
  geom_hex() + geom_gate("/nonDebris/Lymphocytes") + geom_stats() +
  ggtitle("Lymphocytes: FSC-A vs SSC-A")

# 2. Singlets: FSC-A vs FSC-H
p2 <- ggcyto(gs, subset = "/nonDebris/Lymphocytes", aes(x = "FSC-A", y = "FSC-H")) +
  geom_hex() + geom_gate("/nonDebris/Lymphocytes/SingletCells") + geom_stats() +
  ggtitle("Singlets: FSC-A vs FSC-H")

# 3. Live: L_D vs FSC-A
p3 <- ggcyto(gs, subset = "/nonDebris/Lymphocytes/SingletCells", aes(x = "BUV 496-A", y = "FSC-A")) + geom_hex() + 
  geom_gate("/nonDebris/Lymphocytes/SingletCells/LiveCells") + geom_stats() + ggtitle("Live: L_D vs FSC-A")

# 4. T Cell inclusion B Cell exclusion: CD19 vs CD3
p4 <- ggcyto(gs, subset = "/nonDebris/Lymphocytes/SingletCells/LiveCells", aes(x = "PE-Cy7-A", y = "APC-A")) + geom_hex() + 
  geom_gate("/nonDebris/Lymphocytes/SingletCells/LiveCells/Tcells") + geom_stats() + ggtitle("CD3+CD19-: CD3 vs CD19")

# 5. NKT-like: CD56 vs CD3
p5 <- ggcyto(gs, subset = "/nonDebris/Lymphocytes/SingletCells/LiveCells", aes(x = "BV510-A", y = "APC-A")) + geom_hex() +
  geom_gate("/nonDebris/Lymphocytes/SingletCells/LiveCells/NKTlike") + geom_stats() + ggtitle("NKT-like Cells: CD3 vs CD56")

# 6. Senescent T Cell: CD57 vs CD3
p6 <- ggcyto(gs, subset = "/nonDebris/Lymphocytes/SingletCells/LiveCells", aes(x = "PE-CF594-A", y = "APC-A")) + geom_hex() +
  geom_gate("/nonDebris/Lymphocytes/SingletCells/LiveCells/SenescentT") + geom_stats() + ggtitle("Senescent T Cells: CD3 vs CD57")


print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
