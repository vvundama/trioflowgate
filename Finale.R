library(flowWorkspace)
library(ggcyto)
library(ggplot2)

gs_path <- "/Users/vishnuvundamati/Desktop...REPLACE WITH DIRECTORY CONTAINING APPENDED GATING SETS"
gs <- load_gs(gs_path)

sampleNames(gs)
gh <- gs[[1]]
getNodes(gh)
autoplot(gs[[1]], "NKTlike")
autoplot(gs[[1]], "SenescentT")
autoplot(gs[[1]], "DoublePositive")

