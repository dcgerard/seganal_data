## Fit competing models to norm model output
library(segtest)
library(future)
load("./output/sprep_f1.RData")
ploidy <- 6

## Set up cluster
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  nc <- 1
} else {
  eval(parse(text = args[[1]]))
}

## autodr
plan(multisession, workers = nc)
segout_gl_autodr <- seg_multi(
  g = sprep_seg_gl$g,
  p1_ploidy = ploidy,
  p2_ploidy = ploidy,
  p1 = sprep_seg_gl$p1,
  p2 = sprep_seg_gl$p2,
  model = "auto_dr"
  )
plan(sequential)

plan(multisession, workers = nc)
segout_gl_auto <- seg_multi(
  g = sprep_seg_gl$g,
  p1_ploidy = ploidy,
  p2_ploidy = ploidy,
  p1 = sprep_seg_gl$p1,
  p2 = sprep_seg_gl$p2,
  model = "auto"
)
plan(sequential)

## Save
save(
  segout_gl_auto,
  segout_gl_autodr,
  file = "./output/segout_f1_competing.RData")
