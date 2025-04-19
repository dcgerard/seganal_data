## Fit stuff to norm model
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

## Segtest with genotype likelihoods
plan(multisession, workers = nc)
segout_gl <- seg_multi(
  g = sprep_seg_gl$g,
  p1_ploidy = ploidy,
  p2_ploidy = ploidy,
  p1 = sprep_seg_gl$p1,
  p2 = sprep_seg_gl$p2,
  model = "seg"
  )
plan(sequential)

## polymapR with posterior probabilities
polymapr_out <- data.frame(
  p_value = rep(NA_real_, length.out = nrow(sprep_polymapr_pp$g)),
  p_invalid = rep(NA_real_, length.out = nrow(sprep_polymapr_pp$g))
)
for (i in seq_len(nrow(sprep_polymapr_pp$g))) {
  pnow <- polymapr_test(
    x = sprep_polymapr_pp$g[i, , ],
    g1 = sprep_polymapr_pp$p1[[i]],
    g2 = sprep_polymapr_pp$p2[[i]]
  )
  polymapr_out$p_value[[i]] <- pnow$p_value
  polymapr_out$p_invalid[[i]] <- pnow$p_invalid
}

## segtest with genotypes
plan(multisession, workers = nc)
segout_g <- seg_multi(
  g = sprep_seg_g$g,
  p1_ploidy = ploidy,
  p2_ploidy = ploidy,
  p1 = sprep_seg_g$p1,
  p2 = sprep_seg_g$p2,
  model = "seg"
)
plan(sequential)

## moppoly like test with genotypes
mapout_g <- seg_multi(
  g = sprep_mappoly_g$g,
  p1_ploidy = ploidy,
  p2_ploidy = ploidy,
  p1 = sprep_mappoly_g$p1,
  p2 = sprep_mappoly_g$p2,
  model = "auto",
  outlier = FALSE
)

## Save
save(
  segout_gl,
  segout_g,
  polymapr_out,
  mapout_g,
  file = "./output/segout_f1.RData")
