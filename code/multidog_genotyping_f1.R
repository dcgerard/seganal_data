library(updog)
load("./output/raw_counts.RData")
ploidy <- 6

## Set up cluster
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  nc <- 1
} else {
  eval(parse(text = args[[1]]))
}

## genotype first population
uout_f1_p1 <- multidog(
  refmat = refmat_p1,
  sizemat = sizemat_p1,
  ploidy = ploidy,
  model = "f1",
  p1_id = "Beauregard_BT",
  p2_id = "Tanzania_BT",
  nc = nc
)

save(uout_f1_p1, file = "./output/updog_output_f1_p1.RData")
