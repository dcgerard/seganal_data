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
doParallel::registerDoParallel()
uout_norm_p1 <- multidog(
  refmat = refmat,
  sizemat = sizemat,
  ploidy = ploidy,
  model = "norm",
  nc = nc
)

save(uout_norm_p1, file = "./output/updog_output_norm_p1.RData")
