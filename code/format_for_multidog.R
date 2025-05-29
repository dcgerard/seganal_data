## Read them in
library(VariantAnnotation)
vd <- readVcf("./data/trifida_chr8.vcf")
nrow(vd)

## Filter out multiallelci SNPs
vd <- vd[sapply(fixed(vd)$ALT, length) < 2, ]
nrow(vd)

## Split parents and offspring
vd_parents <- vd[, grepl("^Beauregard\\_BT\\_\\d+", colnames(vd)) | grepl("^Tanzania\\_BT\\_\\d+", colnames(vd))]
vd <- vd[, grep("^BT", colnames(vd))]

refmat <- matrix(NA_real_, nrow = nrow(vd), ncol = ncol(vd))
sizemat <- matrix(NA_real_, nrow = nrow(vd), ncol = ncol(vd))
refmat_parents <- matrix(NA_real_, nrow = nrow(vd_parents), ncol = ncol(vd_parents))
sizemat_parents <- matrix(NA_real_, nrow = nrow(vd_parents), ncol = ncol(vd_parents))

## Do the offspring

ad <- geno(vd)$AD
dp <- geno(vd)$DP
for (i in seq_len(nrow(vd))) {
  for (j in seq_len(ncol(vd))) {
    refmat[i, j] <- ad[i, j][[1]][[1]]
    sizemat[i, j] <- sum(ad[i, j][[1]])
    stopifnot(sizemat[i, j] == dp[i, j])
  }
}
rownames(refmat) <- rownames(vd)
colnames(refmat) <- colnames(vd)
rownames(sizemat) <- rownames(vd)
colnames(sizemat) <- colnames(vd)

## Do the parents

ad_parents <- geno(vd_parents)$AD
dp_parents <- geno(vd_parents)$DP
for (i in seq_len(nrow(vd_parents))) {
  for (j in seq_len(ncol(vd_parents))) {
    refmat_parents[i, j] <- ad_parents[i, j][[1]][[1]]
    sizemat_parents[i, j] <- sum(ad_parents[i, j][[1]])
    stopifnot(sizemat_parents[i, j] == dp_parents[i, j])
  }
}
rownames(refmat_parents) <- rownames(vd_parents)
colnames(refmat_parents) <- colnames(vd_parents)
rownames(sizemat_parents) <- rownames(vd_parents)
colnames(sizemat_parents) <- colnames(vd_parents)

## Combine parents with offspring
refmat <- cbind(rowSums(refmat_parents[, grepl("^Beauregard\\_BT\\_\\d+", colnames(refmat_parents))]), rowSums(refmat_parents[, grepl("^Tanzania\\_BT\\_\\d+", colnames(refmat_parents))]), refmat)
sizemat <- cbind(rowSums(sizemat_parents[, grepl("^Beauregard\\_BT\\_\\d+", colnames(sizemat_parents))]), rowSums(sizemat_parents[, grepl("^Tanzania\\_BT\\_\\d+", colnames(sizemat_parents))]), sizemat)
colnames(refmat)[1:2] <- c("Beauregard_BT", "Tanzania_BT")
colnames(sizemat)[1:2] <- c("Beauregard_BT", "Tanzania_BT")

sites <- rowRanges(vd)

save(refmat, sizemat, sites, file = "./output/raw_counts.RData")
