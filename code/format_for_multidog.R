library(VariantAnnotation)
vdat <- readVcf("./data/myGBSGenos_trifida_mergedSNPs_mergedTaxa_chr8.vcf")
p1 <- vdat[, grep("^BT", colnames(vdat))]
p2 <- vdat[, grep("^TB", colnames(vdat))]
p3 <- vdat[, grep("^NKB", colnames(vdat))]

p1_parents <- vdat[, grepl("^Beauregard\\_BT\\_\\d+", colnames(vdat)) | grepl("^Tanzania\\_BT\\_\\d+", colnames(vdat))]
p2_parents <- vdat[, grepl("^Tanzania\\_TB\\_\\d+", colnames(vdat)) | grepl("^Beauregard\\_TB\\_\\d+", colnames(vdat))]
p3_parents <- vdat[, grepl("^NewKawogo\\_NKB\\_\\d+", colnames(vdat)) | grepl("^Beauregard\\_NKB\\_\\d+", colnames(vdat))]

any(colnames(p1) %in% colnames(p2))
any(colnames(p1) %in% colnames(p3))
any(colnames(p2) %in% colnames(p3))

refmat_p1 <- matrix(NA_real_, nrow = nrow(p1), ncol = ncol(p1))
sizemat_p1 <- matrix(NA_real_, nrow = nrow(p1), ncol = ncol(p1))
refmat_p2 <- matrix(NA_real_, nrow = nrow(p2), ncol = ncol(p2))
sizemat_p2 <- matrix(NA_real_, nrow = nrow(p2), ncol = ncol(p2))
refmat_p3 <- matrix(NA_real_, nrow = nrow(p3), ncol = ncol(p3))
sizemat_p3 <- matrix(NA_real_, nrow = nrow(p3), ncol = ncol(p3))

refmat_p1_parents <- matrix(NA_real_, nrow = nrow(p1_parents), ncol = ncol(p1_parents))
sizemat_p1_parents <- matrix(NA_real_, nrow = nrow(p1_parents), ncol = ncol(p1_parents))
refmat_p2_parents <- matrix(NA_real_, nrow = nrow(p2_parents), ncol = ncol(p2_parents))
sizemat_p2_parents <- matrix(NA_real_, nrow = nrow(p2_parents), ncol = ncol(p2_parents))
refmat_p3_parents <- matrix(NA_real_, nrow = nrow(p3_parents), ncol = ncol(p3_parents))
sizemat_p3_parents <- matrix(NA_real_, nrow = nrow(p3_parents), ncol = ncol(p3_parents))

## Do the offspring

ad_p1 <- geno(p1)$AD
dp_p1 <- geno(p1)$DP
for (i in seq_len(nrow(p1))) {
  for (j in seq_len(ncol(p1))) {
    refmat_p1[i, j] <- ad_p1[i, j][[1]][[1]]
    sizemat_p1[i, j] <- sum(ad_p1[i, j][[1]])
    stopifnot(sizemat_p1[i, j] == dp_p1[i, j])
  }
}
rownames(refmat_p1) <- rownames(p1)
colnames(refmat_p1) <- colnames(p1)
rownames(sizemat_p1) <- rownames(p1)
colnames(sizemat_p1) <- colnames(p1)

ad_p2 <- geno(p2)$AD
dp_p2 <- geno(p2)$DP
for (i in seq_len(nrow(p2))) {
  for (j in seq_len(ncol(p2))) {
    refmat_p2[i, j] <- ad_p2[i, j][[1]][[1]]
    sizemat_p2[i, j] <- sum(ad_p2[i, j][[1]])
    stopifnot(sizemat_p2[i, j] == dp_p2[i, j])
  }
}
rownames(refmat_p2) <- rownames(p2)
colnames(refmat_p2) <- colnames(p2)
rownames(sizemat_p2) <- rownames(p2)
colnames(sizemat_p2) <- colnames(p2)

ad_p3 <- geno(p3)$AD
dp_p3 <- geno(p3)$DP
for (i in seq_len(nrow(p3))) {
  for (j in seq_len(ncol(p3))) {
    refmat_p3[i, j] <- ad_p3[i, j][[1]][[1]]
    sizemat_p3[i, j] <- sum(ad_p3[i, j][[1]])
    stopifnot(sizemat_p3[i, j] == dp_p3[i, j])
  }
}
rownames(refmat_p3) <- rownames(p3)
colnames(refmat_p3) <- colnames(p3)
rownames(sizemat_p3) <- rownames(p3)
colnames(sizemat_p3) <- colnames(p3)

## Do the parents

ad_p1_parents <- geno(p1_parents)$AD
dp_p1_parents <- geno(p1_parents)$DP
for (i in seq_len(nrow(p1_parents))) {
  for (j in seq_len(ncol(p1_parents))) {
    refmat_p1_parents[i, j] <- ad_p1_parents[i, j][[1]][[1]]
    sizemat_p1_parents[i, j] <- sum(ad_p1_parents[i, j][[1]])
    stopifnot(sizemat_p1_parents[i, j] == dp_p1_parents[i, j])
  }
}
rownames(refmat_p1_parents) <- rownames(p1_parents)
colnames(refmat_p1_parents) <- colnames(p1_parents)
rownames(sizemat_p1_parents) <- rownames(p1_parents)
colnames(sizemat_p1_parents) <- colnames(p1_parents)

ad_p2_parents <- geno(p2_parents)$AD
dp_p2_parents <- geno(p2_parents)$DP
for (i in seq_len(nrow(p2_parents))) {
  for (j in seq_len(ncol(p2_parents))) {
    refmat_p2_parents[i, j] <- ad_p2_parents[i, j][[1]][[1]]
    sizemat_p2_parents[i, j] <- sum(ad_p2_parents[i, j][[1]])
    stopifnot(sizemat_p2_parents[i, j] == dp_p2_parents[i, j])
  }
}
rownames(refmat_p2_parents) <- rownames(p2_parents)
colnames(refmat_p2_parents) <- colnames(p2_parents)
rownames(sizemat_p2_parents) <- rownames(p2_parents)
colnames(sizemat_p2_parents) <- colnames(p2_parents)

ad_p3_parents <- geno(p3_parents)$AD
dp_p3_parents <- geno(p3_parents)$DP
for (i in seq_len(nrow(p3_parents))) {
  for (j in seq_len(ncol(p3_parents))) {
    refmat_p3_parents[i, j] <- ad_p3_parents[i, j][[1]][[1]]
    sizemat_p3_parents[i, j] <- sum(ad_p3_parents[i, j][[1]])
    stopifnot(sizemat_p3_parents[i, j] == dp_p3_parents[i, j])
  }
}
rownames(refmat_p3_parents) <- rownames(p3_parents)
colnames(refmat_p3_parents) <- colnames(p3_parents)
rownames(sizemat_p3_parents) <- rownames(p3_parents)
colnames(sizemat_p3_parents) <- colnames(p3_parents)

## Combine parents with offspring
refmat_p1 <- cbind(rowSums(refmat_p1_parents[, grepl("^Beauregard\\_BT\\_\\d+", colnames(refmat_p1_parents))]), rowSums(refmat_p1_parents[, grepl("^Tanzania\\_BT\\_\\d+", colnames(refmat_p1_parents))]), refmat_p1)
sizemat_p1 <- cbind(rowSums(sizemat_p1_parents[, grepl("^Beauregard\\_BT\\_\\d+", colnames(sizemat_p1_parents))]), rowSums(sizemat_p1_parents[, grepl("^Tanzania\\_BT\\_\\d+", colnames(sizemat_p1_parents))]), sizemat_p1)
colnames(refmat_p1)[1:2] <- c("Beauregard_BT", "Tanzania_BT")
colnames(sizemat_p1)[1:2] <- c("Beauregard_BT", "Tanzania_BT")

refmat_p2 <- cbind(rowSums(refmat_p2_parents[, grepl("^Tanzania\\_TB\\_\\d+", colnames(refmat_p2_parents))]), rowSums(refmat_p2_parents[, grepl("^Beauregard\\_TB\\_\\d+", colnames(refmat_p2_parents))]), refmat_p2)
sizemat_p2 <- cbind(rowSums(sizemat_p2_parents[, grepl("^Tanzania\\_TB\\_\\d+", colnames(sizemat_p2_parents))]), rowSums(sizemat_p2_parents[, grepl("^Beauregard\\_TB\\_\\d+", colnames(sizemat_p2_parents))]), sizemat_p2)
colnames(refmat_p2)[1:2] <- c("Tanzania_TB", "Beauregard_TB")
colnames(sizemat_p2)[1:2] <- c("Tanzania_TB", "Beauregard_TB")

refmat_p3 <- cbind(rowSums(refmat_p3_parents[, grepl("^NewKawogo\\_NKB\\_\\d+", colnames(refmat_p3_parents))]), rowSums(refmat_p3_parents[, grepl("^Beauregard\\_NKB\\_\\d+", colnames(refmat_p3_parents))]), refmat_p3)
sizemat_p3 <- cbind(rowSums(sizemat_p3_parents[, grepl("^NewKawogo\\_NKB\\_\\d+", colnames(sizemat_p3_parents))]), rowSums(sizemat_p3_parents[, grepl("^Beauregard\\_NKB\\_\\d+", colnames(sizemat_p3_parents))]), sizemat_p3)
colnames(refmat_p3)[1:2] <- c("NewKawogo_NKB", "Beauregard_NKB")
colnames(sizemat_p3)[1:2] <- c("NewKawogo_NKB", "Beauregard_NKB")

sites <- rowRanges(vdat)

save(refmat_p1, refmat_p2, refmat_p3, sizemat_p1, sizemat_p2, sizemat_p3, sites, file = "./output/raw_counts.RData")
