## Configure the f1 model for input into segtest
library(updog)
library(segtest)
library(tidyverse)

load("./output/updog_output_f1_p1.RData")
ploidy <- 6

## Filter threshholds ----------------------------------------------------------
## Average read-depth. Same as in Mollinari et al (2020)
##     (https://doi.org/10.1534/g3.119.400620) p282
bound_rd <- 20

## Number of missing individuals bound. Same as in Mollinari et al (2020)
##     (https://doi.org/10.1534/g3.119.400620) p282
bound_miss <- 0.25

## Maxpostprob bound. Same as default in mappoly::read_geno_prob()
bound_pp <- 0.95

## Simple filter for genotype likelihoods --------------------------------------
## Monomorphic filter
uout_f1_p1$snpdf|>
  filter(!(p1geno %in% c(0, ploidy) & p2geno %in% c(0, ploidy))) ->
  filter_df_1

## Read depth filter
uout_f1_p1$inddf |>
  group_by(snp) |>
  summarize(ave_rd = mean(size)) |>
  filter(ave_rd > bound_rd) ->
  filter_df_2

## Subset
snp_list <- intersect(filter_df_1$snp, filter_df_2$snp)
uout_sub_gl <- filter_snp(x = uout_f1_p1, expr = snp %in% snp_list)

## prepare for segtest
sprep_seg_gl <- multidog_to_g(
  mout = uout_sub_gl,
  type = "off_gl",
  ploidy = ploidy)
names(sprep_seg_gl$p1) <- uout_sub_gl$snpdf$snp
names(sprep_seg_gl$p2) <- uout_sub_gl$snpdf$snp
stopifnot(sprep_seg_gl$p1 == uout_sub_gl$snpdf$p1geno)
stopifnot(sprep_seg_gl$p2 == uout_sub_gl$snpdf$p2geno)

## prepare for polymapR
sprep_polymapr_pp <- list()
sprep_polymapr_pp$g <- format_multidog(x = uout_sub_gl, varname = paste0("Pr_", 0:ploidy))
sprep_polymapr_pp$p1 <- uout_sub_gl$snpdf$p1geno
names(sprep_polymapr_pp$p1) <- uout_sub_gl$snpdf$snp
sprep_polymapr_pp$p2 <- uout_sub_gl$snpdf$p2geno
names(sprep_polymapr_pp$p2) <- uout_sub_gl$snpdf$snp
stopifnot(names(sprep_polymapr_pp$p1) == rownames(sprep_polymapr_pp$g))
stopifnot(names(sprep_polymapr_pp$p2) == rownames(sprep_polymapr_pp$g))

## More extensive filter for known genotypes -----------------------------------
## Monomorphic filter
uout_f1_p1$snpdf|>
  filter(!(p1geno %in% c(0, ploidy) & p2geno %in% c(0, ploidy))) ->
  filter_df_1

## maxpostprob filter, missing data filter, average read-depth filter
uout_f1_p1$inddf |>
  mutate(size = if_else(maxpostprob < bound_pp, NA, size)) |>
  group_by(snp) |>
  summarize(ave_rd = mean(size, na.rm = TRUE), pna = mean(is.na(size))) |>
  filter(ave_rd > bound_rd, pna < bound_miss) ->
  filter_df_2

snp_list <- intersect(filter_df_1$snp, filter_df_2$snp)
uout_sub_g <- filter_snp(x = uout_f1_p1, expr = snp %in% snp_list)

uout_sub_g$inddf$geno[uout_sub_g$inddf$maxpostprob < bound_pp] <- NA

## Prepare for segtest
sprep_seg_g <- multidog_to_g(
  mout = uout_sub_g,
  type = "off_g",
  ploidy = ploidy)

## Prepare for mappoly
## Filter out "invalid genotypes"
uout_sub_g_mp <- uout_sub_g
TOL <- sqrt(.Machine$double.eps)
ninvalid <- rep(NA_real_, length.out = nrow(uout_sub_g_mp$snpdf))
names(ninvalid) <- uout_sub_g_mp$snpdf$snp
for (i in seq_len(nrow(uout_sub_g_mp$snpdf))) {
  snp <- uout_sub_g_mp$snpdf$snp[[i]]
  p1_geno <- uout_sub_g_mp$snpdf$p1geno[[i]]
  p2_geno <- uout_sub_g_mp$snpdf$p2geno[[i]]
  gf <- segtest::gf_freq(
    p1_g = p1_geno,
    p1_ploidy = ploidy,
    p1_alpha = 0,
    p1_type = "polysomic",
    p2_g = p2_geno,
    p2_ploidy = ploidy,
    p2_alpha = 0,
    p2_type = "polysomic",
    pi = 0,
    nudge = 0)
  invalid_genos <- (0:ploidy)[gf < TOL]
  badind <- uout_sub_g_mp$inddf$snp == snp & uout_sub_g_mp$inddf$geno %in% invalid_genos
  ninvalid[[i]] <- sum(badind)
  if (ninvalid[[i]] > 1) {
    uout_sub_g_mp$inddf[badind, ]$geno <- NA
  }
}

sprep_mappoly_g <- multidog_to_g(
  mout = uout_sub_g_mp,
  type = "off_g",
  ploidy = ploidy)
sprep_mappoly_g$ninvalid <- ninvalid

## Save results
save(
  uout_sub_gl,
  uout_sub_g,
  uout_sub_g_mp,
  sprep_seg_gl,
  sprep_polymapr_pp,
  sprep_seg_g,
  sprep_mappoly_g,
  file = "./output/sprep_f1.RData")
