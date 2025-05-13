####################
## Analyze output of tests
####################
library(tidyverse)
library(updog)
library(segtest)
library(xtable)
load("./output/segout_f1.RData")
load("./output/segout_f1_competing.RData")
load("./output/sprep_f1.RData")

## Compare segtest models ----
phigh <- segout_gl$p_value > 0.1
x <- sum(segout_gl$null_bic[phigh] - segout_gl_auto$null_bic[phigh] > 0)
n <- sum(phigh)
binom.test(x = x, n = n, p = 0.5)

x <- sum(segout_gl$null_bic[phigh] - segout_gl_autodr$null_bic[phigh] > 0)
n <- sum(phigh)
binom.test(x = x, n = n, p = 0.5)

phigh <- segout_gl_autodr$p_value > 0.1
x <- sum(segout_gl_autodr$null_bic[phigh] - segout_gl_auto$null_bic[phigh] < 0)
n <- sum(phigh)
binom.test(x = x, n = n, p = 0.5)

## p-value plots ----
select(mapout_g, snp, mappoly = p_value) |>
  full_join(select(segout_gl_autodr, snp, segtest = p_value), by = join_by(snp)) |>
  full_join(select(polymapr_out, snp, polymapr = p_value), by = join_by(snp)) |>
  left_join(uout_sub_gl$snpdf, by = join_by(snp)) ->
  df_p

df_p |>
  select(snp, mappoly, segtest, polymapr) |>
  pivot_longer(cols = -snp, names_to = "method", values_to = "p_value") |>
  ggplot(aes(x = p_value)) +
  geom_histogram(bins = 30, color = "black", fill = "white") +
  facet_wrap(.~method, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) ->
  pl

ggsave(
  filename = "./output/figs/sp_pval_hist.pdf",
  plot = pl,
  height = 3,
  width = 6,
  family = "Times")

alpha <- 0.05
sum(p.adjust(df_p$mappoly, method = "BH") > alpha, na.rm = TRUE)
sum(p.adjust(df_p$polymapr, method = "BH") > alpha, na.rm = TRUE)
sum(p.adjust(df_p$segtest, method = "BH") > alpha, na.rm = TRUE)

df_p |>
  filter(!is.na(mappoly), !is.na(segtest), !is.na(polymapr)) |>
  select(snp, mappoly, segtest, polymapr) |>
  pivot_longer(cols = -snp, names_to = "method", values_to = "p_value") |>
  ggplot(aes(x = p_value)) +
  geom_histogram(bins = 30, color = "black", fill = "white") +
  facet_wrap(.~method) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"))

cor(df_p[, c("mappoly", "polymapr", "segtest")], use = "pairwise.complete.obs")

## mappoly and segtest differ ----
df_p |>
  filter(mappoly < 0.01, segtest > 0.5)

df_p |>
  filter(snp %in% c("S8_2301244", "S8_2301256", "S8_16038595", "S8_4395910")) ->
  df_m_low_s_high

df_p |>
  filter(snp %in% c("S8_4866982", "S8_18370562", "S8_460771", "S8_2212186")) ->
  df_s_low_m_high

uout_ps <- filter_snp(uout_sub_g_mp, snp %in% df_m_low_s_high$snp | snp %in% df_s_low_m_high$snp)

snplist <- c(df_m_low_s_high$snp, df_s_low_m_high$snp)

uout_ps$snpdf |>
  select(snp, bias, seq) |>
  expand_grid(g = 0:6) |>
  mutate(xi = updog:::xi_fun(p = g / 6, eps = seq, h = bias),
         slope = xi / (1 - xi),
         intercept = 0,
         snp = factor(snp, levels = snplist)) ->
  slope_df


uout_ps$inddf |>
  mutate(
    alt = size - ref,
    geno = factor(geno, levels = 0:6),
    snp = factor(snp, levels = snplist)) ->
  dat_df

dat_df |>
  group_by(snp) |>
  summarize(max_ref = max(ref), max_alt = max(alt)) |>
  mutate(ref = pmax(max_ref, max_alt),
         alt = ref) ->
  lim_df

colvec <- rev(palette.colors(n = 7))
ggplot() +
  geom_point(data = dat_df, mapping = aes(x = alt, y = ref, color = geno)) +
  xlab("Alternative Count") +
  ylab("Reference Count") +
  facet_wrap(.~ snp, scales = "free", ncol = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  scale_color_manual(values = colvec, name = "Posterior Mode\nGenotype") +
  geom_abline(data = slope_df, mapping = aes(slope = slope, intercept = intercept, group = g), linetype = "dashed", color = "grey") +
  geom_point(data = lim_df, aes(x = alt, y = ref), alpha = 0) ->
  pl

ggsave(filename = "./output/figs/snps_mp_seg.pdf", plot = pl, height = 7, width = 5, family = "Times")

## rerun mappoly approach with all genotypes
uout_m_low_s_high <- filter_snp(uout_sub_gl, snp %in% df_m_low_s_high$snp)
sprep_m_low_s_high <- multidog_to_g(mout = uout_m_low_s_high, type = "off_g", ploidy = 6)
sout_m_low_s_high <- seg_multi(g = sprep_m_low_s_high$g, p1_ploidy = 6, p2_ploidy = 6, p1 = sprep_m_low_s_high$p1, p2 = sprep_m_low_s_high$p2, model = "auto", outlier = FALSE)

df_m_low_s_high |>
  select(snp, mappoly, segtest) |>
  left_join(select(sout_m_low_s_high, snp, mappoly_new = p_value)) ->
  tab_m_low_s_high

## rerun segtest approach without low counts
uout_s_low_m_high <- filter_snp(uout_sub_gl, snp %in% df_s_low_m_high$snp)
uout_s_low_m_high$inddf |>
  filter(ref > 30) ->
  uout_s_low_m_high$inddf
sprep_s_low_m_high <- multidog_to_g(mout = uout_s_low_m_high, type = "off_gl", ploidy = 6)
sout_s_low_m_high <- seg_multi(g = sprep_s_low_m_high$g, p1_ploidy = 6, p2_ploidy = 6, p1 = sprep_s_low_m_high$p1, p2 = sprep_s_low_m_high$p2, model = "auto", outlier = FALSE)

df_s_low_m_high |>
  select(snp, mappoly, segtest) |>
  left_join(select(sout_s_low_m_high, snp, segtest_new = p_value)) ->
  tab_s_low_m_high

bind_rows(tab_m_low_s_high, tab_s_low_m_high) |>
  xtable(display = rep("G", 6), label = "tab:snps.mp.seg") |>
  print(include.rownames = FALSE, file = "./output/figs/snps_mp_seg.txt")


table(filter_snp(uout_sub_gl, snp == "S8_2212186")$inddf$geno)

## polymapR and segtest differ ----
df_p |>
  filter(snp %in% c("S8_18370562", "S8_4866982", "S8_2212186", "S8_17348718")) ->
  df_s_low_p_high

df_p |>
  filter(snp %in% c("S8_15463274", "S8_4572598", "S8_41319", "S8_10014838")) ->
  df_p_low_s_high
