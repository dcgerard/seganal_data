####################
## Analyze output of tests
####################
library(tidyverse)
library(GGally)
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
  theme(strip.background = element_rect(fill = "white"))

df_p |>
  filter(!is.na(mappoly), !is.na(segtest), !is.na(polymapr)) |>
  select(snp, mappoly, segtest, polymapr) |>
  pivot_longer(cols = -snp, names_to = "method", values_to = "p_value") |>
  ggplot(aes(x = p_value)) +
  geom_histogram(bins = 30, color = "black", fill = "white") +
  facet_wrap(.~method) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"))

ggpairs(data = select(df_p, mappoly, polymapr, segtest), mapping = aes(alpha = I(1 / 10)))

