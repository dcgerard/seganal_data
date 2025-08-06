
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Reproduce the data analysis of Gerard et al (2025)

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15784738.svg)](https://doi.org/10.5281/zenodo.15784738)
<!-- badges: end -->

This repo contains the analysis scripts needed to reproduce the
real-data analyses of Gerard et al. (2025). Data come from Mollinari et
al. (2020) via Figshare
([doi:10.25387/g3.10255844](https://doi.org/10.25387/g3.10255844)).
These analyses were evaluated using R version 4.5.0. The exact versions
of all packages used are listed in the “renv.lock” file.

1.  Download the data from Figshare:
    <https://gsajournals.figshare.com/ndownloader/files/18517700>

2.  Rename the downloaded file “trifida_chr8.vcf”.

3.  Place “trifida_chr8.vcf” in the “data” folder.

4.  Edit the `nc` argument in the `Makefile` to be the number of cores
    you wish to use for the analysis.

5.  Open up R and restore the `renv`:

    ``` r
    renv::restore()
    ```

6.  Run `make` in the terminal:

    ``` bash
    make
    ```

7.  Get coffee:

    - [Aslin Coffee DC](https://maps.app.goo.gl/n8vVbjkwwrC9fiyy5)
    - [Doubles](https://maps.app.goo.gl/CXNaN1HpgVxZDk9h6)
    - [Bar Americano](https://maps.app.goo.gl/U6XJmTazJssadUS4A)
    - [BREATHE CO](https://maps.app.goo.gl/CpVTvioWjSbm8zWx5)
    - [Sad Coffee Co.](https://maps.app.goo.gl/KYKTVSi57dWizNTQA)

# References

Gerard, D, Ambrosano, GB, Pereira, GdS, & Garcia, AAF (2025). “Tests for
segregation distortion in higher ploidy F1 populations.” *bioRxiv*,
p. 1–20.
[bioRxiv:2025.06.23.661114](https://doi.org/10.1101/2025.06.23.661114)

Mollinari, M., Olukolu, B. A., Pereira, G. D. S., Khan, A., Gemenet, D.,
Yencho, G. C., & Zeng, Z. B. (2020). Unraveling the hexaploid
sweetpotato inheritance using ultra-dense multilocus mapping. *G3:
Genes, Genomes, Genetics*, 10(1), 281-292.
[doi:10.1534/g3.119.400620](https://doi.org/10.1534/g3.119.400620)
