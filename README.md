
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Reproduce the real-data analysis of Gerard et al (2025)

<!-- badges: start -->

<!-- badges: end -->

This repo contains the analysis scripts needed to reproduce the
real-data analyses of Gerard et al. (2025).

1.  Open up R and restore the `renv`:

        renv::restore()

2.  Run make in the terminal.

        make

3.  Get coffee.

    - [Aslin Coffee DC](https://maps.app.goo.gl/n8vVbjkwwrC9fiyy5)
    - [Doubles](https://maps.app.goo.gl/CXNaN1HpgVxZDk9h6)
    - [Bar Americano](https://maps.app.goo.gl/U6XJmTazJssadUS4A)

Note that running `make` uses `wget` to download the raw data from
Figshare (<https://doi.org/10.25387/g3.10255844>). You can instead
download it manually
(<https://gsajournals.figshare.com/ndownloader/files/18517700>). Just
make sure the VCF file is called “trifida_chr8.vcf” and is placed in the
“data” folder.

# References

Gerard D, Ambrosano GB, Pereira GDS, & Garcia AAF (2025). “Tests for
segregation distortion in higher ploidy F1 populations.” *Unpublished
Manuscript*.
