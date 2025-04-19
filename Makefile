nc = 12
rexec = R CMD BATCH --no-save --no-restore
rout = ./output/rout

.PHONY: all
all: sp

.PHONY: sp
sp: ./output/sprep_f1.RData ./output/sprep_norm.RData

## Format raw data for use in multidog
## Output: refmats and sizemats for three populations, and rowRanges object for site locations
./output/raw_counts.RData: ./code/format_for_multidog.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) $< $(rout)/$(basename $(<F)).Rout

## Genotype all three populations using multidog method, assuming model = "f1"
./output/updog_output_f1_p1.RData: ./code/multidog_genotyping_f1.R ./output/raw_counts.RData
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout

## Genotyping all three populations using multidog, assuming model = "norm"
./output/updog_output_norm_p1.RData: ./code/multidog_genotyping_norm.R ./output/raw_counts.RData
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout

## Prepare for segtest, polymapr, mappoly using norm data
./output/sprep_norm.RData: ./code/format_for_segtest_norm.R ./output/updog_output_norm_p1.RData
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) $< $(rout)/$(basename $(<F)).Rout

## Prepare for segtest, polymapr, mappoly using f1 data
./output/sprep_f1.RData: ./code/format_for_segtest_f1.R ./output/updog_output_f1_p1.RData
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) $< $(rout)/$(basename $(<F)).Rout

