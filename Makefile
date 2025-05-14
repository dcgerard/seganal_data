nc = 12
rexec = R CMD BATCH --no-save --no-restore
rout = ./output/rout

spout = ./output/figs/sp_pval_hist.pdf \
        ./output/figs/snps_mp_seg.pdf \
        ./output/figs/snps_mp_seg.txt \
        ./output/figs/snps_polymapr_seg.pdf \
        ./output/figs/snps_polymapr_seg.txt \
        ./output/figs/tab_polymapr.txt

.PHONY: all
all: sp

.PHONY: sp
sp: $(spout) ./output/segout_f1.RData ./output/segout_norm.RData ./output/segout_norm_competing.RData ./output/segout_f1_competing.RData

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

## Fit segtest on norm data
./output/segout_norm.RData: ./code/seg_fit_norm.R ./output/sprep_norm.RData
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout

## Fit segtest on f1 data
./output/segout_f1.RData: ./code/seg_fit_f1.R ./output/sprep_f1.RData
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout

## Fit segtest with competing models on norm data
./output/segout_norm_competing.RData: ./code/seg_fit_competing_norm.R ./output/sprep_norm.RData
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout

## Fit segtest with competing models on f1 data
./output/segout_f1_competing.RData: ./code/seg_fit_competing_f1.R ./output/sprep_f1.RData
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout

## Plot results
$(spout): ./code/plot_seg.R ./output/segout_f1.RData ./output/segout_f1_competing.RData ./output/sprep_f1.RData
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) $< $(rout)/$(basename $(<F)).Rout
