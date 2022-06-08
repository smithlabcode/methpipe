# This file is part of methpipe
#
# Copyright (C) 2010-2019: Andrew D. Smith
#
# Authors: Andrew D. Smith
#
# This is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

EXTRA_DIST = README.md
ACLOCAL_AMFLAGS = -I m4

SUBDIRS := src/smithlab_cpp
install installdirs: SUBDIRS := $(filter-out src/smithlab_cpp, $(SUBDIRS))
AM_CPPFLAGS = -I $(top_srcdir)/src/common -I $(top_srcdir)/src/smithlab_cpp

CXXFLAGS = -O3 # default has optimization on

noinst_LIBRARIES = libmethpipe.a
libmethpipe_a_SOURCES = \
	src/common/BetaBin.cpp \
	src/common/Distro.cpp \
	src/common/EmissionDistribution.cpp \
	src/common/Epiread.cpp \
	src/common/EpireadStats.cpp \
	src/common/LevelsCounter.cpp \
	src/common/MethpipeSite.cpp \
	src/common/Smoothing.cpp \
	src/common/ThreeStateHMM.cpp \
	src/common/TwoStateHMM.cpp \
	src/common/TwoStateHMM_PMD.cpp \
	src/common/bsutils.cpp \
	src/common/numerical_utils.cpp

libmethpipe_a_SOURCES += \
	src/common/BetaBin.hpp \
	src/common/Distro.hpp \
	src/common/EmissionDistribution.hpp \
	src/common/Epiread.hpp \
	src/common/EpireadStats.hpp \
	src/common/LevelsCounter.hpp \
	src/common/MethpipeFiles.hpp \
	src/common/MethpipeSite.hpp \
	src/common/Smoothing.hpp \
	src/common/ThreeStateHMM.hpp \
	src/common/TwoStateHMM.hpp \
	src/common/TwoStateHMM_PMD.hpp \
	src/common/bsutils.hpp \
	src/common/numerical_utils.hpp

LDADD = libmethpipe.a src/smithlab_cpp/libsmithlab_cpp.a

### ANALYSIS SUBDIR
bin_PROGRAMS = pmd methcounts bsrate hmr hypermr levels roimethstat \
	methstates methentropy hmr_rep multimethstat

pmd_SOURCES = src/analysis/pmd.cpp
methstates_SOURCES =  src/analysis/methstates.cpp
bsrate_SOURCES = src/analysis/bsrate.cpp
methentropy_SOURCES = src/analysis/methentropy.cpp
methcounts_SOURCES = src/analysis/methcounts.cpp
roimethstat_SOURCES = src/analysis/roimethstat.cpp
multimethstat_SOURCES = src/analysis/multimethstat.cpp
hmr_SOURCES = src/analysis/hmr.cpp
hmr_rep_SOURCES = src/analysis/hmr_rep.cpp
levels_SOURCES = src/analysis/levels.cpp
hypermr_SOURCES = src/analysis/hypermr.cpp

### UTILS SUBDIR
bin_PROGRAMS += lc_approx fast-liftover lift-filter merge-bsrate \
	merge-methcounts duplicate-remover symmetric-cpgs \
	clean-hairpins selectsites guessprotocol format_reads

clean_hairpins_SOURCES = src/utils/clean-hairpins.cpp
guessprotocol_SOURCES = src/utils/guessprotocol.cpp
duplicate_remover_SOURCES = src/utils/duplicate-remover.cpp
merge_bsrate_SOURCES = src/utils/merge-bsrate.cpp
format_reads_SOURCES = src/utils/format_reads.cpp
lc_approx_SOURCES = src/utils/lc_approx.cpp
selectsites_SOURCES = src/utils/selectsites.cpp
symmetric_cpgs_SOURCES = src/utils/symmetric-cpgs.cpp
merge_methcounts_SOURCES = src/utils/merge-methcounts.cpp
lift_filter_SOURCES = src/utils/lift-filter.cpp
fast_liftover_SOURCES = src/utils/fast-liftover.cpp

### AMRFINDER SUBDIR
bin_PROGRAMS += allelicmeth amrfinder amrtester

allelicmeth_SOURCES = src/amrfinder/allelicmeth.cpp
amrfinder_SOURCES = src/amrfinder/amrfinder.cpp
amrtester_SOURCES = src/amrfinder/amrtester.cpp

### RADMETH SUBDIR
bin_PROGRAMS += radmeth methdiff dmr

dmr_SOURCES = src/radmeth/dmr.cpp
methdiff_SOURCES = src/radmeth/methdiff.cpp
radmeth_SOURCES = \
	src/radmeth/radmeth.cpp \
	src/radmeth/regression.cpp \
	src/radmeth/combine_pvals.cpp
radmeth_SOURCES += \
	src/radmeth/regression.hpp \
	src/radmeth/combine_pvals.hpp \
	src/radmeth/merge.hpp

### RSEG SUBDIR
bin_PROGRAMS += rseg rseg-diff deadzones

rseg_SOURCES = \
	src/rseg/rseg.cpp \
	src/rseg/numerical_utils.cpp \
	src/rseg/rseg_utils.cpp \
	src/rseg/Distro.cpp \
	src/rseg/SplitDistro.cpp \
	src/rseg/Smoothing.cpp \
	src/rseg/BAMFile.cpp \
	src/rseg/SortGenomicRegion.cpp \
	src/rseg/LoadReadsByRegion.cpp \
	src/rseg/SelectBinSize.cpp \
	src/rseg/ReadCounts.cpp \
	src/rseg/TwoStateScaleHMM.cpp \
	src/rseg/TwoStateScaleResolveMixture.cpp \
	src/rseg/EvaluateBoundaries.cpp

rseg_SOURCES += \
	src/rseg/numerical_utils.hpp \
	src/rseg/rseg_utils.hpp \
	src/rseg/Distro.hpp \
	src/rseg/SplitDistro.hpp \
	src/rseg/Smoothing.hpp \
	src/rseg/BAMFile.hpp \
	src/rseg/SortGenomicRegion.hpp \
	src/rseg/LoadReadsByRegion.hpp \
	src/rseg/SelectBinSize.hpp \
	src/rseg/ReadCounts.hpp \
	src/rseg/TwoStateScaleHMM.hpp \
	src/rseg/TwoStateScaleResolveMixture.hpp \
	src/rseg/EvaluateBoundaries.hpp

rseg_diff_SOURCES = \
	src/rseg/rseg-diff.cpp \
	src/rseg/numerical_utils.cpp \
	src/rseg/rseg_utils.cpp \
	src/rseg/Distro.cpp \
	src/rseg/SplitDistro.cpp \
	src/rseg/Smoothing.cpp \
	src/rseg/BAMFile.cpp \
	src/rseg/SortGenomicRegion.cpp \
	src/rseg/LoadReadsByRegion.cpp \
	src/rseg/SelectBinSize.cpp \
	src/rseg/ReadCounts.cpp \
	src/rseg/TwoStateScaleSplitHMM.cpp \
	src/rseg/TwoStateScaleSplitResolveMixture.cpp \
	src/rseg/ThreeStateScaleSplitHMM.cpp \
	src/rseg/ThreeStateScaleSplitResolveMixture.cpp \
	src/rseg/EvaluateBoundaries.cpp

rseg_diff_SOURCES += \
	src/rseg/numerical_utils.hpp \
	src/rseg/rseg_utils.hpp \
	src/rseg/Distro.hpp \
	src/rseg/SplitDistro.hpp \
	src/rseg/Smoothing.hpp \
	src/rseg/BAMFile.hpp \
	src/rseg/SortGenomicRegion.hpp \
	src/rseg/LoadReadsByRegion.hpp \
	src/rseg/SelectBinSize.hpp \
	src/rseg/ReadCounts.hpp \
	src/rseg/TwoStateScaleSplitHMM.hpp \
	src/rseg/TwoStateScaleSplitResolveMixture.hpp \
	src/rseg/ThreeStateScaleSplitHMM.hpp \
	src/rseg/ThreeStateScaleSplitResolveMixture.hpp \
	src/rseg/EvaluateBoundaries.hpp

deadzones_SOURCES = src/rseg/deadzones.cpp

### MLML SUBDIR
bin_PROGRAMS += mlml

mlml_SOURCES = src/mlml/mlml.cpp
