# Simon et al 2025

This repository contains the data and code associated with image
analysis in the manuscript Simon et. al., (2025) XXXX Suppression of ERK
signalling promotes pluripotent epiblast in the human blastocyst. DOI
XXXX

# Repository structure and useage

## /ImageJ

ImageJ \*.ijm macro scripts for batch running StarDist segmentation on
confocal microscopy images of pre-implantation human embryos for later
quantification in Cell Profiler

For use with StarDist 2D plug-in (Schmidt et. al., 2018).

Modified pipeline from
<https://github.com/todd-fallesen/Niakan_Lab_KLF17> (Lea et. al., 2021)

## /Cell_Profiler

Pipelines for the segmentation and measurement of cytoplasmic and
nuclear fluorescence intensities of ImageJ processed confocal microscopy
images of pre-implantation human embryos.

-   [[FGF_ERKi_NANOG_GATA4_GATA3_nuclear_intensity.cppipe]{.underline}](https://gitlab.developers.cam.ac.uk/ctr/ctr-bioinformatics/niakan-lab/simon-et-al-2025/-/blob/main/Cell_Profiler/FGF_ERKi_NANOG_GATA4_GATA3_nuclear_intensity.cppipe?ref_type=heads "FGF_ERKi_NANOG_GATA4_GATA3_nuclear_intensity.cppipe") -
    pipeline used for nuclear intensity quantification for FGF and ERKi
    experiments

-   [[pERK_SOX2_OTX2_cyto_nuclear_intensity.cppipe]{.underline}](https://gitlab.developers.cam.ac.uk/ctr/ctr-bioinformatics/niakan-lab/simon-et-al-2025/-/blob/main/Cell_Profiler/pERK_SOX2_OTX2_cyto_nuclear_intensity.cppipe?ref_type=heads "pERK_SOX2_OTX2_cyto_nuclear_intensity.cppipe") -
    pipeline used for cytoplasmic and nuclear intensities for pERK
    experiments

These scripts were compiled in CellProfiler v4.1.3 (McQuin et. al.,
2018).

Modified pipeline from
<https://github.com/todd-fallesen/Niakan_Lab_KLF17> (Lea et. al., 2021)

## /data

Pre-processed and processed data for image analysis in R

/embryo_example - pre-processed data to be used with ../R/processing.R -
examples of Cell Profiler data outputs, and manual annotations for ICM
cells and corrections for under segmentation.

-   \*.csv - all processed data containing experiment metadata, raw and
    Z-corrected fluorescence intensitities for each of the FGF, ERKi and
    pERK experiments.

For raw imaging data, please see Figshare repository DOI dx.doi.org/XXXX

## /R

R scripts for data processing, analysis and figure generation

-   [[processing.R]{.underline}](https://gitlab.developers.cam.ac.uk/ctr/ctr-bioinformatics/niakan-lab/simon-et-al-2025/-/blob/main/R/processing.R?ref_type=heads "processing.R") -
    processing image quantification from Cell Profiler data output,
    combining with manual annotations and corrections, and performing
    Z-correction on fluorescence intensities

    -   [[ebcor.R]{.underline}](https://gitlab.developers.cam.ac.uk/ctr/ctr-bioinformatics/niakan-lab/simon-et-al-2025/-/blob/main/R/ebcor.R?ref_type=heads "ebcor.R") -
        function for Empirical Bayes correction of Z-axis associated
        fluorescence decay. Modified from
        <https://github.com/nestorsaiz/saiz-et-al_2016> (Saiz et. al.,
        2016)

-   \*\_analysis.R - analysis and figures of FGF or ERK treatments, or
    pERK experiments

    -   \*\_data.R - load processed experimental data

    -   \*\_hclust.R - hierarchical clustering of nuclear fluorescence
        intensities to assign ICM lineage identities. Modified from
        <https://github.com/nestorsaiz/saiz-et-al_2020> (Saiz et. al.
        2020)

These scripts were compiled in R Studio version 2023.12.1+402, and R
version 4.3.3
