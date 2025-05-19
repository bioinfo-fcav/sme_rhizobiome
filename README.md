# sme_rhizobiome
Collection of R and bash scripts for exploring metabarcoding data related to the influence of Spirulina maxima extract on the rhizosphere microbiome of Solanum lycopersicum.

## Usage

The base directory (by default) is *~/sme_rhizobiome/*, which must contain the raw data files (\*_R1_\*.fastq and \*_R2_\*.fastq) stored in the *data/* subdirectory. The output of the analysis will be stored in the *analysis/* subdirectory.

Shell scripts are executed by passing command-line parameters, while R scripts are run within the R environment. The number in each script's filename indicates the order in which they should be executed.
