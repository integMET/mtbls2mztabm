# mtbls2mztabm

Turn the ISA-TAB files of a [MetaboLights](https://www.ebi.ac.uk/metabolights/)
study into a Bioconductor `SummarizedExperiment` and an
[mzTab-M 2.1](https://github.com/HUPO-PSI/mzTab-M) file.

```r
library(mtbls2mztabm)

## study factors and their per-sample values
mtbls_factors("MTBLS103")
head(mtbls_factor_values("MTBLS103"))

## SummarizedExperiment of one assay (MAF abundances x samples)
se <- mtbls_summarized_experiment(
    "MTBLS103",
    assayName = "a_MTBLS103_ibanez_02_HILIC_mass_spectrometry.txt")

## the whole way in one call: ISA-TAB -> SummarizedExperiment -> mzTab-M
mtbls_to_mztabm("MTBLS103")          # writes MTBLS103.mztab

## a GC-TOF-MS study deposited as raw netCDF files works the same way
mtbls_to_mztabm("MTBLS40")           # writes MTBLS40.mztab
```

All ISA-TAB files are retrieved from the MetaboLights FTP server with the
[MsBackendMetaboLights](https://bioconductor.org/packages/MsBackendMetaboLights)
package. Every function also accepts a local directory of ISA-TAB files
through `path`; the files of MTBLS103 and MTBLS40 are shipped in
`inst/extdata`.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("MsBackendMetaboLights", "SummarizedExperiment",
                       "MsCoreUtils", "S4Vectors"))
## RmzTabM is not yet on Bioconductor
BiocManager::install("RforMassSpectrometry/RmzTabM")
## this package
BiocManager::install("mtbls2mztabm")
```

## What goes into the SummarizedExperiment

The object follows the layout of the `SummarizedExperiment` produced by the
xcms preprocessing workflow for MetaboLights data (MTBLS8735), so both can
be handled by the same code:

| Part | Content |
|---|---|
| `assays` | `raw` and `raw_filled`, the abundance matrix of the metabolite assignment file |
| `colData` | `sample_name`, `derived_spectra_data_file`, `metabolite_asssignment_file`, `source_name`, `organism`, `blood_sample_type`, `sample_type`, `age`, `unit`, `phenotype`, `injection_index`, `species`, `tissue`, `polarity`, `instrument` |
| `rowData` | `mzmed`, `mzmin`, `mzmax`, `rtmed`, `rtmin`, `rtmax`, `npeaks`, one count column per phenotype level, `ms_level`, `feature_id`, then the MAF annotations |

Fields a study does not deposit are `NA`.
