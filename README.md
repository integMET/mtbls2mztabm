# mtbls2mztabm

Turn the ISA-TAB files of a [MetaboLights](https://www.ebi.ac.uk/metabolights/)
study into an [mzTab-M 2.1](https://github.com/HUPO-PSI/mzTab-M) file.

The conversion runs through two intermediate files that can be inspected,
edited and reused: the study factors (sample groups) and study metadata go to
a [TOML](https://toml.io) file, the metabolite assignment file (MAF) to an
[Apache Parquet](https://parquet.apache.org) file. The *metadata* (`MTD`)
section of the mzTab-M file is built from the first and the *small molecule*
(`SML`) section from the second, both with the
[RmzTabM](https://github.com/RforMassSpectrometry/RmzTabM) package.

```
                   ┌─ factors, samples, runs, metadata ─▶ <study>.toml    ─▶ MTD ─┐
MetaboLights ID ───┤                                                              ├─▶ <study>.mztab
                   └─ metabolite assignment file (MAF) ─▶ <study>.parquet ─▶ SML ─┘
```

```r
library(mtbls2mztabm)

## the whole way in one call: writes MTBLS103.mztab and, next to it, the
## intermediate files MTBLS103.toml and MTBLS103.parquet
mtbls_to_mztabm("MTBLS103")

## a GC-TOF-MS study deposited as raw netCDF files works the same way
mtbls_to_mztabm("MTBLS40")
```

All ISA-TAB files are retrieved from the MetaboLights FTP server with the
[MsBackendMetaboLights](https://bioconductor.org/packages/MsBackendMetaboLights)
package. Every function also accepts a local directory of ISA-TAB files
through `path`; the files of MTBLS103 and MTBLS40 are shipped in
`inst/extdata`.

## Step by step

Each step of `mtbls_to_mztabm()` is a function of its own:

```r
## 1. study factors, samples, runs and study metadata of an assay -> TOML
mtbls_factors_to_toml("MTBLS103", "MTBLS103.toml")

## 2. metabolite assignment file of the assay -> Parquet
mtbls_maf_to_parquet("MTBLS103", "MTBLS103.parquet")

## 3. MTD section from the TOML file, 4. SML section from the Parquet file
mtd <- mtbls_toml_to_mtd("MTBLS103.toml")
sml <- mtbls_parquet_to_sml("MTBLS103.parquet", "MTBLS103.toml", mtd)

## 5. combine both and write the mzTab-M file
m <- mtd_sml_to_mztabm(mtd, sml, file = "MTBLS103.mztab")
```

The TOML file is meant to be read and edited, for example to add a missing
scan polarity or to correct an ontology term, before the metadata section is
built from it. With several assays in a study, `assayName` selects the assay;
`factor` selects the study factors that become study variables.

## What goes into the files

| File | Content |
|---|---|
| TOML | `[study]` identifier, title, description, URI; `[source]` assay and MAF file; `[instrument]`; `[[publication]]`; `[[contact]]` (with ORCID); `[[factor]]` the study factors with their groups; `[[sample]]` organism, organism part and the group of the sample in every factor; `[[run]]` per assay-table row: assay name, sample, data file, scan polarity and the MAF column with its abundances |
| Parquet | the MAF as it is, with stable column types (text annotations, numeric abundances) and the study identifier and MAF name as file metadata |
| mzTab-M `MTD` | samples, MS runs, assays, study variable groups and study variables, instrument, publications, URI, software, compound databases |
| mzTab-M `SML` | one row per MAF row: identifiers, formula, structures, name, reliability, the abundance of every assay and the average and variation of every study variable |

Fields a study does not deposit are left out. mzTab-M needs a few values the
ISA-TAB files rarely provide: an identification level (`reliability`) for every
small molecule, which is reported as the lowest level (`4`) unless the MAF
gives one (argument `reliability`), and an ORCID for a contact, so contacts are
only written if the investigation file has one.

## Validation

`mtbls_to_mztabm()` checks the file it writes with the official mzTab-M
validator (`validate = TRUE`, the default), which `RmzTabM` downloads on first
use. Use `validate = FALSE` to skip the check, for example without internet
access.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("MsBackendMetaboLights", "MsCoreUtils",
                       "nanoparquet", "RcppTOML"))
## RmzTabM is not yet on Bioconductor
BiocManager::install("RforMassSpectrometry/RmzTabM")
## this package
BiocManager::install("mtbls2mztabm")
```

## Changes from 0.99.0

The `SummarizedExperiment` interface (`mtbls_summarized_experiment()` and
`se_to_mztabm()`) was removed together with the dependencies on
*SummarizedExperiment* and *S4Vectors*. The mzTab-M files now hold a
small molecule summary (`SML`) instead of a small molecule feature table
(`SMF`): a MAF reports metabolites rather than features, and the mandatory
`exp_mass_to_charge` of an `SMF` row, which MAFs mostly lack, made the files
fail the mzTab-M validation. See `NEWS.md`.
