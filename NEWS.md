# mtbls2mztabm 0.99.1

## SIGNIFICANT USER-VISIBLE CHANGES

- `mtbls_to_mztabm()` writes one mzTab-M file for each metabolite assignment
  file (MAF) of a study instead of only one for its first mass spectrometry
  assay. With several MAFs, the file names get the part of the MAF name that
  tells them apart (`MTBLS1066_NEG_RP_LCMS.mztab`) and the value holds the paths
  of all files. A MAF that fails gives a warning and the others are still
  written. A MAF shared by several assays gets one file from the first of them;
  a MAF no assay refers to and the MAFs of assays of another technology get
  none. `assayName` still writes the single file of one assay.
- `mtbls_to_mztabm()` no longer goes through a `SummarizedExperiment`. It runs
  the workflow study identifier -> study description (TOML) and metabolite
  assignment file (Parquet) -> `MTD` section (from the TOML file) and `SML`
  section (from the Parquet file) -> mzTab-M file. The two intermediate files
  are kept next to the mzTab-M file (arguments `toml` and `parquet`).
- The mzTab-M files hold a small molecule summary (`SML`, profile `M+S`)
  instead of a small molecule feature table (`SMF`, profile `M+F`). A MAF
  reports metabolites and mostly no m/z, which `SMF` rows require, so the old
  files did not pass the mzTab-M validation of *RmzTabM* (>= 0.99.3). The
  files of the bundled studies MTBLS103 and MTBLS40 now validate without a
  message.
- Removed `mtbls_summarized_experiment()` and `se_to_mztabm()`, and the
  dependencies on *SummarizedExperiment* and *S4Vectors*.
- The attributes of the value of `mtbls_to_mztabm()` are `"toml"`, `"parquet"`
  and `"mztabm"` (`"se"` is gone). The study variable group is named after the
  study factor (e.g. `HIAE`) rather than `phenotype`.

## NEW FEATURES

- `mtbls_factors_to_toml()` saves the study factors (sample groups), the
  samples, the runs and the study metadata of an assay as a TOML file.
- `mtbls_maf_to_parquet()` saves the metabolite assignment file of an assay as a
  Parquet file with stable column types.
- `mtbls_toml_to_mtd()` builds the mzTab-M `MTD` section from the TOML file
  with *RmzTabM*. Several study factors can be used as study variables
  (`factor`), samples without ontology terms are supported and runs without a
  data file get the location of the study.
- `mtbls_parquet_to_sml()` builds the `SML` section from the Parquet file with
  *RmzTabM*, including the average and variation of every study variable.
  Retention times and m/z of the MAF and annotations that the `SML` columns
  cannot hold (for example the names of unidentified peaks) are kept in
  optional columns.
- `mtd_sml_to_mztabm()` combines both sections into an `RmzTabM::MzTabM`
  object, adds the compound databases found in the identifiers and, with
  `file`, writes the mzTab-M file.
- `mtbls_to_mztabm()` gained the arguments `polarity` (for assay tables
  without scan polarity), `reliability`, `toml`, `parquet` and `validate`;
  `verbose = FALSE` now silences all messages.

## BUG FIXES

- A scan polarity that is missing for only some of the runs no longer raises
  an "Unsupported scan polarity" error.
- The MS runs of the `MTD` section refer to the instrument
  (`ms_run[n]-instrument_ref`). The current mzTab-M validator rejects an
  instrument without such a reference, so the files of all studies with an
  instrument failed the validation.
- The abundance of a study variable with the average function `mean` is
  calculated from the available abundances of its assays
  (`mean(x, na.rm = TRUE)`). It was `null` if any of them was missing and is
  now `null` only if all of them are.

# mtbls2mztabm 0.99.0

## NEW FEATURES

- `mtbls_to_mztabm()` writes an mzTab-M 2.1 file for a MetaboLights study
  from its ISA-TAB files, taking only the study identifier.
- `mtbls_summarized_experiment()` builds a `SummarizedExperiment` of one
  assay in the layout of the xcms MetaboLights workflow (assays `raw` and
  `raw_filled`, the same 15 `colData()` and 12 leading `rowData()` columns).
- `se_to_mztabm()` maps such an object onto an `RmzTabM::MzTabM` object,
  filling instrument, database, title, description and publication from the
  ISA-TAB files.
- `mtbls_investigation()`, `mtbls_factors()`, `mtbls_factor_values()` and
  `mtbls_maf()` give access to the investigation file, the study factors
  and the metabolite assignment file.
- All functions read from the MetaboLights FTP server through
  MsBackendMetaboLights by default and accept a local ISA-TAB directory
  through `path`; the ISA-TAB files of MTBLS103 (LC-MS) and MTBLS40
  (GC-TOF-MS) ship with the package.
- Where an assay table lists no derived spectral data files, the raw files
  (e.g. the netCDF files of GC-MS studies such as MTBLS40) are used as the
  mzTab-M `ms_run`s.
