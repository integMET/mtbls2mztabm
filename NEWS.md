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
