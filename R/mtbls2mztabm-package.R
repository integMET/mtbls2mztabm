#' @keywords internal
"_PACKAGE"

#' @importFrom MsBackendMetaboLights mtbls_ftp_path mtbls_list_files
#'     mtbls_metadata mtbls_sample_data
#' @importFrom MsCoreUtils retry
#' @importFrom RmzTabM MzTabM mtdSkeleton mtdFromSampleData mtdSample mtdSort
#'     setMtdField setMtdInstrument setMtdDatabase setMtdContact smlCreate
#'     smlAddStudyVariableColumns writeMzTabM
#' @importFrom RcppTOML parseTOML
#' @importFrom nanoparquet read_parquet write_parquet
#' @importFrom utils packageVersion read.table URLencode
#' @importFrom stats setNames
#' @importFrom tools file_ext file_path_sans_ext
NULL
