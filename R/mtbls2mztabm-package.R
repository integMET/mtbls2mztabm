#' @keywords internal
"_PACKAGE"

#' @importFrom methods is
#' @importFrom SummarizedExperiment SummarizedExperiment colData rowData
#'     assayNames assay colData<- rowData<-
#' @importFrom S4Vectors DataFrame
#' @importFrom MsBackendMetaboLights mtbls_ftp_path mtbls_list_files
#'     mtbls_metadata mtbls_sample_data
#' @importFrom MsCoreUtils retry
#' @importFrom RmzTabM MzTabM sampleCols msRunCols assayCols smfCols
#'     setMtdField setMtdInstrument setMtdDatabase setMtdContact writeMzTabM
#' @importFrom utils packageVersion read.table
#' @importFrom stats setNames
#' @importFrom tools file_ext file_path_sans_ext
NULL
