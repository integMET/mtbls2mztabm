ext <- system.file("extdata", "MTBLS103", package = "mtbls2mztabm")
hilic <- "a_MTBLS103_ibanez_02_HILIC_mass_spectrometry.txt"
ref_coldata <- c("sample_name", "derived_spectra_data_file",
                 "metabolite_asssignment_file", "source_name", "organism",
                 "blood_sample_type", "sample_type", "age", "unit",
                 "phenotype", "injection_index", "species", "tissue",
                 "polarity", "instrument")
ref_rowdata <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax",
                 "npeaks")
