library(data.table)
library(ukbnmr)


# Read in raw data
#nmr_file <- "/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb48298.tab.gz"
nmr_file <- "/humgen/florezlab/UKBB_app27892/UKB_RAP/nmr.10-11-2023.AKM.txt"
nmr_raw <- fread(nmr_file, data.table=FALSE, stringsAsFactors=FALSE)

# Preprocessing (log-transformation, adjustments, etc.; see ukbnmr package documentation)
processed <- remove_technical_variation(nmr_raw)

# Write preprocessed metabolite values and various metadata
fwrite(processed$biomarkers, file="../data/processed/nmr/nmr_data.csv")
fwrite(processed$biomarker_qc_flags, file="../data/processed/nmr/nmr_biomarker_qc_flags.csv")
fwrite(processed$sample_processing, file="../data/processed/nmr/nmr_sample_qc_flags.csv")
fwrite(processed$log_offset, file="../data/processed/nmr/nmr_biomarker_log_offset.csv")
fwrite(processed$outlier_plate_detection, file="../data/processed/nmr/nmr_outlier_plate_info.csv")
