library(readr)
library(dplyr)


### Read in phenotypes ---------------------------------------------------------

main_phenos <- read_tsv("/WGHS/data/bsln/bsln_pheno_plink_28346.txt") %>%
  filter(RACE ==  1) %>%
  select(id=IID, age=AGERND, bmi=BMI, hdl=HDL, 
	 educ=EDUC, income=INCOME, 
	 alc=ALCUSE, currsmk=CURRSMK, pastsmk=PASTSMK) %>%
  mutate(hdl = ifelse(hdl == -9, NA, hdl),
         bmi = ifelse(bmi == -9, NA, bmi),
         log_hdl = log(hdl),
         age_sq = age^2,
	 smk = case_when(currsmk == 1 ~ "current", 
			 pastsmk == 1 ~ "former", 
			 (currsmk == 0) & (pastsmk == 0) ~ "never",
			 TRUE ~ "other"),
	 across(c(educ, income, alc), ~ as.factor(ifelse(.x == -9, "other", .x))))

gPC_phenos <- read_tsv("/WGHS/data/genotype/final/eigenvectors/eigenvectors_for_plink.txt") %>%
  select(id=IID, all_of(paste0("E", 1:10))) %>%
  rename_with(~ gsub("^E", "gPC", .x))

inflammation_phenos <- read_tsv("*****") %>%
  select(id, hscrp, fibrinogen, icam1)

migraine_phenos <- read_table("/WGHS/proj/migraine_latent_classes_R21/data/whs_migr08_v3_w_features_recoded_b_migid.txt") %>% 
  select(id=IID, migraine=actmig.2, headache=nomighead06.1) %>%
  mutate(headache = ifelse(is.na(headache), 0, headache))

diet_phenos <- read_csv("/WGHS/proj/KennyW_mediation/ahei_dataset.csv") %>%
  mutate(ahei = ifelse(is.na(ahei), median(ahei, na.rm=TRUE), ahei)) %>%
  select(id, fish, oily_fish, lcn3fa)

### Read in and merge genotypes ------------------------------------------------

rs295849_genos <- read_tsv("/WGHS/proj/snps_Kenny/extract_tabix/chr_17_36804493_dose.txt") %>%
  select(id, rs295849=dose)

genos <- full_join(rs2862183_genos, rs295849_genos, by="id")

### Read in and munge metabolite data (NMR and Metabolon) ----------------------

# NMR

nmr_df <- read_csv("/WGHS/data/bsln/WGHS_NMR_LP4_20161228.csv") %>%
  select(id=faux_id, everything()) %>%
  select(-comment_pl4_17)

glyca_df <- select(nmr_df, id, glyca = GLYCA_17)  # For use as an inflammatory biomarker

# Metabolon

load("/WGHS/proj/migraine_biomarker_signature/metabolon/qc/BRIG-01-22PHML+_DATA_TABLES_CORRECTED_batch_normalized.RData")
missingness_vec <- sapply(names(metabolon.data), function(m) sum(!is.na(metabolon.data[[m]])) / nrow(metabolon.data))
exclude_vec <- names(missingness_vec)[missingness_vec < 0.25]

metabolon_linker <- read_tsv("/WGHS/proj/migraine_biomarker_signature/metabolon/data/BRIG-01-22PHMLplus_sample_meta_data_WGHS.txt") %>% 
  select(PARENT_SAMPLE_NAME, id=CLIENT_SAMPLE_ID)

load("/WGHS/proj/migraine_biomarker_signature/metabolon/qc/BRIG-01-22PHML+_DATA_TABLES_CORRECTED_batch_normalized_imputed.RData")
n3fa_metab_ids <- c(epa_lcms = "X2050", dha_lcms = "X100000665", dpa_lcms = "X100001181")
metabolon_df <- metabolon.data %>%
  select(-all_of(exclude_vec)) %>%
  select(all_of(n3fa_metab_ids)) %>%
  mutate(n3fa_lcms = epa_lcms + dha_lcms + dpa_lcms) %>%
  rownames_to_column(var = "PARENT_SAMPLE_ID") %>%
  inner_join(metabolon_linker, by = "PARENT_SAMPLE_NAME") %>%
  select(-PARENT_SAMPLE_NAME)

load("/WGHS/proj/migraine_biomarker_signature/metabolon/qc/BRIG-01-22PHMLplus_chemical_annotation.RData")
metabolon_annot <- metabolites

# Merge

metabolite_df <- full_join(metabolon_df, nmr_df, by = "id")
  
# metabolon_analysis_df$income[metabolon_analysis_df$income == 1] <- NA  # Only a single person in the Metabolon subset

### Merge non-metabolite data --------------------------------------------------

pheno_geno_df <- main_phenos %>%
  left_join(gPC_phenos, by="id") %>%
  left_join(inflammation_phenos, by="id") %>%
  left_join(glyca_df, by = "id") %>%
  left_join(migraine_phenos, by="id") %>%
  left_join(diet_phenos, by="id") %>%
  left_join(genos, by="id")

### Save all objects to be used for analysis -----------------------------------

save("pheno_geno_df", 
     "metabolite_df", "metabolon_annot",
     file="analysis_objects.RData")
