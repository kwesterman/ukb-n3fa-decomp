library(tidyverse)
library(data.table)


### Longitudinal data retrieval ------------------------------------------------

wrangle_long_data <- function(
    df, 
    field_ids, 
    merge_array_func = function(x) mean(x, na.rm=TRUE)
) {
  df %>%
    select(id = f.eid, matches(paste0("f\\.", field_ids, "\\.", collapse = "|"))) %>%
    filter(if_any(-id, ~ !is.na(.))) %>%
    pivot_longer(-id, names_to = "field_colname", values_to = "value") %>%
    mutate(field_colname = gsub("^f\\.", "", field_colname)) %>%
    mutate(split_cols = str_split(field_colname, "\\.")) %>%
    mutate(field = sapply(split_cols, "[[", 1),
           instance = sapply(split_cols, "[[", 2),
           array_idx = sapply(split_cols, "[[", 3)) %>%
    select(-field_colname, -split_cols) %>%
    group_by(id, field, instance) %>%
    summarise(value = merge_array_func(value), .groups = "drop") %>%
    mutate(field = names(field_ids)[match(field, field_ids)]) %>%
    pivot_wider(names_from = "field", values_from = "value")
}

### Basic variables ------------------------------------------------------------

base_pheno_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb10528.tab.gz",
                       data.table=FALSE, stringsAsFactors=FALSE)

ac_fields <- c(ac = 54, ac_date = 53)
centers <- c(
  "11012" = "Barts", "11021" = "Birmingham", "11011" = "Bristol", "11008" = "Bury",
  "11003" = "Cardiff", "11024" = "Cheadle_revisit", "11020" = "Croydon",
  "11005" = "Edinburgh", "11004" = "Glasgow", "11018" = "Hounslow", "11010" = "Leeds",
  "11016" = "Liverpool", "11001" = "Manchester", "11017" = "Middlesborough",
  "11009" = "Newcastle", "11013" = "Nottingham", "11002" = "Oxford",
  "11007" = "Reading", "11014" = "Sheffield", "10003" = "Stockport_pilot",
  "11006" = "Stoke", "11022" = "Swansea", "11023" = "Wrexham",
  "11025" = "Cheadle_pilot", "11027" = "Newcastle_pilot"
)
ac_long_df <- base_pheno_df %>%
  select(f.eid, contains("f.54."), contains("f.53.")) %>%
  mutate(across(everything(), as.character)) %>%
  wrangle_long_data(ac_fields, merge_array_func = function(x) x[1]) %>%
  mutate(id = as.integer(id),
         ac = centers[ac])  # Recode numbers to location names

basic_fields <- c(sex = 31, age = 21003,  bmi = 21001, 
                  sbp = 4080, dbp = 4079, fasting_hrs = 74)
basic_phenos_long_df <- base_pheno_df %>%
  wrangle_long_data(basic_fields) %>%
  group_by(id) %>%
  mutate(sex = ifelse(is.na(sex), sex[instance == 0], sex)) %>%
  ungroup() %>%
  mutate(age_squared = age^2,
         ageBySex = age * sex)

### Covariates and confounders -------------------------------------------------

income_fields <- c(income = 738)
income_coding <- c(
  "1" = "Less than 18,000",
  "2" = "18,000 to 30,999",
  "3" = "31,000 to 51,999",
  "4" = "52,000 to 100,000",
  "5" = "Greater than 100,000",
  "-1" = "Do not know",
  "-3" = "Prefer not to answer"
)
income_long_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_may_2023/ukb672750.tab.gz",
                        data.table=FALSE, stringsAsFactors=FALSE) %>%
  wrangle_long_data(income_fields) %>%
  mutate(income = income_coding[as.character(income)])

education_fields <- c(education = 6138)
education_coding <- c(
  "1" = "College or University degree",
  "2" = "A levels/AS levels or equivalent",
  "3" = "O levels/GCSEs or equivalent",
  "4" = "CSEs or equivalent",
  "5" = "NVQ or HND or HNC or equivalent",
  "6" = "Other professional qualifications eg: nursing, teaching",
  "-7" = "None of the above",
  "-3" = "Prefer not to answer"
)
education_long_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_aug_2022/ukb669148.tab.gz",
                           data.table=FALSE, stringsAsFactors=FALSE) %>%
  wrangle_long_data(education_fields) %>%
  mutate(education = education_coding[as.character(education)])

lifestyle_fields <- c(smoking = 20116)
smoking_coding <- c(
  "-3" = "Prefer not to answer",
  "0" =	"Never",
  "1" =	"Previous",
  "2" =	"Current"
)
lifestyle_long_df <- base_pheno_df %>%
  wrangle_long_data(lifestyle_fields) %>%
  mutate(smoking = smoking_coding[as.character(smoking)])

alcohol_fields <- c(alcohol = 1558)
alcohol_coding <- c(
  "1" = "Daily or almost daily",
  "2" = "Three or four times a week",
  "3" = "Once or twice a week",
  "4" = "One to three times a month",
  "5" = "Special occasions only",
  "6" = "Never",
  "-3" = "Prefer not to answer"
)
alcohol_long_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb40167.tab.gz",
                         data.table=FALSE, stringsAsFactors=FALSE) %>%
  wrangle_long_data(alcohol_fields) %>%
  mutate(alcohol = alcohol_coding[as.character(alcohol)])

covariate_long_df <- income_long_df %>%
  full_join(education_long_df, by=c("id", "instance")) %>%
  full_join(lifestyle_long_df, by=c("id", "instance")) %>%
  full_join(alcohol_long_df, by=c("id", "instance"))

saveRDS(covariate_long_df, "covariate_long_df.rds")

### Biomarkers -----------------------------------------------------------------

bm_fields <- c(
  alt = 30620, alb = 30600, apoB = 30640, hscrp = 30710, chol = 30690, glu = 30740, 
  hba1c = 30750, hdl = 30760, ldl = 30780, shbg = 30830, tg = 30870, vitD = 30890
)

biomarker_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb28679.tab.gz", 
                      data.table=FALSE, stringsAsFactors=FALSE)
biomarker_long_df <- biomarker_df %>%
  wrangle_long_data(bm_fields)

# Collect medication data for biomarker adjustments
drug_df <- base_pheno_df %>%
  select(id=f.eid, contains("f.20003.0"), contains("f.6177.0"))

statin_ids <- c(
  1140861958, 1140888594, 1140888648, 1141146234, 1141192410, 1140861922, 1141146138
)
drug_df$num_statins <- 0
for (statin in statin_ids) {
  drug_df$num_statins <- drug_df$num_statins + rowSums(drug_df[, grep("20003", names(drug_df), value=TRUE)] == statin, na.rm=TRUE)
}
drug_df$bp_med <- rowSums(drug_df[, (
  grep("6177", names(drug_df), value=TRUE)
)] == 2, na.rm=TRUE)
drug_df <- select(drug_df, id, num_statins, bp_med)

# Make medication-based biomarker adjustments
biomarker_long_df <- left_join(biomarker_long_df, 
                               select(drug_df, id, num_statins),
                               by="id")  # Not currently joining by instance!
statin_adj_bms <- c("chol", "ldl", "apoB")
statin_adj_factors <- c(
  chol = 0.749,
  ldl = 0.684,
  apoB = 0.719
)
for (bm in statin_adj_bms) {
  adj_factor <- ifelse(biomarker_long_df$num_statins > 0, statin_adj_factors[bm], 1)
  biomarker_long_df[[paste0(bm, "_statinadj")]] <- biomarker_long_df[[bm]] / adj_factor
}

basic_phenos_long_df <- left_join(basic_phenos_long_df, 
                                  select(drug_df, id, bp_med), 
                                  by="id")
bp_adj_factors <- c(sbp = 15, dbp = 10)
for (bm in c("sbp", "dbp")) {
  adj_factor <- ifelse(basic_phenos_long_df$bp_med, bp_adj_factors[bm], 0)
  basic_phenos_long_df[[paste0(bm, "_medsadj")]] <- basic_phenos_long_df[[bm]] + adj_factor
}

saveRDS(biomarker_long_df, "biomarker_long_df.rds")
saveRDS(basic_phenos_long_df, "basic_phenos_long_df.rds")

### Outcomes for sample exclusion ----------------------------------------------

medical_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb38040.tab.gz",
                    data.table=FALSE, stringsAsFactors=FALSE) %>%
  select(id=f.eid, everything())
medical_df$diabetes <- rowSums(medical_df[, grepl("f\\.2443\\.0\\.", names(medical_df)), drop=FALSE] == 1, na.rm=T) > 0
medical_df$MI <- rowSums(medical_df[, grepl("f\\.6150\\.0\\.", names(medical_df)), drop=FALSE] == 1, na.rm=T) > 0
medical_df$angina <- rowSums(medical_df[, grepl("f\\.6150\\.0\\.", names(medical_df)), drop=FALSE] == 2, na.rm=T) > 0
cirrhosis_codes <- c(paste0("K70", 2:4), "K717", paste0("K74", 0:6))
cirrhosis_primary_ids <- c()
for (f in grep("f\\.41202\\.0\\.", names(medical_df), value=T)) {
  cirrhosis_primary_ids <- c(cirrhosis_primary_ids, medical_df$id[medical_df[[f]] %in% cirrhosis_codes])
}
cirrhosis_secondary_ids <- c()
for (f in grep("f\\.41204\\.0\\.", names(medical_df), value=T)) {
  cirrhosis_secondary_ids <- c(cirrhosis_secondary_ids, medical_df$id[medical_df[[f]] %in% cirrhosis_codes])
}
medical_df$pregnant <- rowSums(medical_df[, grepl("f\\.3140\\.0\\.", names(medical_df)), drop=FALSE] == 1, na.rm=T) > 0
cancer_tmp <- select(medical_df, 1, contains("f.40005."))
cancer <- cancer_tmp[, 1:7] %>%
  inner_join(select(filter(ac_long_df, instance == 0), id, ac_date), by="id")  # Add assessment center dates
cancer$ac_date = as.Date(cancer$ac_date)
for (i in 2:7) {
  x <- ifelse(abs(difftime(cancer[, i, drop=TRUE], cancer$ac_date, units="days")) <= 365, TRUE, FALSE)  # TRUE if cancer diagnosis within a year of assessment center visit
  cancer <- cbind(cancer, x)
}
cancer$cancer_within_1yearac = apply(cancer[, 9:14], 1, function(x) {
  ifelse(any(x == TRUE, na.rm=TRUE), TRUE, FALSE)
})
cancer[names(cancer) == "x"] <- NULL

medical_df <- medical_df %>%
  left_join(select(cancer, id, cancer_within_1yearac), by="id") %>%
  mutate(CHD = MI | angina,
         cirrhosis = id %in% c(cirrhosis_primary_ids, cirrhosis_secondary_ids)) %>%
  select(id, diabetes, CHD, cirrhosis, pregnant, cancer_within_1yearac)

saveRDS(medical_df, "medical_df.rds")

### Dietary data ---------------------------------------------------------------

# 24-hour recall

calc_diet_fields <- function(field_ids, df, coding=FALSE) {
  # Given a list of fields constituting a food group:
  # - Determine the set of 24HR that are valid for that food group
  # - Recode the relevant variables based on their codings if necessary
  # - Sum over all fields for that food group
  valid_24hr <- (findInterval(df$TCALS_qc / 4.18, c(600, 4800)) == 1) &
    df$typical_diet_qc == 1
  df <- select(df, all_of(field_ids))
  if (coding) {  # Recode the variable if necessary (for food groups)
    df <- mutate_all(df, ~codings[as.character(.)])
  }
  ifelse(valid_24hr,  # Sum over fields if valid 24HR, else NA
         rowSums(df, na.rm=TRUE),
         NA)
}

calc_med_score <- function(diet_df) {
  # For additional details, see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6842574/
  med_score_df <- diet_df %>%
    mutate(veg_score = VEG > median(VEG, na.rm=TRUE),
           legume_score = LEGUMES > median(LEGUMES, na.rm=TRUE),
           fruit_score = FRUIT > median(FRUIT, na.rm=TRUE),
           nut_score = NUTS > median(NUTS, na.rm=TRUE),
           fish_score = FISH > median(FISH, na.rm=TRUE),
           whgrain_score = WHGRAIN > median(WHGRAIN, na.rm=TRUE),
           mufa2sfa_score = MUFA2SFA > median(MUFA2SFA, na.rm=TRUE),
           redprmeat_score = REDPRMEAT < median(REDPRMEAT, na.rm=TRUE),
           alc_score = (ALC > 5) & (ALC < 25)) %>%
    mutate_at(vars(veg_score, legume_score, fruit_score, nut_score, fish_score, whgrain_score, mufa2sfa_score, redprmeat_score, alc_score),
              ~ifelse(is.na(.), mean(., na.rm=TRUE), .)) %>%
    mutate(mds = veg_score + legume_score + fruit_score + nut_score + fish_score + whgrain_score + mufa2sfa_score + redprmeat_score + alc_score)
  med_score_df$mds
}

winsorize <- function(x, SDs=3) {
  lims <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x[x < lims[1]] <- lims[1]
  x[x > lims[2]] <- lims[2]
  x
}

diet_qc_fields_list <- list(typical_diet = "100020", TCALS = "26002")
diet_nutrient_fields_list <- list(
  CHO = "26013", SUCROSE = "26059", 
  FAT = "26008", SFA = "26014", MUFA = "26032",
  N3FA = "26015", N6FA = "26016",
  PRO = "26005", ALC = "26030", FIBER = "26017", 
  NUT_FE = "26019", NUT_NA = "26025", NUT_K = "26024", NUT_FOL = "26022",
  NUT_MG = "26025", NUT_VITC = "26023", NUT_VITD = "26029"
)
# diet_nutrient_fields_list <- list(
#   CHO = "100005", SUGARS = "100008", 
#   FAT = "100004", SFA = "100006", PUFA = "100007",
#   PRO = "100003", ALC = "100022", FIBER = "100009", 
#   NUT_FE = "100011", NUT_K = "100016", NUT_FOL = "100014",
#   NUT_MG = "100017", NUT_VITC = "100015", NUT_VITD = "100021"
# )
diet_food_group_fields_list <- list(
  VEG =  as.character(seq(104060, 104380, 10)),  
  LEGUMES = as.character(c(104000, 104010)),  # Broad beans (104110) and green beans (104120) here instead of veg?
  FRUIT = as.character(seq(104410, 104590, 10)),
  NUTS = as.character(seq(102410, 102440, 10)),
  FISH = as.character(seq(103150, 103230, 10)),
  OILY_FISH = "103160",
  WHGRAIN = as.character(c(101260, 102720, 102740, 100800, 100840, 100850)),
  REDPRMEAT = as.character(c(103010, 103020, 103030, 103040, 103050, 103070, 103080))
)
diet_fields_list <- c(diet_qc_fields_list, diet_nutrient_fields_list, 
                      diet_food_group_fields_list)
all_diet_vars <- c(names(diet_fields_list), "MUFA")
codings <- c(
  "1"=1, "2"=2, "3"=3, "4"=4, "5"=5, 
  "100"=1, "200"=2, "300"=3, "400"=4, "500"=5, "600"=6,
  "444"=0.25, "555"=0.5
)

diet_24hr_primary_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb22861.tab.gz", 
                      data.table=FALSE, stringsAsFactors=FALSE)
diet_24hr_nutrients_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_aug_2022/ukb670995.tab.gz", 
                                data.table=FALSE, stringsAsFactors=FALSE)
diet_24hr_df <- full_join(diet_24hr_primary_df, diet_24hr_nutrients_df, 
                          by=c("f.eid"))
rm(diet_24hr_primary_df, diet_24hr_nutrients_df)

# Do the 24HR wrangling group-by-group to keep memory usage under control
diet_24hr_qc_fields <- setNames(c("100020", "26002"), 
                                c("typical_diet_qc", "TCALS_qc"))
diet_24hr_qc_df <- diet_24hr_df %>%  # Get fields for "typical diet" and total calories to be used for QC
  wrangle_long_data(diet_24hr_qc_fields)

diet_24hr_long_list <- lapply(names(diet_fields_list), function(grp) {
  print(grp)
  grp_fields <- diet_fields_list[[grp]]
  is_food_group <- (grp %in% names(diet_food_group_fields_list))
  diet_24hr_df %>%
    filter(if_any(contains("f.26002."), ~ !is.na(.))) %>%
    wrangle_long_data(setNames(grp_fields, grp_fields)) %>%
    inner_join(diet_24hr_qc_df, by=c("id", "instance")) %>%
    mutate(!!grp := calc_diet_fields(grp_fields, ., coding = is_food_group)) %>%
    select(id, instance, !!grp)
})
diet_24hr_long_df <- reduce(diet_24hr_long_list, function(x, y) {
  full_join(x, y, by=c("id", "instance"))
}) %>%
  mutate(PUFA = FAT - SFA - MUFA,
         TCALS = TCALS / 4.18) %>%  # Energy from kJ to kcals
  mutate_at(vars(CHO, PRO), ~. * 4) %>%  # Nutrients from g to kcals (other than alcohol)
  mutate_at(vars(FAT, SFA, MUFA, PUFA), ~. * 9) %>%
  mutate_at(vars(all_of(all_diet_vars)), winsorize) %>%
  select(id, instance, all_of(all_diet_vars)) %>%
  mutate(MUFA2SFA = MUFA / SFA) %>%
  mutate(mds = calc_med_score(.))

supp_24hr_long_df <- diet_24hr_df %>%  # Supplements reported on 24HR
  wrangle_long_data(c(fish_oil_24hr = 20084),
                    function(x) ifelse(all(is.na(x)), NA, 
                                       as.integer(any(x == 472, na.rm = TRUE))))

diet_supp_24hr_long_df <- diet_24hr_long_df %>%
  left_join(supp_24hr_long_df, supp_24hr_long_df, by=c("id", "instance"))

saveRDS(diet_supp_24hr_long_df, "diet_supp_24hr_long_df.rds")

write_csv(diet_supp_24hr_long_df, "../data/processed/all_24hr_data.csv")

diet_supp_24hr_collapsed_df <- diet_supp_24hr_long_df %>%  # Collapse b/c these responses don't correspond to exam visits
  filter(!is.na(TCALS)) %>%
  group_by(id) %>%
  summarise(across(everything(), ~ mean(., na.rm = TRUE)),
            num_recalls = n()) %>%
  mutate(instance = 0)

# Food frequency questionnaire

ffq_cat_to_qt <- function(x) {
  case_when(  # Data-coding 100377
    x == 5 ~ 1,  # "Once or more daily"
    x == 4 ~ 5.5 / 7,  # "5-6 times a week"
    x == 3 ~ 3 / 7,  # "2-4 times a week"
    x == 2 ~ 1 / 7,  # "Once a week"
    x == 1 ~ 0.5 / 7,  # "Less than once a week"
    x == 0 ~ 0,  # "Never"
    TRUE ~ as.numeric(NA)
  )
}

ffq_fields <- c(oily_fish = 1329, nonoily_fish = 1339,
                prmeat = 1349, poultry = 1359,
                beef = 1369, lamb = 1379)
ffq_long_df <- base_pheno_df %>%
  wrangle_long_data(ffq_fields) %>%
  mutate(across(-c(id, instance), ffq_cat_to_qt))

saveRDS(ffq_long_df, "ffq_long_df.rds")

# Supplements (from touchscreen or verbal interview)

extra_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_may_2023/ukb672750.tab.gz",
                  data.table=FALSE, stringsAsFactors=FALSE)
supp_touchscreen_long_df <- extra_df %>%  # Supplements reported on touchscreen
  wrangle_long_data(c(fish_oil = 6179),
                    function(x) ifelse(all(is.na(x)), NA, 
                                       as.integer(any(x == 1, na.rm=TRUE))))
supp_verbal_long_df <- extra_df %>%  # Supplements reported in verbal interview
  wrangle_long_data(c(fish_oil_verbal = 20003),
                    function(x) ifelse(all(is.na(x)), NA, 
                                       as.integer(any(x == 1193, na.rm=TRUE))))

supp_long_df <- supp_touchscreen_long_df %>%
  inner_join(supp_verbal_long_df, by=c("id", "instance"))

saveRDS(supp_long_df, "supp_long_df.rds")

# Merge dietary data from various sources

full_diet_long_df <- diet_supp_24hr_collapsed_df %>%
  full_join(ffq_long_df, by=c("id", "instance")) %>%
  full_join(supp_long_df, by=c("id", "instance"))

### Add genetic PCs and relatedness --------------------------------------------

gPC_df <- base_pheno_df %>%
  select(id=f.eid, contains("f.22009.0"), used_in_PCA=f.22020.0.0) %>%
  rename_with(.fn=~gsub("f.22009.0.", "gPC", .)) %>%
  mutate(unrelated = (used_in_PCA == 1)) %>%  # An unrelated subset was used in the central PCA
  select(id, all_of(paste0("gPC", 1:20)), unrelated)

### Merge and write "raw" phenotypes -------------------------------------------

# ac_long_df <- readRDS("ac_long_df.rds")
# biomarker_long_df <- readRDS("biomarker_long_df.rds")
# basic_phenos_long_df <- readRDS("basic_phenos_long_df.rds")
# diet_24hr_long_df <- readRDS("diet_24hr_long_df.rds")
# ffq_long_df <- readRDS("ffq_long_df.rds")
# supp_extra_long_df <- readRDS("supp_extra_long_df.rds")

withdrawn_consent <- scan("/humgen/florezlab/UKBB_app27892/withdraw27892_232_14_Nov_2022.txt", what=character())

phenos <- ac_long_df %>%
  left_join(basic_phenos_long_df, by=c("id", "instance")) %>%
  left_join(covariate_long_df, by=c("id", "instance")) %>%
  left_join(biomarker_long_df, by=c("id", "instance")) %>%
  left_join(medical_df, by="id") %>%
  left_join(full_diet_long_df, by=c("id", "instance")) %>%
  left_join(gPC_df, by="id") %>%
  filter(!(id %in% withdrawn_consent)) %>%
  mutate(id = format(id, scientific=FALSE)) %>%
  mutate(across(contains("gPC"), ~. * mds, .names="mdsBy{.col}"))

write_csv(phenos, "../data/processed/ukb_phenos_longitudinal_raw.csv")

### Phenotype processing and exclusions ----------------------------------------

logged_risk_factors <- c("alt", "tg", "hscrp")
risk_factors <- c(
  "bmi",
  "sbp", "dbp", "sbp_medsadj", "dbp_medsadj",
  "alt_log", 
  "chol", "ldl", "hdl", "apoB",
  "tg_log",
  "hba1c", "glu",
  "hscrp_log",
  "vitD",
  "chol_statinadj", "ldl_statinadj", "apoB_statinadj"
)

processed_phenos <- phenos %>%
  filter(!diabetes & !CHD & !cirrhosis & !cancer_within_1yearac & !pregnant) %>%
  mutate(across(one_of(logged_risk_factors), list(log=log))) %>%
  mutate(across(one_of(risk_factors), 
                ~ifelse(findInterval(., mean(., na.rm=TRUE) + c(-5, 5) * sd(., na.rm=TRUE)) != 1, 
                        as.numeric(NA), .)))

### Write processed phenotypes -------------------------------------------------

write_csv(processed_phenos, "../data/processed/ukb_phenos_longitudinal.csv")

processed_phenos %>%
  filter(unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_longitudinal_unrelated.csv")

### Add Pan-UKBB data to generate European subset ------------------------------

anc_rel_df <- fread("/humgen/florezlab/UKBB_app27892/ukbreturn2442/all_pops_non_eur_pruned_within_pop_pc_covs_app27892.csv",
                    data.table=FALSE, stringsAsFactors=FALSE) %>%
  mutate(f.eid = as.character(f.eid),
         unrelated = !related_return2442) %>%
  select(id=f.eid, ancestry=pop_return2442, unrelated,
         one_of(paste0("PC", 1:10, "_return2442"))) %>%
  rename_at(vars(contains("PC")), ~gsub("_return2442", "", .)) %>%
  rename_at(vars(contains("PC")), ~gsub("PC", "gPC", .))

processed_phenos_panUKBB <- processed_phenos %>%
  select(-contains("gPC"), -unrelated) %>%
  inner_join(anc_rel_df, by="id") %>%
  mutate(across(contains("gPC"), ~. * mds, .names="mdsBy{.col}"))

processed_phenos_panUKBB %>%
  filter(ancestry == "EUR") %>%
  write_csv("../data/processed/ukb_phenos_longitudinal_EUR.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "EUR", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_longitudinal_EUR_unrelated.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "AFR", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_longitudinal_AFR_unrelated.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "AMR", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_longitudinal_AMR_unrelated.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "CSA", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_longitudinal_CSA_unrelated.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "EAS", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_longitudinal_EAS_unrelated.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "MID", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_longitudinal_MID_unrelated.csv")
