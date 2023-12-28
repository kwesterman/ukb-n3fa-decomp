library(dplyr)


### Load analysis data ---------------------------------------------------------

load("analysis_objects.RData")
load("imputed_metabolite_data.RData")

analysis_df <- left_join(pheno_geno_df, imputed_metabolite_df, by = "id")

all_rsids <- grep("^rs", names(analysis_df), value = TRUE)

covar_sets <- list(
  basic = c("age", "age_sq"),
  primary = c("age", "age_sq", "educ", "income", "alc", "smk", "ahei"),
  with_gPCs = c("age", "age_sq", "educ", "income", "alc", "smk", "ahei", paste0("gPC", 1:5))
)

### Variable summaries ---------------------------------------------------------

pheno_geno_summary <- summary(analysis_df)

### Reproduce primary gene-PA interactions -------------------------------------

test_gxe <- function(g, e, y, df) {
  form_str <- paste0(y, " ~ ", g, " * ", e, 
                     paste(covar_sets$primary, collapse = " + "))
  lm_fit <- lm(as.formula(form_str), data = df)
  
}

replication_res_df <- expand_grid(
  g = all_rsids,
  e = c("fish", "etc", "n3fa_lcms", "n3fa_pred"),
  y = c("n3fa_lcms", "n3fa_pred", "hscrp", "etc")
) %>%
  rowwise() %>%
  mutate(model_fit = list(test_gxe(g, e, y, analysis_df)),
         model_res = filter(broom::tidy(model_fit), term == paste0(g, ":", e))) %>%
  unnest(model_res)

### Save results ---------------------------------------------------------------

save("pheno_geno_summary", "replication_res_df",
     file = "replication_results.RData")
