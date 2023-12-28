# library(readr)
# library(dplyr)
library(caret)
library(glmnet)


### Load metabolite data -------------------------------------------------------

load("analysis_objects.RData")

### Train predictive model for total N3FA --------------------------------------

# Set up training / prediction split

full_df <- select(metabolite_df, -id, -c(epa, dha, dpa))  # Leaving only the n3fa_lcms outcome of interest and NMR variables

training_df <- filter(full_df, !is.na(n3fa_lcms))

# prediction_df <- filter(full_df, is.na(n3fa_lcms))

# Train model using cross-validation

set.seed(123)
cv5 <- trainControl(method = "cv", number = 5)

enet_model <- train(
  n3fa_lcms ~ ., 
  data = training_df,
  method = "glmnet", 
  trControl = cv5, 
  tuneLength = 5  # Re-visit glmnet tuning with caret
)

imputed_metabolite_df <- bind_cols(id = metabolite_df$id, full_df) %>%
  mutate(n3fa_pred = predict(enet_model, newdata = full_df)) %>%
  select(id, n3fa_lcms, n3fa_pred)

### Save all objects to be used for analysis -----------------------------------

save("imputed_metabolite_df", file="imputed_metabolite_data.RData")
