library(data.table)

PAS_price_rds_files <- dir(path = "./4_Output/PSA-PAS-price", pattern = "PSA_output.*\\.rds")
List_price_rds_files <- dir(path = "./4_Output/PSA-list-price", pattern = "PSA_output.*\\.rds")

PSA_results_PAS_price <- PAS_price_rds_files |>
  lapply(function(filename) file.path("./4_Output/PSA-PAS-price", filename)) |>
  lapply(readRDS) |>
  rbindlist()

setorder(PSA_results_PAS_price, iteration, oo_pop, trt_n)
PSA_results_PAS_price[, dd_drug_price_options := "PAS price"]

PSA_results_List_price <- List_price_rds_files |>
  lapply(function(filename) file.path("./4_Output/PSA-list-price", filename)) |>
  lapply(readRDS) |>
  rbindlist()

setorder(PSA_results_List_price, iteration, oo_pop, trt_n)
PSA_results_List_price[, dd_drug_price_options := "List price"]

combined_results <- rbind(PSA_results_PAS_price, PSA_results_List_price)

combined_results[, total_costs := mol_0 + mol_1 + mol_2 + mol_3 + mol_4 + mol_5 + mol_6 + mol_7 + mol_8 + mol_9 + mol_10 + mol_11 + mol_12 + mol_999 + other_costs]

saveRDS(combined_results, "./4_Output/PSA-combined.rds")
