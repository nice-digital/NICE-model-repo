# Convert some of the Excel inputs into CSV files saved into "data/"

# Load functions for data extraction and cleaning
source("3_Functions/excel/extract.R")

# Define path to excel file
excel_path <- "1_Data/ID6184_RCC_model inputs FAD version [UK RWE unredacted, ACIC redacted, cPAS redacted].xlsm"

# Extract and clean data from excel
i <- f_excel_extract(excel_path, verbose = FALSE)
i <- c(i,f_excel_cleanParams(i$R_table_param))

# Extract desired parameters to input into a csv
# subset <- i[grep("List_pop|List_comparators", names(i))]
subset <- i[grep("List_pop", names(i))]

# Find the longest list, then pad the shorter lists with NA values
max_length <- max(sapply(subset, length))
padded <- lapply(subset, function(x) {
  length(x) <- max_length
  x
})

# Combine into a dataframe then save to csv
df <- t(as.data.frame(do.call(rbind, padded)))
write.csv(df, "shinyapp/data/lists.csv", row.names=FALSE)