
On this page, the hazard ratios (HR) from some pre-run network meta-analyses are applied to the extrapolated RWE survival curves generated on the previous page.


::: {.cell}

```{.r .cell-code}
# 3.3.5 Comparative efficacy propagation (NMA) ---------------------------------------------------------------

# Pull in the data and calculate means by pop line mol endpoint reftrt and reftrial

# First read in RDS file containing the PH NMA coda samples


# 3.3.5.1.1 PH NMA data -----------------------------------------------------

# Option to read in PH NMA CODA from local machine, uncomment this and comment out the line below to use
RDS_path2 <- file.path(d_path, "PH_NMA_CODA.rds")
if (file.exists(RDS_path2)) {
  i$PHNMA <- readRDS(RDS_path2)
} else {
  i$PHNMA <- readRDS(rstudioapi::selectFile(
    caption = "Please select 'PH_NMA_CODA.rds'",
    label = "PH_NMA_CODA.rds",
    path = "./1_Data/",
    filter = "R Files (*.rds)",
    existing = TRUE
  ))
}


colnames(i$PHNMA$data) <- c("Run", "Population", "Line", "Molecule", "Endpoint", "Reference.treatment", "Reference.trial", "HR")
i$PHNMA$data$Reference.endpoint <- i$PHNMA$data$Endpoint
```
:::


