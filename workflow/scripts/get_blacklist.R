# First try to download Cut & Run specific blacklists from
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03027-3#MOESM2
# If not available for given genome, download the ChIP-seq blacklist

# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(openxlsx)

genome <- snakemake@params["genome"]
blacklist_file <- snakemake@output[[1]]

# Download blacklist xlsx file
download.file(
  "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-023-03027-3/MediaObjects/13059_2023_3027_MOESM2_ESM.xlsx",
  "resources/blacklists.xlsx"
)

# Load blacklist xlsx file
sheet_names <- getSheetNames("resources/blacklists.xlsx")
blacklists <- list()
for (i in sheet_names) {
  blacklist <- read.xlsx(
    "resources/blacklists.xlsx",
    sheet = i,
    colNames = FALSE,
  )
  blacklists[[i]] <- blacklist
}

# Select the blacklist for the genome
if (genome %in% sheet_names) {
  print(paste("Fetching", genome, "Cut & Run blacklist..."))
  blacklist <- blacklists[[genome]]
  write.table(
    blacklist,
    blacklist_file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
} else {
  writeLines(paste(
    "Cut & Run blacklist not available for",
    genome,
    "\nFetching ChIP-seq blacklist instead..."
  ))
  url <- snakemake@params["url"]
  download.file(url, blacklist_file)
}
