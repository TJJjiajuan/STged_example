## Load necessary libraries
library(Matrix)
library(readr)
library(dplyr)
library(readxl)

# Read barcode data
file_path <- "GSE113576_barcodes.tsv.gz"
barcodes_data <- read_tsv(file_path, col_names = FALSE, show_col_types = FALSE)

# Load expression matrix (ensure the correct file is being used)
expression_matrix_file <- "GSE113576_matrix.mtx.gz"  # Ensure you have the correct matrix file
expression_matrix <- tryCatch(
  readMM(expression_matrix_file),  # Use the correct file for expression data
  error = function(e) {
    message("An error occurred while reading the file: ", e)
    return(NULL)
  }
)

# Check if the matrix was read successfully
if (!is.null(expression_matrix)) {
  print(dim(expression_matrix))  # Print the dimensions of the matrix
  print(expression_matrix[1:5, 1:5])  # View a small subset of the matrix
} else {
  message("Failed to load the matrix.")
}

# Load the Excel data for cell labels
file_path <- "NIHMS1024025-supplement-Table_S1.xlsx"
lable_data <- readxl::read_excel(file_path)

# Rename columns based on the first row, then remove the first row
colnames(lable_data) <- lable_data[1, ]
lable_data <- lable_data[-1, ]

# Extract and sort unique cell types
cell_type <- sort(unique(lable_data$`Cell class (determined from clustering of all cells)`))

# Exclude specific cell types
sec_lable <- cell_type[-c(1, 6, 9, 13, 14)]  # Adjust indices to your needs

# Filter the relevant rows from lable_data and select specific columns
lable_data_sec <- lable_data[lable_data$`Cell class (determined from clustering of all cells)` %in% sec_lable, c(1, 3, 4)]

# Rename columns
colnames(lable_data_sec) <- c("name", "replicate", "class")

# Check the distribution of the cell types for replicate "1"
table(lable_data_sec[lable_data_sec$replicate %in% "1", ]$class)

# Recoding cell types with standardized names
lable_data_sec <- lable_data_sec %>%
  mutate(class = recode(class,
                        "Astrocytes" = "Astrocyte",
                        "Immature oligodendrocyte" = "OD Immature",
                        "Mature oligodendrocyte" = "OD Mature",
                        " Mural" = "Pericyte"))

# Downsample each cell type to a maximum of 3,000 cells
cell_type_col <- "class"

sampled_data <- lable_data_sec %>%
  group_by(!!sym(cell_type_col)) %>%
  mutate(n_cells = n()) %>%
  filter(n_cells <= 3000 | row_number() <= 3000) %>%
  ungroup() %>%
  select(-n_cells)

# Check the distribution of the sampled cell types
table(sampled_data$class)
print(sampled_data)

# Filter the barcodes based on sampled data
barcode_Sec <- barcodes_data$X1 %in% sampled_data$name

# Subset the expression matrix using the filtered barcodes
scRNA_data <- expression_matrix[, barcode_Sec]

# Print the result
print(dim(scRNA_data))
print(scRNA_data[1:5, 1:5])  # View a subset of the filtered matrix
 

scRNA = list(scexp = scRNA_data,
             sclabel = sampled_data[,-2])
saveRDS(scRNA,"scRNA_MPOA.RDS")