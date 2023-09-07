library(tidyr)

resmarkers <- read.table("resistance_marker_module/resmarker_table_2.txt", header=T)

resmarkers$resmarker <- paste(resmarkers$Gene, resmarkers$CodonID, sep = "_") #resmarkers names
resmarkers$resmarker_sampleID <- paste(resmarkers$SampleID, resmarkers$resmarker, sep = "_") #column to check for multiple markers in a single sample
resmarkers$contents <- paste(resmarkers$AA, " [", resmarkers$Reads, "]", sep = "")

# Create a factor variable for AARefAlt with the desired order
resmarkers$AARefAlt <- factor(resmarkers$AARefAlt, levels = c("REF", "ALT"), ordered = TRUE)

# Sort resmarkers by resmarker_sampleID and AARefAlt
resmarkers <- resmarkers[order(resmarkers$resmarker_sampleID, resmarkers$AARefAlt, decreasing = F), ]

#aggregate contents by resmarker_sampleID 
agg_contents <- aggregate(contents ~ resmarker_sampleID, data = resmarkers, FUN = paste, collapse = ", ")

# Split the agg_contents$resmarker_sampleID by the penultimate "_" into SampleID and resmarker columns
split_cols <- strsplit(agg_contents$resmarker_sampleID, "_")
agg_contents$SampleID <- sapply(split_cols, function(x) paste(x[1:(length(x) - 2)], collapse = "_"))
agg_contents$resmarker <- sapply(split_cols, function(x) paste(x[(length(x) - 1):length(x)], collapse = "_"))
sampleIDs_unique<-unique(agg_contents$SampleID)
resmarker_unique<-unique(agg_contents$resmarker)

#create df and fill it with the aggregated contets
new_df <- as.data.frame(matrix(agg_contents$contents, nrow = length(sampleIDs_unique), ncol = length(resmarker_unique), byrow = TRUE))

# Assign row and column names
rownames(new_df) <- sampleIDs_unique
colnames(new_df) <- resmarker_unique

# Extract the column names and split into name and number
col_names <- colnames(new_df)
col_names_split <- strsplit(col_names, "_")
col_names_name <- sapply(col_names_split, function(x) x[1])
col_names_number <- as.integer(sapply(col_names_split, function(x) x[2]))

# Sort the column names by name and then by number
sorted_col_names <- col_names[order(col_names_name, col_names_number)]

# Reorder the columns of new_df
new_df <- new_df[, sorted_col_names]