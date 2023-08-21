# ----------------------------------------------------------------------------
# Copyright (c) 2022--, Panpiper development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]

# The unitig file to be filtered
file <- read.delim(input, header = TRUE, check.names = FALSE)

# Perform Holm's correction on the p-values
file$pv_adj <- p.adjust(file$`filter-pvalue`, method = "holm")

# Filter unitigs based on Holm's corrected p-values
threshold <- 0.1  # Set the significance level (alpha)
file <- file[file$pv_adj <= threshold, ]
file <- file[order(file$pv_adj), ]

# Write the filtered unitig file with adjusted p-values
write.table(file, output, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")