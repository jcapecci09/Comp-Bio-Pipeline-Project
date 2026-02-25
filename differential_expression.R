library(sleuth)
library(dplyr)

# Read in table
stab = read.table('table.txt', header=TRUE)

# Make sure condition is a factor
stab$condition <- as.factor(stab$condition)

# intialize sleuth object
so = sleuth_prep(stab)

# fit model comapring conditions
so = sleuth_fit(so, ~condition, 'full')

# fit reduced model to compare in the likelihood ratio test
so = sleuth_fit(so, ~1, 'reduced')

# Perform likelihood ratio test for different expression
so = sleuth_lrt(so, 'reduced', 'full')

# extract test results
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all=FALSE)

# find sigfnficant results
sleuth_significant = dplyr::filter(sleuth_table, qval <=0.05) |> dplyr::arrange(pval)

# Filter out unwated columns
sleuth_results = sleuth_significant[, c('target_id', 'test_stat', 'pval', 'qval')]

# Write to result text file
write.table(sleuth_results, file='sleuth_results.txt', sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
