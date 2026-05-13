# Title:   Script for analyzing the Rocker DMS data from Liposome screen
# Author:  Fatlum Hajredini 
# Contact: fatlum.hajredini@gmail.com
# Date:    09/05/2026

# setwd("/Volumes/DATA/Liposome_final/Data")

# --- Load Data ---
sort <- read.csv(file = "Rocker_sorted.txt", head = FALSE)    # Read data from sort pool
unsort <- read.csv(file = "Rocker_unsorted.txt", head = FALSE) # Read data from unsorted pool

unique_seqs <- unique(rbind(sort, unsort))

# --- Filter out stop codons ---
# Find indexes of sequences with internal stop codons 
contains_stop <- c() 

for (i in unique_seqs[, 1]) {
  hold_seq <- unlist(strsplit(i, ""))
  # Check if sequences contain stop codons besides a potentially terminal one
  contains_stop <- append(contains_stop, sum(hold_seq == "*") > 1)
}

# Remove sequences with internal stop codons 
unique_seqs_filtered <- unique_seqs[!contains_stop, ]

# --- Count unique sequences in the sorted and unsorted pools ---
count_sort <- c()
count_unsort <- c()

for (i in unique_seqs_filtered) {
  # Compute occurrences in the sort pool
  count_sort <- append(count_sort, sum(sort[, 1] == i)) 
  # Compute occurrences in the unsort pool
  count_unsort <- append(count_unsort, sum(unsort[, 1] == i)) 
}

count_total <- count_sort + count_unsort  # Total number of reads per sequence
p_sort <- count_sort / count_total        # Probability of being sorted

# --- Clean Sequences ---
# Drop starting Met (M) and terminal Stop codon (*)
final_seqs <- substr(unique_seqs_filtered, 2, nchar(unique_seqs_filtered) - 1)

# --- Identify Mutations ---
wt_seq <- unlist(strsplit("YYKEIAHALFSALFALSELYIAVRY", "")) # WT Rocker sequence
index  <- seq(1, 25, 1)                                   # Sequence indexing

n <- length(final_seqs)

# Initialize data vectors
mutation_1    <- character(n) # Identity of the first mutation
mutation_2    <- character(n) # Identity of the second mutation (if present)
index_1       <- numeric(n)   # Location of the first mutation 
index_2       <- numeric(n)   # Location of the second mutation (if present)
num_mutations <- numeric(n)   # Count of mutations per sequence

# Compare sequences to WT to find substitutions and positions
for (i in seq_along(final_seqs)) {
  seq_hold <- unlist(strsplit(final_seqs[i], ""))
  diff_idx <- which(seq_hold != wt_seq)
  mutations <- length(diff_idx)
  num_mutations[i] <- mutations
  
  if (mutations == 1) {
    mutation_1[i] <- seq_hold[diff_idx]
    mutation_2[i] <- "0"
    index_1[i]    <- diff_idx
    index_2[i]    <- 0
  } else if (mutations == 2) {
    mutation_1[i] <- seq_hold[diff_idx[1]]
    mutation_2[i] <- seq_hold[diff_idx[2]]
    index_1[i]    <- diff_idx[1]
    index_2[i]    <- diff_idx[2]
  }
}

# --- Consolidate and Export Data ---
experiment_all <- data.frame(
  "sequence"       = final_seqs, 
  "count_sorted"   = count_sort, 
  "count_unsorted" = count_unsort, 
  "total_reads"    = count_total, 
  "p_sort"         = p_sort,
  "n_mutations"    = num_mutations,
  "first_mutation" = mutation_1,
  "second_mutation"= mutation_2,
  "first_index"    = index_1,
  "second_index"   = index_2
)

# Order by total reads (descending)
experiment_all <- experiment_all[order(experiment_all$total_reads, decreasing = TRUE), ]

# Export results
write.csv(experiment_all, file = "Rocker_DMS.csv", quote = FALSE)

# --- Categorical Filtering ---
single_mutants <- experiment_all[experiment_all$n_mutations == 1, ]
double_mutants <- experiment_all[experiment_all$n_mutations == 2, ]

# Filter double mutants by Proline (P) substitutions
proline_filter <- (double_mutants$first_mutation == "P") + (double_mutants$second_mutation == "P")

double_mutants_noP     <- double_mutants[proline_filter == 0, ] # No Prolines 
double_mutants_singleP <- double_mutants[proline_filter == 1, ] # One Proline
double_mutants_doubleP <- double_mutants[proline_filter == 2, ] # Two Prolines

# --- Visualization ---

# Figure 4D: Double mutants vs Single mutants
plot(log(double_mutants_noP$total_reads), double_mutants_noP$p_sort, pch = 19, cex = 0.5)
points(log(single_mutants$total_reads), single_mutants$p_sort, pch = 19, cex = 0.5, col = "darkgreen")

# Figure 6A: Impact of Proline substitutions on double mutants
plot(log(double_mutants_noP$total_reads), double_mutants_noP$p_sort, pch = 19, cex = 0.5)
points(log(double_mutants_singleP$total_reads), double_mutants_singleP$p_sort, pch = 19, cex = 0.5, col = "red")
points(log(double_mutants_doubleP$total_reads), double_mutants_doubleP$p_sort, pch = 19, cex = 0.5, col = "blue")

# --- Range Filtering (100-200 reads) ---
# Logic: sum == 2 implies both conditions (min and max) are met
doubleNoP_100_200     <- (double_mutants_noP$total_reads > 100) + (double_mutants_noP$total_reads < 200)
double_singleP_100_200 <- (double_mutants_singleP$total_reads > 100) + (double_mutants_singleP$total_reads < 200)
double_doubleP_100_200 <- (double_mutants_doubleP$total_reads > 100) + (double_mutants_doubleP$total_reads < 200)

# Figure 6B: Histograms of sorting odds
par(mfrow = c(1, 3))
hist(double_mutants_noP$p_sort[doubleNoP_100_200 == 2], col = 'black', xlim = c(0, 1), breaks = 10, ylim = c(0, 60))
hist(double_mutants_singleP$p_sort[double_singleP_100_200 == 2], col = 'red', ylim = c(0, 200))
hist(double_mutants_doubleP$p_sort[double_doubleP_100_200 == 2], col = "blue", xlim = c(0, 1), breaks = 5)

# --- Proline Position Analysis ---
filter_SingleP <- double_mutants_singleP[double_singleP_100_200 == 2, ]

# Initialize with base indices to ensure barplot shows all 25 positions
position_greater <- index 
position_less    <- index

# Categorize proline positions by sorting odds (thresh = 0.5)
for (i in seq_along(filter_SingleP$sequence)) {
  if (filter_SingleP$p_sort[i] >= 0.5) {
    if (filter_SingleP$first_mutation[i] == "P") {
      position_greater <- append(position_greater, filter_SingleP$first_index[i])
    } else {
      position_greater <- append(position_greater, filter_SingleP$second_index[i])
    }
  } else {
    if (filter_SingleP$first_mutation[i] == "P") {
      position_less <- append(position_less, filter_SingleP$first_index[i])
    } else {
      position_less <- append(position_less, filter_SingleP$second_index[i])
    }
  }
}

# Figure 6C: Position distribution barplots
par(mfrow = c(2, 1))
barplot(table(position_greater) - 1, main = "Sorting Odds >= 0.5") # Subtract 1 to remove artificial base entries
barplot(table(position_less) - 1, main = "Sorting Odds < 0.5")