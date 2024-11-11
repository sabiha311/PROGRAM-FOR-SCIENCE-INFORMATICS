# Question 1 - Reading the sequence and printing specific letters

# Define the file path
file_path <- "D:/Python 3.12/chr1_GL383518v1_alt.fa"

# Read the sequence, ignoring the header line
sequence <- paste(readLines(file_path)[-1], collapse = "")

# Print the 10th letter (index 9 in 0-indexing)
cat("10th letter of the sequence:", substr(sequence, 10, 10), "\n")

# Print the 758th letter (index 757 in 0-indexing)
cat("758th letter of the sequence:", substr(sequence, 758, 758), "\n")

# Question 2: Create the reverse complement of the sequence

# Create a named vector for complement bases
complement <- c(A = "T", T = "A", C = "G", G = "C", 
                a = "t", t = "a", c = "g", g = "c")

# Create the reverse complement of the sequence
reverse_complement <- rev(strsplit(sequence, NULL)[[1]])
reverse_complement <- sapply(reverse_complement, function(base) complement[base])
reverse_complement <- paste(reverse_complement, collapse = "")

# Print the 79th letter of the reverse complement (index 78 in 0-indexing)
cat(sprintf("79th letter of the reverse complement: %s\n", substr(reverse_complement, 79, 79)))

# Print the 500th through the 800th letters of the reverse complement (index 499 to 799)
cat(sprintf("500th to 800th letters of the reverse complement: %s\n", substr(reverse_complement, 500, 800)))

# Question 3: Count nucleotides in each kilobase

# Initialize a list to hold the counts per kilobase
my_list <- list()

# Iterate through the sequence in steps of 1000 (kilobases)
for (i in seq(1, nchar(sequence), by = 1000)) {
  # Define the current kilobase
  kilobase_range <- substr(sequence, i, min(i + 999, nchar(sequence)))
  
  # Count nucleotides
  nucleotide_count <- table(strsplit(kilobase_range, NULL)[[1]])
  
  # Store counts in the list
  my_list[[length(my_list) + 1]] <- c(A = as.numeric(nucleotide_count["A"]),
                                      T = as.numeric(nucleotide_count["T"]),
                                      C = as.numeric(nucleotide_count["C"]),
                                      G = as.numeric(nucleotide_count["G"]))
}

# Print the counts for the kilobases
for (i in seq_along(my_list)) {
  cat(sprintf("Nucleotide counts for kilobase starting at position %d: %s\n", 
              (i - 1) * 1000 + 1, toString(my_list[[i]])))
}

# 4a: Counts in the first 1000 base pairs

# Assuming 'my_list' is already created from Question 3
first_kb_list <- my_list[[1]]  # Get the first kilobase counts
cat(sprintf("Counts in the first 1000 base pairs: %s\n", toString(first_kb_list)))

# 4b: Create a list for each kilobase

# Initialize an empty list to store counts for each kilobase
all_kb_lists <- list()

# Loop through each kilobase in my_list
for (i in seq_along(my_list)) {
  kb_list <- my_list[[i]]
  all_kb_lists[[i]] <- kb_list  # Store the counts for each kilobase
}

# Print the counts for each kilobase
for (i in seq_along(all_kb_lists)) {
  cat(sprintf("Kilobase %d counts: %s\n", i, toString(all_kb_lists[[i]])))
}

# 4c: Calculate the sum of each list

# Calculate the sum for each kilobase list
list_sums <- sapply(all_kb_lists, sum)

# Expected sum for each list (since each kilobase should contain 1000 bases)
expected_sum <- 1000

# Print each list and its sum, and compare to the expected value
for (i in seq_along(all_kb_lists)) {
  cat(sprintf("Kilobase %d: %s, Sum = %d, Expected = %d\n", 
              i, toString(all_kb_lists[[i]]), list_sums[i], expected_sum))
}


# 4d: Check for any lists whose sums are not equal to the expected value

# Calculate the sum for each kilobase list
list_sums <- sapply(all_kb_lists, sum)

# Expected sum for each list (since each kilobase should contain 1000 bases)
expected_sum <- 1000

# Check for any lists whose sums are not equal to the expected value
unequal_sums <- which(list_sums != expected_sum)

if (length(unequal_sums) > 0) {
  cat(sprintf("There are %d kilobases where the sum is not equal to 1000: %s\n", 
              length(unequal_sums), toString(unequal_sums)))
} else {
  cat("All kilobases have the expected sum of 1000.\n")
}
#4e :1. # Q: What is the expected sum for each list?
# A: The expected sum for each list is 1000, as each kilobase should ideally contain 1000 base pairs.

# 2Q: Are there any lists whose sums are not equal to the expected value?
if (length(unequal_sums) > 0) {
  cat(sprintf("There are %d kilobases where the sum is not equal to 1000: %s\n", 
              length(unequal_sums), toString(unequal_sums)))
} else {
  cat("All kilobases have the expected sum of 1000.\n")
}

# 3Q: Explanation for differences in expected and observed results
# A: Differences in expected and observed sums can arise due to:
#    - Non-standard characters in the sequence (e.g., 'N' for unknown bases), which are not counted in A, C, G, or T.
#    - Partial kilobases: If the total sequence length is not a multiple of 1000, the final kilobase may contain fewer than 1000 bases.
#    - These situations would result in kilobase sums that are less than the expected 1000.                            
