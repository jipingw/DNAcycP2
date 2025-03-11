library(tinytest)
seqs = c(
    "AAAAATTTTTAAAAATTTTTAAAAATTTTTAAAAATTTTTAAAAATTTTT",
    "AAAAACCCCCAAAAACCCCCAAAAACCCCCAAAAACCCCCAAAAACCCCC")

preds <- DNAcycP2::cycle(seqs, smooth = TRUE)

expect_equal(DNAcycP2::cycle(seqs[1], smooth=TRUE)[[1]]$C0S_norm,
            -0.4635092, tolerance=1e-4)
expect_equal(DNAcycP2::cycle(seqs[1], smooth=TRUE)[[1]]$C0S_unnorm,
            -0.3132579, tolerance=1e-4)

expect_equal(DNAcycP2::cycle(seqs[2], smooth=FALSE)[[1]]$C0_norm,
             3.776794, tolerance=1e-4)
expect_equal(DNAcycP2::cycle(seqs[2], smooth=FALSE)[[1]]$C0_unnorm,
             1.656955, tolerance=1e-4)

# Test `cycle_fasta()` function
temp_file <- tempfile(fileext = ".fasta")
writeLines(">1\nACTGCTAGTCACTGCTAGTCACTGCTAGTCACTGCTAGTCACTGCTAGTC", temp_file)
expect_equal(cycle_fasta(temp_file, smooth=TRUE)[[1]]$C0S_unnorm,
    -0.03295413,
    tolerance = 1e-4
)
expect_equal(cycle_fasta(temp_file, smooth=FALSE)[[1]]$C0_unnorm,
    -0.1674057,
    tolerance = 1e-4
)

# Test `cycle()` function with `save_path_prefix` argument
save_path_prefix <- tempfile()
DNAcycP2::cycle(seqs, smooth=TRUE, save_path_prefix=save_path_prefix)
expect_true(file.exists(paste0(save_path_prefix, "_C0S_norm.txt")))
expect_true(file.exists(paste0(save_path_prefix, "_C0S_unnorm.txt")))

DNAcycP2::cycle(seqs, smooth=FALSE, save_path_prefix=save_path_prefix)
expect_true(file.exists(paste0(save_path_prefix, "_C0_norm.txt")))
expect_true(file.exists(paste0(save_path_prefix, "_C0_unnorm.txt")))
