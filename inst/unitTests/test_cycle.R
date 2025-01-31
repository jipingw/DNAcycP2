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
