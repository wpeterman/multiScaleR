test_that("shipped Bayesian example objects carry no absolute build paths", {
  # Provenance records name source files. They must be package-relative so the
  # released data never embed a developer's home directory.
  files <- system.file("extdata", c("bayesian-annual-example.rds",
    "bayesian-prediction-example.rds", "bayesian-computation-example.rds"),
    package = "multiScaleR")
  files <- files[nzchar(files)]
  skip_if(!length(files), "example objects are not installed")
  for (file in files) {
    bytes <- serialize(readRDS(file), NULL, ascii = TRUE)
    text <- rawToChar(bytes[bytes != as.raw(0)])
    expect_false(grepl("[A-Za-z]:[/\\]+Users[/\\]|/home/|/Users/", text),
                 info = basename(file))
  }
  annual <- readRDS(files[basename(files) == "bayesian-annual-example.rds"])
  expect_true(all(grepl("^(R|src)/", names(annual$provenance$engine_fingerprint))))
})
