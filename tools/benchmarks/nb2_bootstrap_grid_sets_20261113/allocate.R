# Allocate pending grid cases to two nonoverlapping single-thread jobs.
# Rscript allocate.R ORIGINAL_GRID ALLOCATION_DIR
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
original <- normalizePath(args[1L], mustWork = TRUE)
out <- normalizePath(args[2L], mustWork = FALSE)
stopifnot(!dir.exists(file.path(original, "run.lock")))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
config <- readRDS(file.path(original, "configuration.rds"))
stopifnot(identical(config$panel, "grid"))
plan <- config$plan
done <- file.exists(file.path(original, "cases", paste0(plan$id, ".rds")))
remaining <- plan[!done, , drop = FALSE]
stopifnot(nrow(remaining) >= 2L)
# Preserve the frozen plan order within each shard. Alternating cases
# balances the mix of scenarios and keeps nearby waves close together.
left <- seq_len(nrow(remaining)) %% 2L == 1L
remaining$shard <- ifelse(left, "a", "b")
writeLines(remaining$id[left], file.path(out, "a.txt"))
writeLines(remaining$id[!left], file.path(out, "b.txt"))
utils::write.csv(remaining, file.path(out, "allocation.csv"), row.names = FALSE)
print(with(remaining, table(shard, scenario)))
