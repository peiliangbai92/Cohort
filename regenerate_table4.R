### Regenerate Table 4 LaTeX rows from the aggregated summary.
res <- readRDS("logs/summary_table.rds")

format_row <- function(name, sub) {
    if (nrow(sub) == 1) {
        sprintf("  %s & %d & %.3f & %.3f & %.3f & %.0f \\\\",
                name, sub$cp, sub$truth, sub$mean, sub$sd, sub$pct)
    } else {
        ## multi-CP: use multirow
        rows <- character()
        for (i in seq_len(nrow(sub))) {
            if (i == 1) {
                prefix <- sprintf("  \\multirow{%d}*{%s}", nrow(sub), name)
            } else {
                prefix <- "                  "
            }
            rows <- c(rows, sprintf("%s & %d & %.3f & %.3f & %.3f & %.0f \\\\",
                                    prefix, sub$cp[i], sub$truth[i],
                                    sub$mean[i], sub$sd[i], sub$pct[i]))
        }
        paste(rows, collapse = "\n")
    }
}

scenarios <- unique(res$scenario)
single_block <- character()
multi_block <- character()
for (s in scenarios) {
    sub <- res[res$scenario == s, ]
    line <- format_row(s, sub)
    if (grepl("^S", s)) {
        single_block <- c(single_block, line)
    } else {
        multi_block <- c(multi_block, line)
    }
}

cat("=== Single-CP block (S1-S7) ===\n")
cat(paste(single_block, collapse = "\n"), "\n")
cat("\n=== Multi-CP block (M1-M5) ===\n")
cat(paste(multi_block, collapse = "\n"), "\n")
