#!/usr/bin/env Rscript
# Reproduces the Skillings--Mack analyses reported in Table 5 and Figure 3.
# Usage: Rscript analysis/reproduce.R [results.csv] [output-directory]

required_packages <- c("NSM3")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages)) {
    stop("Missing R package(s): ", paste(missing_packages, collapse = ", "),
         ". Install them with: install.packages(c(\"",
         paste(missing_packages, collapse = "\", \""), "\"))", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[[1]] else "data/instance_generator/all_data.csv"
outdir <- if (length(args) >= 2) args[[2]] else "analysis/output"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

factor_columns <- c("grid", "slope", "wind", "delay", "num_resources",
                    "num_decision_points", "first_release_time", "last_release_time")
required_columns <- c("instance", factor_columns, "algorithm", "seed", "objv", "lb")
dat <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns)) {
    stop("Input is missing required column(s): ", paste(missing_columns, collapse = ", "), call. = FALSE)
}
# The default configuration is shared by the four parameter groups. Keep the
# first recorded result for a repeated (instance, algorithm, seed), matching
# the normalization used for the paper analyses.
dat <- dat[!duplicated(dat[c("instance", "algorithm", "seed")]), , drop = FALSE]

# The four parameter groups in the manuscript. Factors not being varied are
# fixed at their stated default levels.
defaults <- c(grid = "Medium", num_decision_points = "Moderate", delay = "High",
              num_resources = "Moderate", slope = "Moderate", wind = "Light",
              first_release_time = "Early", last_release_time = "VeryLate")
groups <- list(
    instance_size = c("grid", "num_decision_points"),
    suppression_capacity = c("delay", "num_resources"),
    environmental_factors = c("slope", "wind"),
    release_window = c("first_release_time", "last_release_time")
)

instance_seed <- sub(".*_", "", dat$instance)

subset_group <- function(varying) {
    fixed <- setdiff(names(defaults), varying)
    keep <- rep(TRUE, nrow(dat))
    for (name in fixed) keep <- keep & dat[[name]] == defaults[[name]]
    dat[keep, , drop = FALSE]
}

# Returns the omnibus test, treatment scores, and all post-hoc comparisons.
# Each block-treatment cell has the same number of replications in this study.
skillings_mack <- function(df, replication, treatment, block, value, alpha = 0.001) {
    df <- df[order(df[[block]], df[[treatment]], df[[replication]]), , drop = FALSE]
    replication_levels <- unique(df[[replication]])
    treatment_levels <- unique(df[[treatment]])
    block_levels <- unique(df[[block]])
    c_reps <- length(replication_levels)
    k <- length(treatment_levels)
    n <- length(block_levels)

    expected <- c_reps * k * n
    if (nrow(df) != expected || any(table(df[[block]], df[[treatment]]) != c_reps)) {
        stop("The design is not balanced for this analysis.", call. = FALSE)
    }

    values <- array(df[[value]], dim = c(c_reps, k, n))
    # The original analysis used NSM3's default Monte-Carlo procedure.
    omnibus <- NSM3::pMackSkil(values)

    ranks <- ave(df[[value]], df[[block]], FUN = function(x) rank(x, ties.method = "average"))
    mean_ranks <- aggregate(ranks, list(block = df[[block]], treatment = df[[treatment]]), mean)
    names(mean_ranks)[3] <- "mean_rank"
    scores <- aggregate(mean_rank ~ treatment, mean_ranks, sum)
    scores <- scores[order(scores$mean_rank), ]
    names(scores) <- c("treatment", "score")

    delta <- sqrt(k * (nrow(df) + n) / 12) * NSM3::cRangeNor(alpha, k)
    comparisons <- expand.grid(treatment_a = scores$treatment, treatment_b = scores$treatment,
                               stringsAsFactors = FALSE)
    comparisons <- comparisons[comparisons$treatment_a < comparisons$treatment_b, ]
    score_map <- setNames(scores$score, scores$treatment)
    comparisons$difference <- abs(score_map[comparisons$treatment_a] - score_map[comparisons$treatment_b])
    comparisons$significant <- comparisons$difference > delta

    list(p_value = omnibus, scores = scores, delta = delta, comparisons = comparisons,
         replications = c_reps, treatments = k, blocks = n)
}

# NSM3's default omnibus procedure is Monte Carlo; fix its seed so that the
# reported p-values are reproducible.
set.seed(202601)

summary_rows <- list()
algorithm_scores <- list()
algorithm_pairs <- list()
difficulty_scores <- list()
difficulty_pairs <- list()

for (group_name in names(groups)) {
    varying <- groups[[group_name]]
    group_data <- subset_group(varying)
    group_data$instance_seed <- instance_seed[match(group_data$instance, dat$instance)]
    group_data$algorithm_seed <- group_data$seed
    group_data$factor_pair <- do.call(paste, c(group_data[varying], sep = " | "))

    # Table 5: algorithms are treatments; each instance is a block.
    algorithm_data <- group_data
    algorithm_data$block <- paste(algorithm_data$instance_seed,
                                  do.call(paste, c(algorithm_data[varying], sep = " | ")),
                                  sep = " :: ")
    result <- skillings_mack(algorithm_data, "algorithm_seed", "algorithm", "block", "objv")
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
        analysis = "algorithms", group = group_name, p_value = result$p_value$p.val,
        delta = result$delta, replications = result$replications,
        treatments = result$treatments, blocks = result$blocks)
    algorithm_scores[[group_name]] <- transform(result$scores, group = group_name)
    algorithm_pairs[[group_name]] <- transform(result$comparisons, group = group_name, delta = result$delta)

    # Figure 3: factor combinations are treatments; algorithm-instance pairs
    # are blocks; response is relative deviation from the best-known value.
    best <- aggregate(objv ~ instance, group_data, min)
    names(best)[2] <- "bkv"
    difficulty_data <- merge(group_data, best, by = "instance", sort = FALSE)
    difficulty_data$relative_deviation <- 100 * (difficulty_data$objv / difficulty_data$bkv - 1)
    difficulty_data$block <- paste(difficulty_data$algorithm, difficulty_data$instance_seed, sep = " :: ")
    result <- skillings_mack(difficulty_data, "algorithm_seed", "factor_pair", "block", "relative_deviation")
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
        analysis = "difficulty", group = group_name, p_value = result$p_value$p.val,
        delta = result$delta, replications = result$replications,
        treatments = result$treatments, blocks = result$blocks)
    difficulty_scores[[group_name]] <- transform(result$scores, group = group_name)
    difficulty_pairs[[group_name]] <- transform(result$comparisons, group = group_name, delta = result$delta)
}

write.csv(do.call(rbind, summary_rows), file.path(outdir, "skillings_mack_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, algorithm_scores), file.path(outdir, "table5_scores.csv"), row.names = FALSE)
write.csv(do.call(rbind, algorithm_pairs), file.path(outdir, "table5_comparisons.csv"), row.names = FALSE)
write.csv(do.call(rbind, difficulty_scores), file.path(outdir, "difficulty_scores.csv"), row.names = FALSE)
write.csv(do.call(rbind, difficulty_pairs), file.path(outdir, "difficulty_comparisons.csv"), row.names = FALSE)

cat("Wrote Skillings--Mack results to", normalizePath(outdir), "\n")
