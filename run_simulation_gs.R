rm(list = ls())
setwd('C://Users/ya-ping.lin/Documents/IMIN_GS/GS/simulation/')

# ============================================================
# run_simulation_gs_robust_MTfixed_v2.R
# Robust simulation + GS script with:
#   - fixed multi-trait CV (mask only target trait in test lines)
#   - tolerant MT GEBV extraction/alignment
#   - olive boxplot colors
#   - checkpoint saves per scenario
#   - tryCatch so failed replicates do not stop the full run
# ============================================================

source("../gs_pipeline.R")

library(MASS)
library(data.table)
library(parallel)
library(doParallel)
library(foreach)
library(ggplot2)
library(sommer)
library(SNPRelate)
library(gdsfmt)

vcf_file <- "../../MMC_MAF005.vcf"

trait_names <- c("seed_weight", "seed_per_pod", "pod_per_plant")
h2_vec <- c(0.7, 0.5, 0.3)

corr_scenarios <- c(0.1, 0.3, 0.5)

nrep <- 10
kfold <- 5
cv_reps <- 3
causal_prop <- 0.1
ncores <- max(1, parallel::detectCores() - 2)

nugget <- 1e-4
tolParInv <- 1e-2

outdir <- "simulation_results"

olive_fill <- c("Single-trait" = "darkolivegreen4", "Multi-trait" = "darkolivegreen2")

build_X_from_vcf <- function(vcf_fn,
                             gds_fn = sub("\\.vcf(\\.gz)?$", ".gds", vcf_fn, ignore.case = TRUE),
                             method_vcf = "biallelic.only",
                             snpfirstdim = FALSE,
                             verbose = TRUE) {

  gs_load_packages()

  if (!file.exists(vcf_fn)) stop("VCF not found: ", vcf_fn)

  if (!file.exists(gds_fn)) {
    SNPRelate::snpgdsVCF2GDS(vcf_fn, gds_fn, method = method_vcf, verbose = verbose)
  }

  genofile <- SNPRelate::snpgdsOpen(gds_fn)
  on.exit(SNPRelate::snpgdsClose(genofile))

  sample.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id"))
  snp.id    <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id"))

  X <- SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = snpfirstdim)

  rownames(X) <- sample.id
  colnames(X) <- snp.id

  X[X == 3] <- NA
  X
}

simulate_3traits <- function(X,
                             h2 = c(0.7, 0.5, 0.3),
                             trait_names = c("seed_weight", "seed_per_pod", "pod_per_plant"),
                             rg = -0.3,
                             causal_prop = 0.1,
                             seed = 123) {
  set.seed(seed)

  Xs <- scale(X)
  n <- nrow(Xs)
  m <- ncol(Xs)
  k <- length(h2)

  Rg <- matrix(rg, nrow = k, ncol = k)
  diag(Rg) <- 1

  ev <- eigen(Rg, symmetric = TRUE)$values
  if (any(ev <= 0)) stop("Genetic correlation matrix is not positive definite.")

  n_causal <- max(1, round(m * causal_prop))
  causal_idx <- sample(seq_len(m), n_causal)

  B <- matrix(0, nrow = m, ncol = k)
  B[causal_idx, ] <- MASS::mvrnorm(n = n_causal, mu = rep(0, k), Sigma = Rg)

  G_true <- Xs %*% B
  G_true <- scale(G_true)

  E <- matrix(NA, nrow = n, ncol = k)
  for (j in seq_len(k)) {
    var_e <- (1 - h2[j]) / h2[j]
    E[, j] <- rnorm(n, mean = 0, sd = sqrt(var_e))
  }

  Y <- G_true + E

  colnames(Y) <- trait_names
  colnames(G_true) <- trait_names
  colnames(E) <- trait_names
  colnames(B) <- trait_names

  list(
    Y = as.data.frame(Y),
    G_true = as.data.frame(G_true),
    E = as.data.frame(E),
    B = B,
    causal_idx = causal_idx,
    realized_h2 = sapply(seq_len(k), function(j) var(G_true[, j]) / var(Y[, j])),
    realized_rg = cor(G_true)
  )
}

make_cv_folds <- function(ids, k = 5, seed = 1) {
  set.seed(seed)
  ids <- sample(ids)
  split(ids, cut(seq_along(ids), breaks = k, labels = FALSE))
}

# tolerant extractor matching CV/full-fit variations
extract_gebv_from_fit <- function(fit, traits, fallback_ids) {
  Uterm <- fit$U$`u:id`
  fallback_ids <- as.character(fallback_ids)

  if (is.null(Uterm)) stop("fit$U$`u:id` is NULL.")

  nm <- names(Uterm)

  if (!is.null(nm) && all(traits %in% nm)) {
    out <- sapply(traits, function(tr) {
      u <- Uterm[[tr]]
      vals <- if (is.matrix(u)) u[, 1] else as.numeric(u)
      ids  <- rownames(u)
      if (is.null(ids) || length(ids) == 0) ids <- names(u)
      if (is.null(ids) || length(ids) == 0) ids <- fallback_ids
      vals[match(fallback_ids, ids)]
    })
    out <- as.matrix(out)
    rownames(out) <- fallback_ids
    colnames(out) <- traits
    return(out)
  }

  if (!is.null(nm)) {
    clean <- function(x) gsub("[^A-Za-z0-9]", "", x)
    nm_clean <- clean(nm)
    tr_clean <- clean(traits)
    idx <- match(tr_clean, nm_clean)

    if (all(!is.na(idx))) {
      out <- sapply(seq_along(traits), function(i) {
        u <- Uterm[[idx[i]]]
        vals <- if (is.matrix(u)) u[, 1] else as.numeric(u)
        ids  <- rownames(u)
        if (is.null(ids) || length(ids) == 0) ids <- names(u)
        if (is.null(ids) || length(ids) == 0) ids <- fallback_ids
        vals[match(fallback_ids, ids)]
      })
      out <- as.matrix(out)
      rownames(out) <- fallback_ids
      colnames(out) <- traits
      return(out)
    }
  }

  if (is.list(Uterm) && length(Uterm) == length(traits)) {
    out <- sapply(seq_along(traits), function(i) {
      u <- Uterm[[i]]
      vals <- if (is.matrix(u)) u[, 1] else as.numeric(u)
      ids  <- rownames(u)
      if (is.null(ids) || length(ids) == 0) ids <- names(u)
      if (is.null(ids) || length(ids) == 0) ids <- fallback_ids
      vals[match(fallback_ids, ids)]
    })
    out <- as.matrix(out)
    rownames(out) <- fallback_ids
    colnames(out) <- traits
    return(out)
  }

  if (is.matrix(Uterm) && ncol(Uterm) == length(traits)) {
    out <- Uterm
    if (is.null(rownames(out))) rownames(out) <- fallback_ids
    out <- out[match(fallback_ids, rownames(out)), , drop = FALSE]
    colnames(out) <- traits
    return(out)
  }

  stop("Unsupported structure in fit$U$`u:id`. Names found: ",
       paste(nm, collapse = ", "))
}

run_st_cv_trueacc <- function(pheno, gtrue, G, trait, id_col = "id",
                              k = 5, reps = 3, seed = 123,
                              nugget = 1e-4, tolParInv = 1e-2) {
  ids <- as.character(pheno[[id_col]])
  names(ids) <- ids
  Gs <- G + diag(nugget, nrow(G))
  rownames(Gs) <- colnames(Gs) <- rownames(G)

  acc_all <- c()

  for (r in seq_len(reps)) {
    fold_list <- make_cv_folds(ids, k = k, seed = seed + r)

    for (f in seq_along(fold_list)) {
      test_ids <- fold_list[[f]]

      dat_cv <- pheno
      dat_cv[dat_cv[[id_col]] %in% test_ids, trait] <- NA

      train_mask <- !pheno[[id_col]] %in% test_ids
      mu <- mean(pheno[[trait]][train_mask], na.rm = TRUE)
      sdv <- sd(pheno[[trait]][train_mask], na.rm = TRUE)
      if (is.na(sdv) || sdv == 0) sdv <- 1

      dat_fit <- data.frame(id = as.character(dat_cv[[id_col]]), stringsAsFactors = FALSE)
      dat_fit$y <- (dat_cv[[trait]] - mu) / sdv

      fit <- sommer::mmer(
        y ~ 1,
        random = ~ sommer::vsr(id, Gu = Gs[ids, ids]),
        rcov   = ~ units,
        data   = dat_fit,
        tolParInv = tolParInv
      )

      gebv <- extract_gebv_from_fit(fit, traits = c(trait), fallback_ids = dat_fit$id)
      pred <- gebv[match(test_ids, rownames(gebv)), , drop = FALSE]
      rownames(pred) <- test_ids

      tru <- gtrue[gtrue[[id_col]] %in% test_ids, c(id_col, trait)]
      tru <- tru[match(test_ids, tru[[id_col]]), , drop = FALSE]
      acc_all <- c(acc_all, cor(tru[[trait]], pred[, 1], use = "complete.obs"))
    }
  }

  mean(acc_all, na.rm = TRUE)
}

# fixed MT CV: mask only target trait in test lines
run_mt_cv_trueacc <- function(pheno, gtrue, G, traits, id_col = "id",
                              k = 5, reps = 3, seed = 123,
                              nugget = 1e-4, tolParInv = 1e-2) {
  ids <- as.character(pheno[[id_col]])
  names(ids) <- ids
  Gs <- G + diag(nugget, nrow(G))
  rownames(Gs) <- colnames(Gs) <- rownames(G)

  acc_store <- vector("list", length(traits))
  names(acc_store) <- traits
  for (tr in traits) acc_store[[tr]] <- c()

  for (target_trait in traits) {
    for (r in seq_len(reps)) {
      fold_list <- make_cv_folds(ids, k = k, seed = seed + r + match(target_trait, traits) * 1000)

      for (f in seq_along(fold_list)) {
        test_ids <- fold_list[[f]]

        dat_cv <- pheno
        dat_cv[dat_cv[[id_col]] %in% test_ids, target_trait] <- NA

        train_mask <- !pheno[[id_col]] %in% test_ids
        mu <- sapply(traits, function(tr) mean(pheno[train_mask, tr], na.rm = TRUE))
        sdv <- sapply(traits, function(tr) sd(pheno[train_mask, tr], na.rm = TRUE))
        sdv[is.na(sdv) | sdv == 0] <- 1

        dat_scaled <- dat_cv
        for (tr in traits) {
          dat_scaled[[tr]] <- (dat_scaled[[tr]] - mu[tr]) / sdv[tr]
        }
        dat_scaled[[id_col]] <- as.character(dat_scaled[[id_col]])

        form <- as.formula(paste0("cbind(", paste(traits, collapse = ","), ") ~ 1"))

        fit <- sommer::mmer(
          form,
          random = ~ sommer::vsr(id, Gu = Gs[ids, ids], Gtc = unsm(length(traits))),
          rcov   = ~ sommer::vsr(units, Gtc = unsm(length(traits))),
          data   = data.frame(id = dat_scaled[[id_col]], dat_scaled[, traits, drop = FALSE]),
          tolParInv = tolParInv
        )

        gebv <- extract_gebv_from_fit(fit, traits = traits, fallback_ids = dat_scaled[[id_col]])
        pred <- gebv[match(test_ids, rownames(gebv)), , drop = FALSE]
        rownames(pred) <- test_ids

        tru <- gtrue[gtrue[[id_col]] %in% test_ids, c(id_col, target_trait)]
        tru <- tru[match(test_ids, tru[[id_col]]), , drop = FALSE]

        acc_store[[target_trait]] <- c(
          acc_store[[target_trait]],
          cor(tru[[target_trait]], pred[, target_trait], use = "complete.obs")
        )
      }
    }
  }

  sapply(acc_store, mean, na.rm = TRUE)
}

run_one_sim_rep <- function(rep_id, X, G, ids,
                            rg,
                            h2 = c(0.7, 0.5, 0.3),
                            trait_names = c("seed_weight", "seed_per_pod", "pod_per_plant"),
                            causal_prop = 0.1,
                            k = 5,
                            reps = 3,
                            nugget = 1e-4,
                            tolParInv = 1e-2) {
  sim <- simulate_3traits(
    X = X,
    h2 = h2,
    trait_names = trait_names,
    rg = rg,
    causal_prop = causal_prop,
    seed = 1000 + rep_id
  )

  pheno <- data.frame(id = ids, sim$Y, check.names = FALSE, stringsAsFactors = FALSE)
  gtrue <- data.frame(id = ids, sim$G_true, check.names = FALSE, stringsAsFactors = FALSE)

  out_acc <- list()

  for (tr in trait_names) {
    acc_st <- run_st_cv_trueacc(
      pheno = pheno, gtrue = gtrue, G = G, trait = tr, id_col = "id",
      k = k, reps = reps, seed = 10000 + rep_id,
      nugget = nugget, tolParInv = tolParInv
    )

    out_acc[[paste0("ST_", tr)]] <- data.frame(
      replicate = rep_id,
      scenario = paste0("rg_", rg),
      method = "Single-trait",
      trait = tr,
      accuracy = acc_st
    )
  }

  acc_mt <- run_mt_cv_trueacc(
    pheno = pheno, gtrue = gtrue, G = G, traits = trait_names, id_col = "id",
    k = k, reps = reps, seed = 20000 + rep_id,
    nugget = nugget, tolParInv = tolParInv
  )

  for (tr in trait_names) {
    out_acc[[paste0("MT_", tr)]] <- data.frame(
      replicate = rep_id,
      scenario = paste0("rg_", rg),
      method = "Multi-trait",
      trait = tr,
      accuracy = acc_mt[tr]
    )
  }

  realized <- data.frame(
    replicate = rep_id,
    scenario = paste0("rg_", rg),
    h2_seed_weight = sim$realized_h2[1],
    h2_seed_per_pod = sim$realized_h2[2],
    h2_pod_per_plant = sim$realized_h2[3],
    rg_sw_spp = sim$realized_rg[1, 2],
    rg_sw_ppp = sim$realized_rg[1, 3],
    rg_spp_ppp = sim$realized_rg[2, 3]
  )

  list(
    accuracy = data.table::rbindlist(out_acc),
    realized = realized
  )
}

run_parallel_simulation <- function(X, G, ids,
                                    corr_scenarios = c(-0.1, -0.3, -0.5),
                                    nrep = 100,
                                    h2 = c(0.7, 0.5, 0.3),
                                    trait_names = c("seed_weight", "seed_per_pod", "pod_per_plant"),
                                    causal_prop = 0.1,
                                    k = 5,
                                    reps = 3,
                                    nugget = 1e-4,
                                    tolParInv = 1e-2,
                                    ncores = max(1, parallel::detectCores() - 2),
                                    outdir = "simulation_results") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  on.exit({
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  }, add = TRUE)

  all_acc <- list()
  all_realized <- list()
  all_errors <- list()

  for (rg in corr_scenarios) {
    message("Running scenario: rg = ", rg)

    res <- foreach::foreach(
      i = 1:nrep,
      .packages = c("MASS", "sommer", "data.table"),
      .export = c("simulate_3traits",
                  "make_cv_folds",
                  "extract_gebv_from_fit",
                  "run_st_cv_trueacc",
                  "run_mt_cv_trueacc",
                  "run_one_sim_rep")
    ) %dopar% {
      tryCatch(
        run_one_sim_rep(
          rep_id = i,
          X = X,
          G = G,
          ids = ids,
          rg = rg,
          h2 = h2,
          trait_names = trait_names,
          causal_prop = causal_prop,
          k = k,
          reps = reps,
          nugget = nugget,
          tolParInv = tolParInv
        ),
        error = function(e) {
          list(
            accuracy = data.frame(),
            realized = data.frame(),
            error = data.frame(
              replicate = i,
              scenario = paste0("rg_", rg),
              error = conditionMessage(e),
              stringsAsFactors = FALSE
            )
          )
        }
      )
    }

    acc_rg <- data.table::rbindlist(lapply(res, function(x) x$accuracy), fill = TRUE)
    real_rg <- data.table::rbindlist(lapply(res, function(x) x$realized), fill = TRUE)
    err_rg <- data.table::rbindlist(
      lapply(res, function(x) if (!is.null(x$error)) x$error else data.frame()),
      fill = TRUE
    )

    all_acc[[paste0("rg_", rg)]] <- acc_rg
    all_realized[[paste0("rg_", rg)]] <- real_rg
    all_errors[[paste0("rg_", rg)]] <- err_rg

    tag <- gsub("-", "m", as.character(rg), fixed = TRUE)

    data.table::fwrite(acc_rg, file = file.path(outdir, paste0("GS_accuracy_", tag, ".txt")), sep = "\t")
    data.table::fwrite(real_rg, file = file.path(outdir, paste0("Simulation_realized_", tag, ".txt")), sep = "\t")

    if (nrow(err_rg) > 0) {
      data.table::fwrite(err_rg, file = file.path(outdir, paste0("Errors_", tag, ".txt")), sep = "\t")
    }
  }

  list(
    accuracy = data.table::rbindlist(all_acc, fill = TRUE),
    realized = data.table::rbindlist(all_realized, fill = TRUE),
    errors = data.table::rbindlist(all_errors, fill = TRUE)
  )
}

summarize_accuracy <- function(acc_df) {
  acc_df[, .(
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    sd_accuracy = sd(accuracy, na.rm = TRUE)
  ), by = .(scenario, method, trait)]
}

plot_accuracy_boxplots <- function(acc_df, outdir = "simulation_results") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  p1 <- ggplot(acc_df, aes(x = scenario, y = accuracy, fill = method)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 16) +
    facet_wrap(~ trait, scales = "free_y") +
    scale_fill_manual(values = olive_fill) +
    theme_bw(base_size = 12) +
    labs(
      title = "GS accuracy across correlation scenarios",
      x = "Genetic correlation scenario",
      y = "Prediction accuracy"
    )

  ggsave(file.path(outdir, "Boxplot_accuracy_by_trait.png"),
         p1, width = 10, height = 6, dpi = 300)

  p2 <- ggplot(acc_df, aes(x = method, y = accuracy, fill = method)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 16) +
    facet_grid(trait ~ scenario, scales = "free_y") +
    scale_fill_manual(values = olive_fill) +
    theme_bw(base_size = 12) +
    labs(
      title = "Single-trait vs Multi-trait GS accuracy",
      x = "Method",
      y = "Prediction accuracy"
    )

  ggsave(file.path(outdir, "Boxplot_ST_vs_MT_by_scenario.png"),
         p2, width = 11, height = 7, dpi = 300)

  list(p1 = p1, p2 = p2)
}

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

cat("Building GRM from VCF...\n")
grm_out <- build_grm_from_vcf(vcf_file)
G <- grm_out$G
ids <- rownames(G)

cat("Reading SNP matrix X from VCF...\n")
X <- build_X_from_vcf(vcf_file)

common_ids <- intersect(ids, rownames(X))
G <- G[common_ids, common_ids, drop = FALSE]
X <- X[common_ids, , drop = FALSE]
ids <- common_ids

cat("Filtering SNPs...\n")
keep_no_na <- colSums(is.na(X)) == 0
X <- X[, keep_no_na, drop = FALSE]

keep_poly <- apply(X, 2, function(z) var(z, na.rm = TRUE) > 0)
X <- X[, keep_poly, drop = FALSE]

cat("Individuals:", nrow(X), "\n")
cat("Markers:", ncol(X), "\n")
cat("Cores:", ncores, "\n")
cat("tolParInv:", tolParInv, "\n")
cat("nugget:", nugget, "\n")

cat("Running parallel simulation + GS...\n")
result <- run_parallel_simulation(
  X = X,
  G = G,
  ids = ids,
  corr_scenarios = corr_scenarios,
  nrep = nrep,
  h2 = h2_vec,
  trait_names = trait_names,
  causal_prop = causal_prop,
  k = kfold,
  reps = cv_reps,
  nugget = nugget,
  tolParInv = tolParInv,
  ncores = ncores,
  outdir = outdir
)

cat("Saving combined outputs...\n")
data.table::fwrite(result$accuracy,
                   file = file.path(outdir, "GS_accuracy_all_replicates.txt"),
                   sep = "\t")

data.table::fwrite(result$realized,
                   file = file.path(outdir, "Simulation_realized_parameters.txt"),
                   sep = "\t")

if (nrow(result$errors) > 0) {
  data.table::fwrite(result$errors,
                     file = file.path(outdir, "Errors_all_scenarios.txt"),
                     sep = "\t")
}

summary_acc <- summarize_accuracy(data.table::as.data.table(result$accuracy))
data.table::fwrite(summary_acc,
                   file = file.path(outdir, "GS_accuracy_summary.txt"),
                   sep = "\t")

acc_plot <- data.table::as.data.table(result$accuracy)
acc_plot <- acc_plot[is.finite(accuracy) & !is.na(accuracy)]
if (nrow(acc_plot) > 0) {
  plot_accuracy_boxplots(acc_plot, outdir = outdir)
}

cat("Done.\n")
save.image('run_simulation_gs_robust_MTfixed_v2.RData')

