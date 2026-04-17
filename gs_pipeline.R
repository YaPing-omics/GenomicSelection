# ============================================================
# gs_pipeline_empirical_MTfixed_olive_v5.R
# Empirical GS pipeline with:
#   - partial-masking MT CV
#   - proper trait-level covariance structure in sommer
#   - tolerant MT GEBV extraction across CV/full-fit structures
#   - automatic ID checking/filtering between pheno and G
#   - olive plot colors
# ============================================================

gs_load_packages <- function() {
  pkgs_cran <- c("data.table", "sommer", "ggplot2")
  pkgs_bioc <- c("SNPRelate", "gdsfmt")
  for (p in pkgs_cran) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  for (p in pkgs_bioc) {
    if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
  }
  invisible(TRUE)
}

read_pheno_txt <- function(path, header = TRUE) {
  gs_load_packages()
  data.table::fread(path, header = header, data.table = FALSE)
}

build_grm_from_vcf <- function(vcf_fn,
                               gds_fn = sub("\\.vcf(\\.gz)?$", ".gds", vcf_fn, ignore.case = TRUE),
                               method_vcf = "biallelic.only",
                               grm_method = "GCTA",
                               verbose = TRUE) {
  gs_load_packages()
  library(SNPRelate)
  if (!file.exists(vcf_fn)) stop("VCF not found: ", vcf_fn)
  if (file.exists(gds_fn)) file.remove(gds_fn)
  snpgdsVCF2GDS(vcf_fn, gds_fn, method = method_vcf, verbose = verbose)
  genofile <- snpgdsOpen(gds_fn)
  grm <- snpgdsGRM(genofile, method = grm_method, verbose = verbose)
  G <- grm$grm
  ids <- grm$sample.id
  rownames(G) <- colnames(G) <- ids
  snpgdsClose(genofile)
  list(G = G, sample_id = ids, gds_fn = gds_fn)
}

align_pheno_G <- function(pheno, G, id_col = "id", min_overlap = 10, verbose = TRUE) {
  if (!id_col %in% names(pheno)) stop("ID column not found in pheno: ", id_col)
  if (is.null(rownames(G)) || is.null(colnames(G))) stop("G must have rownames and colnames.")

  pheno[[id_col]] <- trimws(as.character(pheno[[id_col]]))
  rownames(G) <- trimws(as.character(rownames(G)))
  colnames(G) <- trimws(as.character(colnames(G)))

  if (!identical(rownames(G), colnames(G))) {
    stop("rownames(G) and colnames(G) must match and be in the same order.")
  }

  pheno_ids <- unique(pheno[[id_col]])
  g_ids <- rownames(G)

  missing_in_G <- setdiff(pheno_ids, g_ids)
  extra_in_G <- setdiff(g_ids, pheno_ids)
  common_ids <- intersect(pheno_ids, g_ids)

  if (verbose) {
    message("Phenotype IDs: ", length(pheno_ids))
    message("G IDs: ", length(g_ids))
    message("Common IDs: ", length(common_ids))
    if (length(missing_in_G) > 0) message("Dropping ", length(missing_in_G), " phenotype IDs not found in G.")
    if (length(extra_in_G) > 0) message(length(extra_in_G), " G IDs have no phenotype records.")
  }

  if (length(common_ids) < min_overlap) {
    stop("Too few overlapping IDs between pheno and G: ", length(common_ids))
  }

  pheno2 <- pheno[pheno[[id_col]] %in% common_ids, , drop = FALSE]
  pheno2[[id_col]] <- trimws(as.character(pheno2[[id_col]]))
  G2 <- G[common_ids, common_ids, drop = FALSE]

  list(pheno = pheno2, G = G2, common_ids = common_ids,
       missing_in_G = missing_in_G, extra_in_G = extra_in_G)
}

olive_fill <- c("Single-trait" = "darkolivegreen4", "Multi-trait" = "darkolivegreen2")
olive_color <- c("Single-trait" = "darkolivegreen4", "Multi-trait" = "darkolivegreen2")

save_cv_summary_plots <- function(cv_single = NULL, cv_multi = NULL, outdir = "GS_results", prefix = "CV") {
  gs_load_packages()
  library(ggplot2)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  parts <- list()
  if (!is.null(cv_single) && nrow(cv_single) > 0) {
    parts[["single"]] <- data.frame(trait = cv_single$trait, mean = cv_single$mean, sd = cv_single$sd,
                                    method = "Single-trait", stringsAsFactors = FALSE)
  }
  if (!is.null(cv_multi) && nrow(cv_multi) > 0) {
    parts[["multi"]] <- data.frame(trait = cv_multi$trait, mean = cv_multi$mean, sd = cv_multi$sd,
                                   method = "Multi-trait", stringsAsFactors = FALSE)
  }
  if (length(parts) == 0) return(invisible(NULL))
  plot_df <- do.call(rbind, parts); rownames(plot_df) <- NULL
  p1 <- ggplot(plot_df, aes(x = trait, y = mean, fill = method)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                  position = position_dodge(width = 0.8), width = 0.2) +
    scale_fill_manual(values = olive_fill) +
    theme_bw(base_size = 12) +
    labs(title = "Cross-validation accuracy by trait", x = "Trait", y = "Mean CV accuracy")
  ggsave(file.path(outdir, paste0(prefix, "_summary_barplot.png")), p1, width = 8, height = 5, dpi = 300)
  p2 <- ggplot(plot_df, aes(x = method, y = mean, color = method)) +
    geom_point(size = 2.8) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.15) +
    facet_wrap(~ trait, scales = "free_y") +
    scale_color_manual(values = olive_color) +
    theme_bw(base_size = 12) +
    labs(title = "CV accuracy comparison", x = "Method", y = "Mean CV accuracy")
  ggsave(file.path(outdir, paste0(prefix, "_summary_comparison.png")), p2, width = 8, height = 5, dpi = 300)
  invisible(plot_df)
}

save_cv_fold_boxplot <- function(cv_per_fold_single = NULL, cv_per_fold_multi = NULL,
                                 traits = NULL, outdir = "GS_results", prefix = "CV_fold") {
  gs_load_packages()
  library(ggplot2)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  pieces <- list()
  if (!is.null(cv_per_fold_multi) && nrow(as.data.frame(cv_per_fold_multi)) > 0) {
    dfm <- as.data.frame(cv_per_fold_multi)
    if (!is.null(traits)) {
      keep <- intersect(traits, colnames(dfm))
      dfm <- dfm[, keep, drop = FALSE]
    }
    if (ncol(dfm) > 0) {
      dfm$FoldRep <- seq_len(nrow(dfm))
      dtm <- data.table::as.data.table(dfm)
      longm <- data.table::melt(dtm, id.vars = "FoldRep", variable.name = "trait", value.name = "accuracy")
      longm$method <- "Multi-trait"
      pieces[["multi"]] <- longm
    }
  }
  if (!is.null(cv_per_fold_single) && length(cv_per_fold_single) > 0) {
    longs <- list()
    for (tr in names(cv_per_fold_single)) {
      dfx <- as.data.frame(cv_per_fold_single[[tr]])
      if (nrow(dfx) > 0) {
        acc <- as.numeric(dfx[, 1])
        longs[[tr]] <- data.frame(FoldRep = seq_along(acc), trait = tr, accuracy = acc,
                                  method = "Single-trait", stringsAsFactors = FALSE)
      }
    }
    if (length(longs) > 0) pieces[["single"]] <- do.call(rbind, longs)
  }
  if (length(pieces) == 0) return(invisible(NULL))
  plot_df <- data.table::rbindlist(pieces, fill = TRUE)
  plot_df <- plot_df[is.finite(accuracy) & !is.na(accuracy), ]
  if (nrow(plot_df) == 0) return(invisible(NULL))
  p <- ggplot(plot_df, aes(x = method, y = accuracy, fill = method)) +
    geom_boxplot(outlier.shape = 16, alpha = 0.8) +
    facet_wrap(~ trait, scales = "free_y") +
    scale_fill_manual(values = olive_fill) +
    theme_bw(base_size = 12) +
    labs(title = "Fold-level CV accuracy comparison", x = "Method", y = "CV accuracy")
  ggsave(file.path(outdir, paste0(prefix, "_boxplot.png")), p, width = 9, height = 5.5, dpi = 300)
  invisible(plot_df)
}

run_gs_sommer <- function(pheno, G, traits, id_col = "id", scale_traits = TRUE,
                          nugget = 1e-6, tolParInv = 1e-2, k = 5, reps = 3, seed = 123,
                          weights = NULL, predict_new = TRUE) {
  gs_load_packages()
  library(sommer)

  stopifnot(id_col %in% names(pheno))
  stopifnot(all(traits %in% names(pheno)))
  stopifnot(!is.null(rownames(G)), !is.null(colnames(G)))

  aligned <- align_pheno_G(pheno = pheno, G = G, id_col = id_col, min_overlap = 10, verbose = FALSE)
  pheno <- aligned$pheno
  G <- aligned$G

  train_ids <- unique(pheno[[id_col]])
  geno_ids  <- rownames(G)
  common_train <- intersect(train_ids, geno_ids)
  new_ids <- setdiff(geno_ids, common_train)

  Gs <- G + diag(nugget, nrow(G))
  rownames(Gs) <- colnames(Gs) <- geno_ids

  extract_gebv <- function(fit, traits, fallback_ids) {
    Uterm <- fit$U$`u:id`
    fallback_ids <- as.character(fallback_ids)

    if (is.null(Uterm)) stop("fit$U$`u:id` is NULL.")

    nm <- names(Uterm)

    # 1) exact trait names
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

    # 2) loose name match after stripping punctuation
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

    # 3) unnamed list with same length as traits; assume order matches traits
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

    # 4) matrix with one column per trait
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

  if (predict_new) {
    dat <- data.frame(id = geno_ids, stringsAsFactors = FALSE)
    for (tr in traits) dat[[tr]] <- NA_real_
    ph_sub <- pheno[match(common_train, pheno[[id_col]]), , drop = FALSE]
    dat[match(common_train, dat$id), traits] <- ph_sub[match(common_train, ph_sub[[id_col]]), traits]
  } else {
    dat <- data.frame(id = common_train, stringsAsFactors = FALSE)
    ph_sub <- pheno[match(common_train, pheno[[id_col]]), , drop = FALSE]
    for (tr in traits) dat[[tr]] <- ph_sub[[tr]]
  }

  if (scale_traits) {
    for (tr in traits) {
      m <- mean(dat[[tr]], na.rm = TRUE)
      s <- sd(dat[[tr]], na.rm = TRUE)
      if (is.na(s) || s == 0) s <- 1
      dat[[tr]] <- (dat[[tr]] - m) / s
    }
  }

  if (length(traits) == 1) {
    dat_fit <- data.frame(id = dat$id)
    dat_fit$y <- dat[[traits[1]]]
    fit_full <- mmer(y ~ 1, random = ~ vsr(id, Gu = Gs), rcov = ~ units,
                     data = dat_fit, tolParInv = tolParInv)
  } else {
    form <- as.formula(paste0("cbind(", paste(traits, collapse = ","), ") ~ 1"))
    fit_full <- mmer(form,
                     random = ~ vsr(id, Gu = Gs, Gtc = unsm(length(traits))),
                     rcov   = ~ vsr(units, Gtc = unsm(length(traits))),
                     data = dat, tolParInv = tolParInv)
  }

  gebv_full <- extract_gebv(fit_full, traits, fallback_ids = dat$id)

  index_tbl <- NULL
  if (!is.null(weights)) {
    stopifnot(all(traits %in% names(weights)))
    w_use <- weights[traits]
    idx <- as.numeric(gebv_full %*% w_use[colnames(gebv_full)])
    id_vec <- rownames(gebv_full)
    if (is.null(id_vec) || length(id_vec) == 0) id_vec <- dat$id
    index_tbl <- data.frame(id = id_vec, gebv_full, Index = idx, row.names = NULL)
    index_tbl <- index_tbl[order(-index_tbl$Index), ]
  }

  all_acc <- NULL
  cv_summary <- NULL

  if (k > 0 && reps > 0) {
    set.seed(seed)
    dat_tr <- data.frame(id = common_train, stringsAsFactors = FALSE)
    ph_tr <- pheno[match(common_train, pheno[[id_col]]), , drop = FALSE]
    for (tr in traits) dat_tr[[tr]] <- ph_tr[[tr]]
    dat_tr$id <- trimws(as.character(dat_tr$id))
    ids <- dat_tr$id

    if (length(traits) == 1) {
      acc_per_rep <- vector("list", reps)
      for (r in seq_len(reps)) {
        fold_id <- sample(rep(1:k, length.out = length(ids))); names(fold_id) <- ids
        acc_mat <- matrix(NA_real_, nrow = k, ncol = 1); colnames(acc_mat) <- traits
        for (f in 1:k) {
          test_ids <- names(fold_id)[fold_id == f]
          train_ids2 <- setdiff(ids, test_ids)
          dat_cv <- dat_tr
          dat_cv[dat_cv$id %in% test_ids, traits] <- NA
          mu  <- sapply(traits, function(tr) mean(dat_tr[dat_tr$id %in% train_ids2, tr], na.rm = TRUE))
          sdv <- sapply(traits, function(tr) sd(dat_tr[dat_tr$id %in% train_ids2, tr], na.rm = TRUE))
          sdv[is.na(sdv) | sdv == 0] <- 1
          dat_cv_scaled <- dat_cv
          for (tr in traits) dat_cv_scaled[[tr]] <- (dat_cv_scaled[[tr]] - mu[tr]) / sdv[tr]
          dat_fit <- data.frame(id = dat_cv_scaled$id); dat_fit$y <- dat_cv_scaled[[traits[1]]]
          fit_cv <- mmer(y ~ 1, random = ~ vsr(id, Gu = Gs[ids, ids]), rcov = ~ units,
                         data = dat_fit, tolParInv = tolParInv)
          gebv_cv <- extract_gebv(fit_cv, traits, fallback_ids = dat_cv_scaled$id)
          pred <- gebv_cv[match(test_ids, rownames(gebv_cv)), , drop = FALSE]
          rownames(pred) <- test_ids
          obs <- dat_tr[match(test_ids, dat_tr$id), traits, drop = FALSE]
          obs[[traits[1]]] <- (obs[[traits[1]]] - mu[traits[1]]) / sdv[traits[1]]
          acc_mat[f, traits[1]] <- cor(pred[, traits[1]], obs[[traits[1]]], use = "complete.obs")
        }
        acc_per_rep[[r]] <- acc_mat
      }
      all_acc <- do.call(rbind, acc_per_rep)
    } else {
      acc_target_list <- vector("list", length(traits)); names(acc_target_list) <- traits
      for (target_trait in traits) {
        acc_per_rep <- vector("list", reps)
        for (r in seq_len(reps)) {
          fold_id <- sample(rep(1:k, length.out = length(ids))); names(fold_id) <- ids
          acc_mat <- matrix(NA_real_, nrow = k, ncol = 1); colnames(acc_mat) <- target_trait
          for (f in 1:k) {
            test_ids <- names(fold_id)[fold_id == f]
            train_ids2 <- setdiff(ids, test_ids)
            dat_cv <- dat_tr
            dat_cv[dat_cv$id %in% test_ids, target_trait] <- NA
            mu  <- sapply(traits, function(tr) mean(dat_tr[dat_tr$id %in% train_ids2, tr], na.rm = TRUE))
            sdv <- sapply(traits, function(tr) sd(dat_tr[dat_tr$id %in% train_ids2, tr], na.rm = TRUE))
            sdv[is.na(sdv) | sdv == 0] <- 1
            dat_cv_scaled <- dat_cv
            for (tr in traits) dat_cv_scaled[[tr]] <- (dat_cv_scaled[[tr]] - mu[tr]) / sdv[tr]
            form <- as.formula(paste0("cbind(", paste(traits, collapse = ","), ") ~ 1"))
            fit_cv <- mmer(form,
                           random = ~ vsr(id, Gu = Gs[ids, ids], Gtc = unsm(length(traits))),
                           rcov   = ~ vsr(units, Gtc = unsm(length(traits))),
                           data = dat_cv_scaled, tolParInv = tolParInv)
            gebv_cv <- extract_gebv(fit_cv, traits, fallback_ids = dat_cv_scaled$id)
            pred <- gebv_cv[match(test_ids, rownames(gebv_cv)), , drop = FALSE]
            rownames(pred) <- test_ids
            obs <- dat_tr[match(test_ids, dat_tr$id), target_trait, drop = TRUE]
            obs <- (obs - mu[target_trait]) / sdv[target_trait]
            pred_vec <- pred[, target_trait]
            acc_mat[f, target_trait] <- cor(pred_vec, obs, use = "complete.obs")
          }
          acc_per_rep[[r]] <- acc_mat
        }
        acc_target_list[[target_trait]] <- do.call(rbind, acc_per_rep)
      }
      all_acc <- do.call(cbind, acc_target_list)
      all_acc <- as.matrix(all_acc)
      colnames(all_acc) <- traits
    }

    cv_summary <- data.frame(
      trait = traits,
      mean  = as.numeric(colMeans(all_acc, na.rm = TRUE)),
      sd    = as.numeric(apply(all_acc, 2, sd, na.rm = TRUE))
    )
  }

  list(train_ids = common_train, geno_ids = geno_ids, new_ids = new_ids,
       GEBV = gebv_full, selection_index = index_tbl,
       cv_per_fold = all_acc, cv_summary = cv_summary, fit_full = fit_full,
       id_alignment = aligned)
}

run_single_and_multi_GS <- function(pheno, G_all, traits, id_col = "id",
                                    outdir = "GS_results", k = 5, reps = 3, seed = 123) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  aligned <- align_pheno_G(pheno = pheno, G = G_all, id_col = id_col, min_overlap = 10, verbose = TRUE)
  pheno_use <- aligned$pheno
  G_use <- aligned$G

  if (length(aligned$missing_in_G) > 0) {
    data.table::fwrite(data.frame(id = aligned$missing_in_G),
                       file.path(outdir, "IDs_missing_in_G.txt"), sep = "\t")
  }
  if (length(aligned$extra_in_G) > 0) {
    data.table::fwrite(data.frame(id = aligned$extra_in_G),
                       file.path(outdir, "IDs_in_G_without_pheno.txt"), sep = "\t")
  }

  cat("Running multi-trait GS...\n")
  res_mt <- run_gs_sommer(pheno = pheno_use, G = G_use, traits = traits, id_col = id_col,
                          k = k, reps = reps, seed = seed, predict_new = TRUE)
  if (!is.null(res_mt$cv_summary)) {
    data.table::fwrite(res_mt$cv_summary, file.path(outdir, "CV_multitrait.txt"), sep="\t")
  }
  gebv_mt <- res_mt$GEBV
  pred_mt <- data.frame(id = rownames(gebv_mt), gebv_mt)
  data.table::fwrite(pred_mt, file.path(outdir, "Pred_multitrait_all.txt"), sep="\t")
  if (length(res_mt$new_ids) > 0) {
    pred_mt_new <- pred_mt[pred_mt$id %in% res_mt$new_ids, ]
    data.table::fwrite(pred_mt_new, file.path(outdir, "Pred_multitrait_new.txt"), sep="\t")
  }

  cat("Running single-trait GS...\n")
  cv_single_list <- list(); pred_single_list <- list(); cv_fold_single_list <- list()
  for (tr in traits) {
    cat("  Trait:", tr, "\n")
    res_st <- run_gs_sommer(pheno = pheno_use, G = G_use, traits = c(tr), id_col = id_col,
                            k = k, reps = reps, seed = seed, predict_new = TRUE)
    if (!is.null(res_st$cv_summary)) {
      cv_single_list[[tr]] <- data.frame(trait = tr, mean = res_st$cv_summary$mean, sd = res_st$cv_summary$sd)
    }
    if (!is.null(res_st$cv_per_fold)) cv_fold_single_list[[tr]] <- as.data.frame(res_st$cv_per_fold)
    gebv_st <- res_st$GEBV
    pred_single_list[[tr]] <- data.frame(id = rownames(gebv_st), trait = tr, GEBV = gebv_st[,1])
  }
  cv_single <- NULL
  if (length(cv_single_list) > 0) {
    cv_single <- do.call(rbind, cv_single_list)
    data.table::fwrite(cv_single, file.path(outdir, "CV_singletrait.txt"), sep="\t")
  }
  pred_single_long <- do.call(rbind, pred_single_list)
  pred_single_long <- data.table::as.data.table(pred_single_long)
  pred_single_wide <- data.table::dcast(pred_single_long, id ~ trait, value.var="GEBV")
  data.table::fwrite(pred_single_wide, file.path(outdir, "Pred_singletrait_all.txt"), sep="\t")

  save_cv_summary_plots(cv_single = cv_single, cv_multi = res_mt$cv_summary, outdir = outdir, prefix = "CV")
  save_cv_fold_boxplot(cv_per_fold_single = cv_fold_single_list, cv_per_fold_multi = res_mt$cv_per_fold,
                       traits = traits, outdir = outdir, prefix = "CV_fold")
  list(multitrait = res_mt, singletrait = pred_single_wide, cv_single = cv_single, id_alignment = aligned)
}

gs_example_run_and_save <- function(pheno_file, vcf_file, traits, id_col = "id",
                                    outdir = "GS_results", k = 5, reps = 3, seed = 123,
                                    tolParInv = 1e-2, nugget = 1e-6,
                                    scale_traits = TRUE, predict_new = TRUE) {
  gs_load_packages(); library(data.table)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  pheno <- read_pheno_txt(pheno_file)
  grm_out <- build_grm_from_vcf(vcf_file)
  G_all <- grm_out$G

  aligned <- align_pheno_G(pheno = pheno, G = G_all, id_col = id_col, min_overlap = 10, verbose = TRUE)
  pheno_use <- aligned$pheno
  G_use <- aligned$G

  res_cv <- run_gs_sommer(pheno = pheno_use, G = G_use, traits = traits, id_col = id_col,
                          k = k, reps = reps, seed = seed, tolParInv = tolParInv,
                          nugget = nugget, scale_traits = scale_traits, predict_new = FALSE)
  if (!is.null(res_cv$cv_summary)) {
    data.table::fwrite(res_cv$cv_summary, file = file.path(outdir, "CV_summary.txt"), sep = "\t")
  }
  if (!is.null(res_cv$cv_per_fold)) {
    cv_fold_df <- as.data.frame(res_cv$cv_per_fold); cv_fold_df$FoldRep <- seq_len(nrow(cv_fold_df))
    data.table::fwrite(cv_fold_df, file = file.path(outdir, "CV_per_fold.txt"), sep = "\t")
  }
  save_cv_summary_plots(cv_single = NULL, cv_multi = res_cv$cv_summary, outdir = outdir, prefix = "Example_CV")
  save_cv_fold_boxplot(cv_per_fold_single = NULL, cv_per_fold_multi = res_cv$cv_per_fold,
                       traits = traits, outdir = outdir, prefix = "Example_CV_fold")
  res_pred <- run_gs_sommer(pheno = pheno_use, G = G_use, traits = traits, id_col = id_col,
                            k = 0, reps = 0, seed = seed, tolParInv = tolParInv,
                            nugget = nugget, scale_traits = scale_traits, predict_new = predict_new)
  gebv_all <- res_pred$GEBV
  pred_all <- data.frame(id = rownames(gebv_all), gebv_all, row.names = NULL)
  fwrite(pred_all, file = file.path(outdir, "Pred_GEBV_all_genotypes.txt"), sep = "\t")
  new_ids <- res_pred$new_ids
  if (length(new_ids) > 0) {
    gebv_new <- gebv_all[new_ids, , drop = FALSE]
    pred_new <- data.frame(id = rownames(gebv_new), gebv_new, row.names = NULL)
    fwrite(pred_new, file = file.path(outdir, "Pred_GEBV_new_genotypes.txt"), sep = "\t")
  } else {
    message("No new genotypes detected (new_ids length = 0).")
  }
  list(pheno = pheno_use, G_all = G_use, cv = res_cv, pred = res_pred, id_alignment = aligned)
}
