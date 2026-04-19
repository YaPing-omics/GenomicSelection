# ============================================================
# gs_pipeline_multitrait_olive_v2.R
# Multi-trait GS pipeline with:
#   - target-trait-only masking in CV
#   - proper trait covariance structure in sommer
#   - marker-number control for GRM construction
#   - marker-subset replication for fair SNP-number comparison
#   - automatic ID checking/filtering between pheno and G
#   - CV summary + fold-level export for replotting
# ============================================================

gs_load_packages <- function() {
  pkgs_cran <- c("data.table", "sommer")
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
                               marker_n = NULL,
                               marker_seed = 123,
                               verbose = TRUE) {
  gs_load_packages()
  library(SNPRelate)
  library(gdsfmt)

  if (!file.exists(vcf_fn)) stop("VCF not found: ", vcf_fn)
  if (file.exists(gds_fn)) file.remove(gds_fn)
  snpgdsVCF2GDS(vcf_fn, gds_fn, method = method_vcf, verbose = verbose)

  genofile <- snpgdsOpen(gds_fn)
  on.exit(snpgdsClose(genofile), add = TRUE)

  sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
  snp.id <- read.gdsn(index.gdsn(genofile, "snp.id"))
  total_markers <- length(snp.id)

  snp.id.use <- snp.id
  if (!is.null(marker_n)) {
    marker_n <- as.integer(marker_n)
    if (length(marker_n) != 1 || is.na(marker_n) || marker_n <= 0) {
      stop("marker_n must be NULL or one positive integer.")
    }
    if (marker_n > total_markers) {
      warning("marker_n > total available markers; using all markers.")
      marker_n <- total_markers
    }
    set.seed(marker_seed)
    snp.id.use <- sort(sample(snp.id, size = marker_n, replace = FALSE))
  }

  grm <- snpgdsGRM(genofile, sample.id = sample.id, snp.id = snp.id.use,
                   method = grm_method, verbose = verbose)

  G <- grm$grm
  ids <- grm$sample.id
  rownames(G) <- colnames(G) <- ids

  list(G = G, sample_id = ids, marker_n_requested = marker_n,
       marker_n_used = length(snp.id.use), total_markers = total_markers,
       marker_seed = marker_seed)
}

align_pheno_G <- function(pheno, G, id_col = "id", min_overlap = 10, verbose = TRUE) {
  pheno[[id_col]] <- trimws(as.character(pheno[[id_col]]))
  rownames(G) <- trimws(as.character(rownames(G)))
  colnames(G) <- trimws(as.character(colnames(G)))

  common_ids <- intersect(unique(pheno[[id_col]]), rownames(G))
  if (verbose) message("Common IDs: ", length(common_ids))
  if (length(common_ids) < min_overlap) stop("Too few overlapping IDs between pheno and G.")

  pheno2 <- pheno[pheno[[id_col]] %in% common_ids, , drop = FALSE]
  G2 <- G[common_ids, common_ids, drop = FALSE]

  list(pheno = pheno2, G = G2, common_ids = common_ids)
}

flatten_cv_per_fold <- function(cv_per_fold, traits, marker_n = NA_integer_, marker_rep = NA_integer_) {
  if (is.null(cv_per_fold)) return(data.table::data.table())
  df <- as.data.frame(cv_per_fold)
  if (nrow(df) == 0 || ncol(df) == 0) return(data.table::data.table())
  keep <- intersect(traits, colnames(df))
  if (length(keep) == 0) return(data.table::data.table())
  df <- df[, keep, drop = FALSE]
  df$FoldRep <- seq_len(nrow(df))
  dt <- data.table::as.data.table(df)
  long <- data.table::melt(dt, id.vars = "FoldRep", variable.name = "trait", value.name = "accuracy")
  long[, method := "Multi-trait"]
  long[, marker_n := marker_n]
  long[, marker_rep := marker_rep]
  data.table::setcolorder(long, c("marker_n", "marker_rep", "method", "trait", "FoldRep", "accuracy"))
  long[]
}

run_gs_multitrait <- function(pheno, G, traits, id_col = "id",
                              scale_traits = TRUE, nugget = 1e-6, tolParInv = 1e-2,
                              k = 5, reps = 3, seed = 123, predict_new = TRUE) {
  gs_load_packages()
  library(sommer)

  aligned <- align_pheno_G(pheno, G, id_col = id_col, min_overlap = 10, verbose = FALSE)
  pheno <- aligned$pheno
  G <- aligned$G

  train_ids <- unique(pheno[[id_col]])
  geno_ids <- rownames(G)
  common_train <- intersect(train_ids, geno_ids)
  new_ids <- setdiff(geno_ids, common_train)

  Gs <- G + diag(nugget, nrow(G))
  rownames(Gs) <- colnames(Gs) <- geno_ids

  extract_gebv <- function(fit, traits, fallback_ids) {
    Uterm <- fit$U$`u:id`
    fallback_ids <- as.character(fallback_ids)
    nm <- names(Uterm)

    if (!is.null(nm) && all(traits %in% nm)) {
      out <- sapply(traits, function(tr) {
        u <- Uterm[[tr]]
        vals <- if (is.matrix(u)) u[, 1] else as.numeric(u)
        ids <- rownames(u)
        if (is.null(ids) || length(ids) == 0) ids <- names(u)
        if (is.null(ids) || length(ids) == 0) ids <- fallback_ids
        vals[match(fallback_ids, ids)]
      })
      out <- as.matrix(out)
      rownames(out) <- fallback_ids
      colnames(out) <- traits
      return(out)
    }

    if (is.list(Uterm) && length(Uterm) == length(traits)) {
      out <- sapply(seq_along(traits), function(i) {
        u <- Uterm[[i]]
        vals <- if (is.matrix(u)) u[, 1] else as.numeric(u)
        ids <- rownames(u)
        if (is.null(ids) || length(ids) == 0) ids <- names(u)
        if (is.null(ids) || length(ids) == 0) ids <- fallback_ids
        vals[match(fallback_ids, ids)]
      })
      out <- as.matrix(out)
      rownames(out) <- fallback_ids
      colnames(out) <- traits
      return(out)
    }

    stop("Unsupported structure in fit$U$`u:id`.")
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

  form <- stats::as.formula(paste0("cbind(", paste(traits, collapse = ","), ") ~ 1"))
  fit_full <- sommer::mmer(form,
                           random = ~ sommer::vsr(id, Gu = Gs, Gtc = sommer::unsm(length(traits))),
                           rcov   = ~ sommer::vsr(units, Gtc = sommer::unsm(length(traits))),
                           data = dat, tolParInv = tolParInv)
  gebv_full <- extract_gebv(fit_full, traits, dat$id)

  all_acc <- NULL
  cv_summary <- NULL

  if (k > 0 && reps > 0) {
    set.seed(seed)
    dat_tr <- data.frame(id = common_train, stringsAsFactors = FALSE)
    ph_tr <- pheno[match(common_train, pheno[[id_col]]), , drop = FALSE]
    for (tr in traits) dat_tr[[tr]] <- ph_tr[[tr]]
    ids <- as.character(dat_tr$id)

    acc_target_list <- vector("list", length(traits))
    names(acc_target_list) <- traits

    for (target_trait in traits) {
      acc_per_rep <- vector("list", reps)
      for (r in seq_len(reps)) {
        fold_id <- sample(rep(1:k, length.out = length(ids)))
        names(fold_id) <- ids
        acc_mat <- matrix(NA_real_, nrow = k, ncol = 1)
        colnames(acc_mat) <- target_trait

        for (f in 1:k) {
          test_ids <- names(fold_id)[fold_id == f]
          train_ids2 <- setdiff(ids, test_ids)

          dat_cv <- dat_tr
          dat_cv[dat_cv$id %in% test_ids, target_trait] <- NA

          mu <- sapply(traits, function(tr) mean(dat_tr[dat_tr$id %in% train_ids2, tr], na.rm = TRUE))
          sdv <- sapply(traits, function(tr) sd(dat_tr[dat_tr$id %in% train_ids2, tr], na.rm = TRUE))
          sdv[is.na(sdv) | sdv == 0] <- 1

          dat_cv_scaled <- dat_cv
          for (tr in traits) dat_cv_scaled[[tr]] <- (dat_cv_scaled[[tr]] - mu[tr]) / sdv[tr]

          fit_cv <- sommer::mmer(form,
                                 random = ~ sommer::vsr(id, Gu = Gs, Gtc = sommer::unsm(length(traits))),
                                 rcov   = ~ sommer::vsr(units, Gtc = sommer::unsm(length(traits))),
                                 data = dat_cv_scaled, tolParInv = tolParInv)

          gebv_cv <- extract_gebv(fit_cv, traits, dat_cv_scaled$id)
          pred <- gebv_cv[match(test_ids, rownames(gebv_cv)), , drop = FALSE]
          rownames(pred) <- test_ids

          obs <- dat_tr[match(test_ids, dat_tr$id), target_trait, drop = TRUE]
          obs <- (obs - mu[target_trait]) / sdv[target_trait]
          acc_mat[f, target_trait] <- cor(pred[, target_trait], obs, use = "complete.obs")
        }
        acc_per_rep[[r]] <- acc_mat
      }
      acc_target_list[[target_trait]] <- do.call(rbind, acc_per_rep)
    }

    all_acc <- do.call(cbind, acc_target_list)
    all_acc <- as.matrix(all_acc)
    colnames(all_acc) <- traits
    cv_summary <- data.frame(trait = traits,
                             mean = as.numeric(colMeans(all_acc, na.rm = TRUE)),
                             sd = as.numeric(apply(all_acc, 2, sd, na.rm = TRUE)))
  }

  list(GEBV = gebv_full, cv_per_fold = all_acc, cv_summary = cv_summary, new_ids = new_ids)
}

run_multitrait_GS <- function(pheno, G_all, traits, id_col = "id",
                              outdir = "GS_results_multitrait", k = 5, reps = 3, seed = 123,
                              marker_n = NA_integer_, marker_rep = NA_integer_) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  out <- run_gs_multitrait(pheno, G_all, traits, id_col = id_col, k = k, reps = reps, seed = seed)

  cv_summary <- data.table::data.table()
  if (!is.null(out$cv_summary) && nrow(out$cv_summary) > 0) {
    cv_summary <- data.table::as.data.table(out$cv_summary)
    cv_summary[, marker_n := marker_n]
    cv_summary[, marker_rep := marker_rep]
    cv_summary[, method := "Multi-trait"]
    data.table::setcolorder(cv_summary, c("marker_n", "marker_rep", "method", "trait", "mean", "sd"))
    data.table::fwrite(cv_summary, file.path(outdir, "CV_multitrait.txt"), sep = "\t")
  }

  cv_each <- flatten_cv_per_fold(out$cv_per_fold, traits, marker_n = marker_n, marker_rep = marker_rep)
  if (nrow(cv_each) > 0) {
    data.table::fwrite(cv_each, file.path(outdir, "CV_multitrait_each_time.txt"), sep = "\t")
  }

  pred <- data.frame(id = rownames(out$GEBV), out$GEBV)
  data.table::fwrite(pred, file.path(outdir, "Pred_multitrait_all.txt"), sep = "\t")

  invisible(list(cv_summary = cv_summary, cv_each_time = cv_each, prediction = pred))
}

run_multitrait_marker_number_experiment <- function(pheno_file, vcf_file, traits, marker_numbers,
                                                    marker_reps = 1,
                                                    id_col = "id", base_outdir = "GS_marker_number_results_multitrait",
                                                    k = 5, reps = 3, seed = 123, marker_seed_base = 1000) {
  if (!dir.exists(base_outdir)) dir.create(base_outdir, recursive = TRUE)
  pheno <- read_pheno_txt(pheno_file)

  all_summary <- list()
  all_each <- list()

  for (mn in marker_numbers) {
    for (mr in seq_len(marker_reps)) {
      outdir <- file.path(base_outdir, paste0("markers_", mn, "_rep_", mr))
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

      marker_seed <- marker_seed_base + mn + mr * 10000

      cat("\nRunning multi-trait marker_n =", mn, " marker_rep =", mr, "\n")
      grm_out <- build_grm_from_vcf(vcf_file, marker_n = mn, marker_seed = marker_seed)
      data.table::fwrite(data.frame(marker_n_requested = mn,
                                    marker_rep = mr,
                                    marker_n_used = grm_out$marker_n_used,
                                    total_markers = grm_out$total_markers,
                                    marker_seed = marker_seed),
                         file.path(outdir, "Marker_number_info.txt"), sep = "\t")

      out <- run_multitrait_GS(pheno, grm_out$G, traits, id_col = id_col, outdir = outdir,
                               k = k, reps = reps, seed = seed,
                               marker_n = mn, marker_rep = mr)

      if (!is.null(out$cv_summary) && nrow(out$cv_summary) > 0) {
        all_summary[[paste0(mn, "_", mr)]] <- out$cv_summary
      }

      if (!is.null(out$cv_each_time) && nrow(out$cv_each_time) > 0) {
        all_each[[paste0(mn, "_", mr)]] <- out$cv_each_time
      }
    }
  }

  if (length(all_summary) > 0) {
    data.table::fwrite(data.table::rbindlist(all_summary, fill = TRUE),
                       file.path(base_outdir, "CV_summary_across_marker_numbers.txt"), sep = "\t")
  }
  if (length(all_each) > 0) {
    data.table::fwrite(data.table::rbindlist(all_each, fill = TRUE),
                       file.path(base_outdir, "CV_each_time_across_marker_numbers.txt"), sep = "\t")
  }
}
