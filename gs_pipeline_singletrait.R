# ============================================================
# gs_pipeline_singletrait_olive_v2.R
# Single-trait GS pipeline with:
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

flatten_cv_per_fold <- function(cv_per_fold, trait, marker_n = NA_integer_, marker_rep = NA_integer_) {
  if (is.null(cv_per_fold)) return(data.table::data.table())
  x <- data.table::as.data.table(as.data.frame(cv_per_fold))
  if (nrow(x) == 0) return(data.table::data.table())
  data.table::data.table(
    marker_n = marker_n,
    marker_rep = marker_rep,
    FoldRep = seq_len(nrow(x)),
    trait = trait,
    accuracy = as.numeric(x[[1]]),
    method = "Single-trait"
  )
}

run_gs_singletrait <- function(pheno, G, trait, id_col = "id",
                               scale_trait = TRUE, nugget = 1e-6, tolParInv = 1e-2,
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

  extract_gebv <- function(fit, fallback_ids) {
    Uterm <- fit$U$`u:id`
    u <- Uterm[[1]]
    vals <- if (is.matrix(u)) u[, 1] else as.numeric(u)
    ids <- rownames(u)
    if (is.null(ids) || length(ids) == 0) ids <- names(u)
    if (is.null(ids) || length(ids) == 0) ids <- fallback_ids
    out <- vals[match(fallback_ids, ids)]
    out <- matrix(out, ncol = 1)
    rownames(out) <- fallback_ids
    colnames(out) <- trait
    out
  }

  if (predict_new) {
    dat <- data.frame(id = geno_ids, stringsAsFactors = FALSE)
    dat[[trait]] <- NA_real_
    ph_sub <- pheno[match(common_train, pheno[[id_col]]), , drop = FALSE]
    dat[match(common_train, dat$id), trait] <- ph_sub[[trait]]
  } else {
    dat <- data.frame(id = common_train, stringsAsFactors = FALSE)
    ph_sub <- pheno[match(common_train, pheno[[id_col]]), , drop = FALSE]
    dat[[trait]] <- ph_sub[[trait]]
  }

  if (scale_trait) {
    m <- mean(dat[[trait]], na.rm = TRUE)
    s <- sd(dat[[trait]], na.rm = TRUE)
    if (is.na(s) || s == 0) s <- 1
    dat[[trait]] <- (dat[[trait]] - m) / s
  }

  dat_fit <- data.frame(id = dat$id, y = dat[[trait]])
  fit_full <- sommer::mmer(y ~ 1, random = ~ sommer::vsr(id, Gu = Gs), rcov = ~ units,
                           data = dat_fit, tolParInv = tolParInv)

  gebv_full <- extract_gebv(fit_full, dat$id)

  all_acc <- NULL
  cv_summary <- NULL
  if (k > 0 && reps > 0) {
    set.seed(seed)
    dat_tr <- data.frame(id = common_train, stringsAsFactors = FALSE)
    ph_tr <- pheno[match(common_train, pheno[[id_col]]), , drop = FALSE]
    dat_tr[[trait]] <- ph_tr[[trait]]
    ids <- as.character(dat_tr$id)

    acc_per_rep <- vector("list", reps)
    for (r in seq_len(reps)) {
      fold_id <- sample(rep(1:k, length.out = length(ids)))
      names(fold_id) <- ids
      acc_mat <- matrix(NA_real_, nrow = k, ncol = 1)
      colnames(acc_mat) <- trait

      for (f in 1:k) {
        test_ids <- names(fold_id)[fold_id == f]
        train_ids2 <- setdiff(ids, test_ids)

        dat_cv <- dat_tr
        dat_cv[dat_cv$id %in% test_ids, trait] <- NA

        mu <- mean(dat_tr[dat_tr$id %in% train_ids2, trait], na.rm = TRUE)
        sdv <- sd(dat_tr[dat_tr$id %in% train_ids2, trait], na.rm = TRUE)
        if (is.na(sdv) || sdv == 0) sdv <- 1

        dat_fit <- data.frame(id = dat_cv$id, y = (dat_cv[[trait]] - mu) / sdv)
        fit_cv <- sommer::mmer(y ~ 1, random = ~ sommer::vsr(id, Gu = Gs), rcov = ~ units,
                               data = dat_fit, tolParInv = tolParInv)
        gebv_cv <- extract_gebv(fit_cv, dat_fit$id)
        pred <- gebv_cv[match(test_ids, rownames(gebv_cv)), , drop = FALSE]
        rownames(pred) <- test_ids

        obs <- dat_tr[match(test_ids, dat_tr$id), trait, drop = TRUE]
        obs <- (obs - mu) / sdv
        acc_mat[f, trait] <- cor(pred[, trait], obs, use = "complete.obs")
      }
      acc_per_rep[[r]] <- acc_mat
    }

    all_acc <- do.call(rbind, acc_per_rep)
    cv_summary <- data.frame(trait = trait,
                             mean = mean(all_acc[, trait], na.rm = TRUE),
                             sd = sd(all_acc[, trait], na.rm = TRUE))
  }

  list(GEBV = gebv_full, cv_per_fold = all_acc, cv_summary = cv_summary, new_ids = new_ids)
}

run_singletrait_GS <- function(pheno, G_all, traits, id_col = "id",
                               outdir = "GS_results_singletrait", k = 5, reps = 3, seed = 123,
                               marker_n = NA_integer_, marker_rep = NA_integer_) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  aligned <- align_pheno_G(pheno, G_all, id_col = id_col, min_overlap = 10, verbose = TRUE)
  pheno_use <- aligned$pheno
  G_use <- aligned$G

  cv_summary_list <- list()
  cv_each_list <- list()
  pred_list <- list()

  for (tr in traits) {
    cat("Running single-trait GS for:", tr, "\\n")
    res <- run_gs_singletrait(pheno_use, G_use, tr, id_col = id_col, k = k, reps = reps, seed = seed)

    if (!is.null(res$cv_summary)) {
      tmp <- res$cv_summary
      tmp$marker_n <- marker_n
      tmp$marker_rep <- marker_rep
      cv_summary_list[[tr]] <- tmp
    }
    cv_each_list[[tr]] <- flatten_cv_per_fold(res$cv_per_fold, tr, marker_n = marker_n, marker_rep = marker_rep)
    pred_list[[tr]] <- data.frame(id = rownames(res$GEBV), trait = tr, GEBV = res$GEBV[, 1])
  }

  cv_summary <- data.table::rbindlist(cv_summary_list, fill = TRUE)
  if (nrow(cv_summary) > 0) data.table::fwrite(cv_summary, file.path(outdir, "CV_singletrait.txt"), sep = "\t")

  cv_each <- data.table::rbindlist(cv_each_list, fill = TRUE)
  if (nrow(cv_each) > 0) data.table::fwrite(cv_each, file.path(outdir, "CV_singletrait_each_time.txt"), sep = "\t")

  pred_long <- data.table::rbindlist(pred_list, fill = TRUE)
  pred_wide <- data.table::dcast(pred_long, id ~ trait, value.var = "GEBV")
  data.table::fwrite(pred_wide, file.path(outdir, "Pred_singletrait_all.txt"), sep = "\t")

  invisible(list(cv_summary = cv_summary, cv_each_time = cv_each, prediction = pred_wide))
}

run_singletrait_marker_number_experiment <- function(pheno_file, vcf_file, traits, marker_numbers,
                                                     marker_reps = 1,
                                                     id_col = "id", base_outdir = "GS_marker_number_results_singletrait",
                                                     k = 5, reps = 3, seed = 123,
                                                     marker_seed_base = 1000) {
  if (!dir.exists(base_outdir)) dir.create(base_outdir, recursive = TRUE)
  pheno <- read_pheno_txt(pheno_file)

  all_summary <- list()
  all_each <- list()

  for (mn in marker_numbers) {
    for (mr in seq_len(marker_reps)) {
      outdir <- file.path(base_outdir, paste0("markers_", mn, "_rep_", mr))
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

      marker_seed <- marker_seed_base + mn + mr * 10000

      cat("\\nRunning single-trait marker_n =", mn, " marker_rep =", mr, "\\n")
      grm_out <- build_grm_from_vcf(vcf_file, marker_n = mn, marker_seed = marker_seed)
      data.table::fwrite(data.frame(marker_n_requested = mn,
                                    marker_rep = mr,
                                    marker_n_used = grm_out$marker_n_used,
                                    total_markers = grm_out$total_markers,
                                    marker_seed = marker_seed),
                         file.path(outdir, "Marker_number_info.txt"), sep = "\t")

      out <- run_singletrait_GS(pheno, grm_out$G, traits, id_col = id_col, outdir = outdir,
                                k = k, reps = reps, seed = seed,
                                marker_n = mn, marker_rep = mr)

      if (nrow(out$cv_summary) > 0) {
        tmp <- data.table::as.data.table(out$cv_summary)
        tmp[, method := "Single-trait"]
        data.table::setcolorder(tmp, c("marker_n", "marker_rep", "method", "trait", "mean", "sd"))
        all_summary[[paste0(mn, "_", mr)]] <- tmp
      }

      if (nrow(out$cv_each_time) > 0) {
        tmp2 <- data.table::as.data.table(out$cv_each_time)
        data.table::setcolorder(tmp2, c("marker_n", "marker_rep", "method", "trait", "FoldRep", "accuracy"))
        all_each[[paste0(mn, "_", mr)]] <- tmp2
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
