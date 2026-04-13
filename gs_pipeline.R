# ============================================================
# gs_pipeline.R
# Functions for: phenotype import, VCF->GDS->GRM, GS (sommer),
# single/multi-trait CV (optional) + prediction for new genotypes
# ============================================================

# ---- Packages ----
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

# ---- Read phenotype from txt/csv/tsv (auto) ----
read_pheno_txt <- function(path, header = TRUE) {
  gs_load_packages()
  data.table::fread(path, header = header, data.table = FALSE)
}

# ---- Build GRM from VCF ----
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

# ---- GS with sommer (single/multi-trait), CV optional, prediction optional ----
run_gs_sommer <- function(pheno, G, traits,
                          id_col = "id",
                          scale_traits = TRUE,
                          nugget = 1e-6,
                          tolParInv = 1e-2,
                          k = 5, reps = 3, seed = 123,
                          weights = NULL,
                          predict_new = TRUE) {
  
  gs_load_packages()
  library(sommer)
  
  stopifnot(id_col %in% names(pheno))
  stopifnot(all(traits %in% names(pheno)))
  stopifnot(!is.null(rownames(G)), !is.null(colnames(G)))
  
  train_ids <- unique(pheno[[id_col]])
  geno_ids  <- rownames(G)
  
  common_train <- intersect(train_ids, geno_ids)
  if (length(common_train) < 10) stop("Too few overlapping IDs between pheno and GRM (G).")
  
  new_ids <- setdiff(geno_ids, common_train)
  
  # stabilize GRM
  Gs <- G + diag(nugget, nrow(G))
  
  # helper: robust GEBV extraction
  extract_gebv <- function(fit, traits, fallback_ids) {
    Uterm <- fit$U$`u:id`
    
    if (!is.null(names(Uterm)) && all(traits %in% names(Uterm))) {
      out <- do.call(cbind, lapply(traits, function(tr) {
        u <- Uterm[[tr]]
        if (is.matrix(u)) u[, 1] else as.numeric(u)
      }))
      colnames(out) <- traits
      
      rid <- NULL
      u0 <- Uterm[[traits[1]]]
      if (is.matrix(u0) && !is.null(rownames(u0))) rid <- rownames(u0)
      if (is.null(rid) || length(rid) == 0) rid <- fallback_ids
      rownames(out) <- rid
      return(out)
    }
    
    if (!is.null(names(Uterm)) && "y" %in% names(Uterm)) {
      u <- Uterm$y
      v <- if (is.matrix(u)) u[, 1] else as.numeric(u)
      out <- matrix(v, ncol = 1)
      colnames(out) <- traits[1]
      
      rid <- NULL
      if (is.matrix(u) && !is.null(rownames(u))) rid <- rownames(u)
      if ((is.null(rid) || length(rid) == 0) && !is.null(names(u))) rid <- names(u)
      if (is.null(rid) || length(rid) == 0) rid <- fallback_ids
      rownames(out) <- rid
      return(out)
    }
    
    u1 <- Uterm[[1]]
    v <- if (is.matrix(u1)) u1[, 1] else as.numeric(u1)
    out <- matrix(v, ncol = 1)
    colnames(out) <- traits[1]
    
    rid <- NULL
    if (is.matrix(u1) && !is.null(rownames(u1))) rid <- rownames(u1)
    if ((is.null(rid) || length(rid) == 0) && !is.null(names(u1))) rid <- names(u1)
    if (is.null(rid) || length(rid) == 0) rid <- fallback_ids
    rownames(out) <- rid
    return(out)
  }
  
  # build data for fitting
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
  
  # scale using available values only
  if (scale_traits) {
    for (tr in traits) {
      m <- mean(dat[[tr]], na.rm = TRUE)
      s <- sd(dat[[tr]], na.rm = TRUE)
      if (is.na(s) || s == 0) s <- 1
      dat[[tr]] <- (dat[[tr]] - m) / s
    }
  }
  
  # fit full model
  if (length(traits) == 1) {
    tr <- traits[1]
    dat_fit <- data.frame(id = dat$id)
    dat_fit$y <- dat[[tr]]
    
    fit_full <- mmer(
      y ~ 1,
      random = ~ vsr(id, Gu = Gs),
      rcov   = ~ units,
      data   = dat_fit,
      tolParInv = tolParInv
    )
  } else {
    form <- as.formula(paste0("cbind(", paste(traits, collapse = ","), ") ~ 1"))
    fit_full <- mmer(
      form,
      random = ~ vsr(id, Gu = Gs),
      rcov   = ~ units,
      data   = dat,
      tolParInv = tolParInv
    )
  }
  
  gebv_full <- extract_gebv(fit_full, traits, fallback_ids = dat$id)
  
  # selection index
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
  
  # CV optional (training only)
  all_acc <- NULL
  cv_summary <- NULL
  
  if (k > 0 && reps > 0) {
    set.seed(seed)
    
    dat_tr <- data.frame(id = common_train, stringsAsFactors = FALSE)
    ph_tr <- pheno[match(common_train, pheno[[id_col]]), , drop = FALSE]
    for (tr in traits) dat_tr[[tr]] <- ph_tr[[tr]]
    
    ids <- dat_tr$id
    acc_per_rep <- vector("list", reps)
    
    for (r in seq_len(reps)) {
      fold_id <- sample(rep(1:k, length.out = length(ids)))
      names(fold_id) <- ids
      
      acc_mat <- matrix(NA_real_, nrow = k, ncol = length(traits))
      colnames(acc_mat) <- traits
      
      for (f in 1:k) {
        test_ids  <- names(fold_id)[fold_id == f]
        train_ids2 <- setdiff(ids, test_ids)
        
        dat_cv <- dat_tr
        dat_cv[dat_cv$id %in% test_ids, traits] <- NA
        
        mu  <- sapply(traits, function(tr) mean(dat_tr[dat_tr$id %in% train_ids2, tr], na.rm = TRUE))
        sdv <- sapply(traits, function(tr) sd(dat_tr[dat_tr$id %in% train_ids2, tr],  na.rm = TRUE))
        sdv[is.na(sdv) | sdv == 0] <- 1
        
        dat_cv_scaled <- dat_cv
        for (tr in traits) dat_cv_scaled[[tr]] <- (dat_cv_scaled[[tr]] - mu[tr]) / sdv[tr]
        
        if (length(traits) == 1) {
          tr <- traits[1]
          dat_fit <- data.frame(id = dat_cv_scaled$id)
          dat_fit$y <- dat_cv_scaled[[tr]]
          
          fit_cv <- mmer(
            y ~ 1,
            random = ~ vsr(id, Gu = Gs[ids, ids]),
            rcov   = ~ units,
            data   = dat_fit,
            tolParInv = tolParInv
          )
        } else {
          form <- as.formula(paste0("cbind(", paste(traits, collapse = ","), ") ~ 1"))
          fit_cv <- mmer(
            form,
            random = ~ vsr(id, Gu = Gs[ids, ids]),
            rcov   = ~ units,
            data   = dat_cv_scaled,
            tolParInv = tolParInv
          )
        }
        
        gebv_cv <- extract_gebv(fit_cv, traits, fallback_ids = dat_cv_scaled$id)
        pred <- gebv_cv[test_ids, , drop = FALSE]
        
        obs <- dat_tr[match(test_ids, dat_tr$id), traits, drop = FALSE]
        for (tr in traits) obs[[tr]] <- (obs[[tr]] - mu[tr]) / sdv[tr]
        
        for (j in seq_along(traits)) {
          tr <- traits[j]
          acc_mat[f, tr] <- cor(pred[, tr], obs[[tr]], use = "complete.obs")
        }
      }
      
      acc_per_rep[[r]] <- acc_mat
    }
    
    all_acc <- do.call(rbind, acc_per_rep)
    cv_summary <- data.frame(
      trait = traits,
      mean  = as.numeric(colMeans(all_acc, na.rm = TRUE)),
      sd    = as.numeric(apply(all_acc, 2, sd, na.rm = TRUE))
    )
  }
  
  list(
    train_ids = common_train,
    geno_ids  = geno_ids,
    new_ids   = new_ids,
    GEBV = gebv_full,
    selection_index = index_tbl,
    cv_per_fold = all_acc,
    cv_summary = cv_summary,
    fit_full = fit_full
  )
}

# ===========================================
# pipeline OF GS plus CV plus prediction plus SAVE OUTPUTS
# ===========================================

run_single_and_multi_GS <- function(pheno, G_all, traits,
                                    id_col = "id",
                                    outdir = "GS_results",
                                    k = 5, reps = 3, seed = 123) {
  
  if (!dir.exists(outdir)) dir.create(outdir)
  
  # ---------------------------
  # 1) MULTI-TRAIT GS
  # ---------------------------
  cat("Running multi-trait GS...\n")
  
  res_mt <- run_gs_sommer(
    pheno = pheno,
    G = G_all,
    traits = traits,
    id_col = id_col,
    k = k, reps = reps,
    predict_new = TRUE
  )
  
  # save CV
  if (!is.null(res_mt$cv_summary)) {
    data.table::fwrite(res_mt$cv_summary,
                       file.path(outdir, "CV_multitrait.txt"),
                       sep="\t")
  }
  
  # save predictions
  gebv_mt <- res_mt$GEBV
  pred_mt <- data.frame(id = rownames(gebv_mt), gebv_mt)
  data.table::fwrite(pred_mt,
                     file.path(outdir, "Pred_multitrait_all.txt"),
                     sep="\t")
  
  if (length(res_mt$new_ids) > 0) {
    pred_mt_new <- pred_mt[pred_mt$id %in% res_mt$new_ids, ]
    data.table::fwrite(pred_mt_new,
                       file.path(outdir, "Pred_multitrait_new.txt"),
                       sep="\t")
  }
  
  # ---------------------------
  # 2) SINGLE-TRAIT GS (loop)
  # ---------------------------
  cat("Running single-trait GS...\n")
  
  cv_single_list <- list()
  pred_single_list <- list()
  
  for (tr in traits) {
    
    cat("  Trait:", tr, "\n")
    
    res_st <- run_gs_sommer(
      pheno = pheno,
      G = G_all,
      traits = c(tr),
      id_col = id_col,
      k = k, reps = reps,
      predict_new = TRUE
    )
    
    # CV
    if (!is.null(res_st$cv_summary)) {
      cv_single_list[[tr]] <- data.frame(
        trait = tr,
        mean = res_st$cv_summary$mean,
        sd   = res_st$cv_summary$sd
      )
    }
    
    # predictions
    gebv_st <- res_st$GEBV
    pred_single_list[[tr]] <- data.frame(
      id = rownames(gebv_st),
      trait = tr,
      GEBV = gebv_st[,1]
    )
  }
  
  # combine single-trait CV
  if (length(cv_single_list) > 0) {
    cv_single <- do.call(rbind, cv_single_list)
    data.table::fwrite(cv_single,
                       file.path(outdir, "CV_singletrait.txt"),
                       sep="\t")
  }
  
  # combine single-trait predictions to wide
  pred_single_long <- do.call(rbind, pred_single_list)
  pred_single_long <- data.table::as.data.table(pred_single_long)
  pred_single_wide <- data.table::dcast(pred_single_long,
                                        id ~ trait,
                                        value.var="GEBV")
  
  data.table::fwrite(pred_single_wide,
                     file.path(outdir, "Pred_singletrait_all.txt"),
                     sep="\t")
  
  list(
    multitrait = res_mt,
    singletrait = pred_single_wide
  )
}

# ===========================================
# pipeline OF GS plus CV plus SAVE OUTPUTS# ===========================================

gs_example_run_and_save <- function(pheno_file,
                                    vcf_file,
                                    traits,
                                    id_col = "id",
                                    outdir = "GS_results",
                                    # CV settings
                                    k = 5, reps = 3, seed = 123,
                                    # model settings
                                    tolParInv = 1e-2, nugget = 1e-6,
                                    scale_traits = TRUE,
                                    # prediction settings
                                    predict_new = TRUE) {
  
  gs_load_packages()
  library(data.table)
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # 1) read phenotype
  pheno <- read_pheno_txt(pheno_file)
  pheno[[id_col]] <- trimws(pheno[[id_col]])
  
  # 2) build GRM from VCF (must include training + new genotypes for prediction)
  grm_out <- build_grm_from_vcf(vcf_file)
  G_all <- grm_out$G
  rownames(G_all) <- trimws(rownames(G_all))
  colnames(G_all) <- trimws(colnames(G_all))
  
  # 3) CV evaluation (training only)
  res_cv <- run_gs_sommer(
    pheno = pheno,
    G = G_all,
    traits = traits,
    id_col = id_col,
    k = k, reps = reps, seed = seed,
    tolParInv = tolParInv, nugget = nugget,
    scale_traits = scale_traits,
    predict_new = FALSE
  )
  
  # save CV outputs
  if (!is.null(res_cv$cv_summary)) {
    fwrite(res_cv$cv_summary,
           file = file.path(outdir, "CV_summary.txt"),
           sep = "\t")
  }
  if (!is.null(res_cv$cv_per_fold)) {
    cv_fold_df <- as.data.frame(res_cv$cv_per_fold)
    cv_fold_df$FoldRep <- seq_len(nrow(cv_fold_df))
    fwrite(cv_fold_df,
           file = file.path(outdir, "CV_per_fold.txt"),
           sep = "\t")
  }
  
  # 4) FINAL FIT + prediction (training + new)
  res_pred <- run_gs_sommer(
    pheno = pheno,
    G = G_all,
    traits = traits,
    id_col = id_col,
    k = 0, reps = 0, seed = seed,
    tolParInv = tolParInv, nugget = nugget,
    scale_traits = scale_traits,
    predict_new = predict_new
  )
  
  # predictions for ALL genotypes
  gebv_all <- res_pred$GEBV
  pred_all <- data.frame(id = rownames(gebv_all), gebv_all, row.names = NULL)
  fwrite(pred_all,
         file = file.path(outdir, "Pred_GEBV_all_genotypes.txt"),
         sep = "\t")
  
  # predictions for NEW genotypes only
  new_ids <- res_pred$new_ids
  if (length(new_ids) > 0) {
    gebv_new <- gebv_all[new_ids, , drop = FALSE]
    pred_new <- data.frame(id = rownames(gebv_new), gebv_new, row.names = NULL)
    fwrite(pred_new,
           file = file.path(outdir, "Pred_GEBV_new_genotypes.txt"),
           sep = "\t")
  } else {
    message("No new genotypes detected (new_ids length = 0).")
  }
  
  # optional: save selection index if weights are used (not in this example)
  # if (!is.null(res_pred$selection_index)) {
  #   fwrite(res_pred$selection_index,
  #          file = file.path(outdir, "Selection_index_rank.txt"),
  #          sep = "\t")
  # }
  
  # return objects in case you want them in memory
  list(
    pheno = pheno,
    G_all = G_all,
    cv = res_cv,
    pred = res_pred
  )
}
