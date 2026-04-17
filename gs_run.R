rm(list=ls())
setwd('C:/Users/ya-ping.lin/Documents/IMIN_GS/GS/MultiENV/')

source("../gs_pipeline.R")

# ---- Example call GS + CV + prediction (edit filenames/traits) ----

pheno <- read_pheno_txt("yield_BLUE.txt")
grm   <- build_grm_from_vcf("../../MMC_MAF005.vcf")
G <- grm$G

traits <- c("PodPerPlant_BLUE","SeedPerPod_BLUE","SW_BLUE")
outdir <- "yield"

out <- run_single_and_multi_GS(
  pheno = pheno,
  G_all = G,
  traits = traits,
  id_col = "id",
  outdir = outdir,
  k = 5,
  reps = 3,
  seed = 123
)

# Example call GS + CV only (edit filenames/traits) 
out <- gs_example_run_and_save(
  pheno_file = "SW_BLUE.txt",
  vcf_file   = "../../MMC_MAF005.vcf",
  traits     = c("SW_BLUE"),
  id_col     = "id",
  outdir     = "SW",
  k = 5, reps = 10, seed = 123
)
