rm(list=ls())

source("gs_pipeline_singletrait.R")

#single traits
# ---- Example call GS + CV + prediction (edit filenames/traits) ----

run_singletrait_marker_number_experiment(
  pheno_file = "5traits_BLUE.txt",
  vcf_file = "MMC_MAF005.vcf",
  traits = c("DF", "DM","SeedPerPod",'PodPerPlant',"SW",'SeedYield'),
  id_col = "id",
  marker_numbers = c(1000, 5000, 10000, 20000),
  marker_reps = 50,
  base_outdir = "GS_singletrait",
  k = 5,
  reps = 3,
  seed = 123,
  marker_seed_base = 1000
)

#multiple traits
# ---- Example call GS + CV + prediction (edit filenames/traits) ----

source("gs_pipeline_multitrait.R")

run_multitrait_marker_number_experiment(
  pheno_file = "5traits_BLUE.txt",
  vcf_file = "MMC_MAF005.vcf",
  traits = c('PodPerPlant',"SW",'SeedYield'),
  id_col = "id",
  marker_numbers = c(5000),
  base_outdir = "GS_multitrait",
  marker_reps = 10,
  k = 5,
  reps = 3,
  seed = 123,
  marker_seed_base = 1000
)
