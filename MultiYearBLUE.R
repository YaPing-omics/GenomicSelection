rm(list=ls())
setwd('C:/Users/ya-ping.lin/Documents/IMIN_GS/GS/MultiENV/India/')

library(lme4)
library(emmeans)

#estimate BLUE for multiple year traits (genotype is fixed and year is random)

get_BLUE_trait_yearRand_emm <- function(df, trait,
                                        id_col="id",
                                        year_col="year") {
  
  d <- df[, c(id_col, year_col, trait)]
  names(d) <- c("id","year","y")
  
  d <- d[complete.cases(d), ]
  d$y <- as.numeric(d$y)
  
  d$id   <- factor(d$id)
  d$year <- factor(d$year)
  
  # Model: genotype fixed, year random
  m <- lmer(y ~ id + (1|year), data = d)
  
  # Genotype BLUEs (estimated marginal means for id)
  emm <- as.data.frame(emmeans(m, ~ id))
  
  out <- data.frame(
    id = as.character(emm$id),
    blue_value = emm$emmean,
    row.names = NULL
  )
  colnames(out)[2] <- paste0(trait, "_BLUE")
  out
}



# =========================
# HOW TO RUN (EDIT PATHS)
# =========================
#convert yield component to a proper table

df <- read.delim('./IndiaYieldCom.txt', sep='\t')

df <- df %>%
  rename(id = Genotype) %>%
  pivot_longer(
    cols = -id,
    names_to = c("trait", "year"),
    names_pattern = "^(PodPerPlant|SW|SeedPerPod)_India_(\\d{4})\\.?\\d*",
    values_to = "value"
  ) %>%
  mutate(
    year = as.integer(year),
    value = as.numeric(value)
  ) %>%
  group_by(id, year, trait) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(value = ifelse(is.nan(value), NA, value)) %>%
  pivot_wider(
    names_from = trait,
    values_from = value
  ) %>%
  arrange(id, year)

head(df)

#convert yield to a proper table

df <- read.delim('./IndiaSeedYield.txt', sep='\t')

head(df)


# read file
df <- read.delim("IndiaSeedYield.txt", check.names = FALSE)

# reshape
df <- df %>%
  rename(id = Genotype) %>%
  pivot_longer(
    cols = -id,
    names_to = c("year", "rep"),
    names_pattern = "SeedYield_India_(\\d{4})_(\\d+)",
    values_to = "SeedYield"
  ) %>%
  mutate(
    year = as.integer(year),
    SeedYield = as.numeric(SeedYield)
  ) %>%
  group_by(id, year) %>%
  summarise(
    SeedYield = mean(SeedYield, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(SeedYield = ifelse(is.nan(SeedYield), NA, SeedYield)) %>%
  arrange(id, year)

# view
head(df)

#convert DMDF into a proper table
# read file
df <- read.delim("IndiaDMDF.txt", check.names = FALSE)

# reshape
df <- df %>%
  rename(id = Genotype) %>%
  pivot_longer(
    cols = -id,
    names_to = c("trait", "year"),
    names_pattern = "^(DM|DF)_India_(\\d{4})\\.?\\d*",
    values_to = "value"
  ) %>%
  mutate(
    year = as.integer(year),
    value = as.numeric(value)
  ) %>%
  group_by(id, year, trait) %>%
  summarise(
    value = mean(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(value = ifelse(is.nan(value), NA, value)) %>%
  pivot_wider(
    names_from = trait,
    values_from = value
  ) %>%
  arrange(id, year)

# view
head(df)

traits <- c("SeedPerPod","PodPerPlant","SW")
traits <- c("SeedYield")
traits <- c("DF","DM")


blue_list <- lapply(traits, function(tr)
  get_BLUE_trait_yearRand_emm(df, tr, id_col="id", year_col="year"))

pheno_BLUE <- Reduce(function(a,b) merge(a,b, by="id", all=TRUE), blue_list)

head(pheno_BLUE)

write.table(pheno_BLUE, 'India_yield_com_BLUE.txt', sep='\t', row.names=F, quote=F)

write.table(pheno_BLUE, 'India_seedyield_BLUE.txt', sep='\t', row.names=F, quote=F)

write.table(pheno_BLUE, 'India_DFDM_BLUE.txt', sep='\t', row.names=F, quote=F)

yield_com=read.delim('./India_yield_com_BLUE.txt')
seed_yield=read.delim('./India_seedyield_BLUE.txt')
dfdm=pheno_BLUE

which(yield_com[,1]==seed_yield[,1])
which(dfdm[,1]==seed_yield[,1])


d1=merge(yield_com, seed_yield, by='id')
d1=merge(dfdm, d1, by='id')

cor(d1[,-1], use = "pairwise.complete.obs")
write.table(d1, '5traits_BLUE.txt', sep='\t', row.names=F, quote=F)


library(corrplot)

# read your data
df <- read.delim("5traits_BLUE.txt", check.names = FALSE)

# remove id column and keep only numeric traits
num_df <- df %>%
  select(-id)

colnames(num_df) <- gsub("_BLUE", "", colnames(num_df))

# compute correlation matrix
cor_mat <- cor(num_df, use = "complete.obs")

# plot
pdf("corrplot_traits.pdf", width = 8, height = 6)

corrplot(
  cor_mat,
  method = "circle",
  type = "upper",
  order = "hclust",   # cluster similar traits
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45
)

dev.off()
