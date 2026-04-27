
library(emmeans)
library(lme4)

library(tidyr)
library(dplyr)
# ============================================================
# Genotype BLUEs with: Year RANDOM, Location FIXED, Genotype FIXED
# Uses emmeans to get clean genotype marginal means (BLUEs)
# Works trait-by-trait, then you can merge traits for multi-trait GS
# ============================================================

get_BLUE_trait_yearRand_locFixed_emm <- function(df, trait,
                                                 id_col="id",
                                                 year_col="year",
                                                 loc_col="loc") {
  d <- df[, c(id_col, year_col, loc_col, trait)]
  names(d) <- c("id","year","loc","y")
  
  d <- d[complete.cases(d), ]
  d$y <- as.numeric(d$y)
  
  d$id   <- factor(d$id)
  d$year <- factor(d$year)
  d$loc  <- factor(d$loc)
  
  # Model: id fixed, loc fixed, year random
  m <- lmer(y ~ id + loc + (1|year), data = d)
  
  # Genotype BLUEs (estimated marginal means for id, averaged over loc levels)
  emm <- as.data.frame(emmeans(m, ~ id))
  
  out <- data.frame(
    id = as.character(emm$id),
    blue_value = emm$emmean,
    row.names = NULL
  )
  colnames(out)[2] <- paste0(trait, "_BLUE")
  out
}

