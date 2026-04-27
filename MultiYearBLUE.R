
library(lme4)
library(emmeans)

#estimate BLUE for multiple year traits (genotype is fixed and year is random)
get_BLUE_trait_yearRand_emm <- function(df, trait,
                                        id_col = "id",
                                        year_col = "year") {
  
  d <- df[, c(id_col, year_col, trait)]
  names(d) <- c("id", "year", "y")
  
  d$y <- as.numeric(d$y)
  d <- d[complete.cases(d), ]
  
  if (nrow(d) == 0) {
    message("Skipping ", trait, ": 0 non-NA cases")
    return(NULL)
  }
  
  d$id <- factor(d$id)
  d$year <- factor(d$year)
  
  if (length(unique(d$id)) < 2) {
    message("Skipping ", trait, ": fewer than 2 genotypes")
    return(NULL)
  }
  
  if (length(unique(d$year)) < 2) {
    message("Only one year for ", trait, ": using lm(y ~ id)")
    m <- lm(y ~ id, data = d)
  } else {
    m <- lmer(y ~ id + (1 | year), data = d)
  }
  
  emm <- as.data.frame(emmeans(m, ~ id))
  
  out <- data.frame(
    id = as.character(emm$id),
    blue_value = emm$emmean,
    row.names = NULL
  )
  
  colnames(out)[2] <- paste0(trait, "_BLUE")
  out
}

