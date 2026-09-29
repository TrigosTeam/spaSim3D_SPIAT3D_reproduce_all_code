library(dplyr)

corr_df <- read.csv("~/R/values_from_figures/showing_structure_corr_df.csv")

structures <- c(
  "mixed_ellipsoid",
  "separated_ellipsoid",
  "ringed_ellipsoid",
  "mixed_network",
  "separated_network",
  "ringed_network"
)

# slope of line of best fit
trend_df <- corr_df %>%
  mutate(
    structure_order = match(structure, structures)
  ) %>%
  group_by(metric) %>%
  summarise(
    slope = coef(lm(corr ~ structure_order))[2],
    r_squared = summary(lm(corr ~ structure_order))$r.squared,
    p_value = summary(lm(corr ~ structure_order))$coefficients[2, 4],
    .groups = "drop"
  )


# spearman's rank correlation

trend_df <- corr_df %>%
  mutate(structure_order = match(structure, structures)) %>%
  group_by(metric) %>%
  summarise(
    slope = coef(lm(corr ~ structure_order))[2],
    spearman_rho = cor(
      structure_order,
      corr,
      method = "spearman",
      use = "complete.obs"
    ),
    .groups = "drop"
  )




