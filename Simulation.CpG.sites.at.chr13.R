################################################################################
##                                                                            ##
##   Experimental design:                                                     ##
##     • Single provirus-integrated HIV transcription sensitized cellular     ##
##       model.                                                               ##
##     • FACS-sorted GFP-Bright and GFP-Dim subpopulations                    ##
##     • Readouts: RRBS (chr13)                                               ##
##                                                                            ##
##   Models:                                                                  ##
##     1.  Beta-binomial RRBS differential methylation (GFP+ vs GFP-)         ##
##   Author : Heng-Chang Chen                                                 ##
##   Date   : 2026-05-25                                                      ##
##                                                                            ##
################################################################################

# 0 ▸ DEPENDENCIES
pkgs <- c(
  # Core
  "dplyr", "tidyr", "tibble", "purrr", "readr", "stringr",
  # Modelling
  "deSolve",       # ODE integration
  "nlme",          # non-linear mixed effects
  "lme4",          # linear mixed-effects
  "lmerTest",      # p-values for lmer
  "betareg",       # beta regression for methylation proportions
  "VGAM",          # beta-binomial distribution
  "survival",      # Kaplan-Meier / Cox
  "survminer",     # ggplot2 survival plots
  "GillespieSSA2", # stochastic SSA
  # Visualisation
  "ggplot2", "patchwork", "ggbeeswarm", "ggridges",
  "viridis", "scales", "ggrepel", "RColorBrewer",
  # Utilities
  "broom", "broom.mixed", "knitr"
)

invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}))

# 1 ▸ EXPERIMENTAL DESIGN & KNOWN BIOLOGY
# CpG sites at chr13

CpG_chr13 # input file

# Model parameters
N_SUBPOPULATIONS     <- 2   # independent GFP-Bright and GFP-Dim
N_CELLS_SORT <- 10000  # cells per sorted fraction per subpopulation
N_CpG        <- nrow(CpG_chr13)
COVERAGE_RRBS <- 28.76  # mean RRBS read depth per CpG in Single provirus-integrated HIV transcription sensitized cellular model (unsorted)

# Subpopulation metadata
subpopulation_meta <- tibble(
  subpopulation_id        = c("Bright", "Dim"),
  lambda = c(log(2)/0.5, log(2)/0.5),
  amplitude_fractions = c(0.22, 0.78), # Biphasic amplitude fractions (A + B = 1)
  # Baseline %GFP+ (prior to sorting) – reflects epigenetic history
  pct_GFP_bulk    = NA_real_   # filled below after simulation
)

# Epigenetic "memory strength" parameter (ξ)
# Higher ξ → stronger epigenetic silencing memory → higher CpG methylation
xi_subpopulation <- setNames(
  runif(N_SUBPOPULATIONS, 0.3, 0.9),
  subpopulation_meta$subpopulation_id
)

# 2 ▸ SIMULATE RRBS DATA (GFP-Bright and GFP-Dim sorted subpopulations)
# =============================================================
# Model: methylation at CpG site j across 2 subpopulations (GFP-Bright versus GFP-Dim)
#
#   θ_{j,p} ~ Beta(α_{j,p}, β_{j,p})
#   observed reads ~ BetaBinomial(n = CpG sites, θ_{j,p}, ρ)
#
# GFP-Dim : high methylation
# GFP-Bright : low methylation at all sites
#
# =============================================================

rbetabinom <- function(n, size, prob, rho = 0.05) {
  # Beta-binomial: overdispersed binomial (rho = overdispersion)
  a <- prob       * (1 - rho) / rho
  b <- (1 - prob) * (1 - rho) / rho
  p <- rbeta(n, shape1 = pmax(a, 0.01), shape2 = pmax(b, 0.01))
  rbinom(n, as.integer(size), p) # size should be an integer
}

simulate_rrbs_subpopulation <- function(subpopulation_id, xi = 0.5) {

  # Per-CpG baseline methylation for this subpopulation
  cpg_sensitivity <- ifelse(CpG_chr13$sig == "sig", 1.4, 1.0)
  
  meth_latent <- pmin(0.05 + xi * 0.85 * cpg_sensitivity *
                        (1 + rnorm(N_CpG, 0, 0.07)), 0.98)
  meth_active <- pmin(0.02 + (1 - xi) * 0.25 * cpg_sensitivity *
                        (1 + rnorm(N_CpG, 0, 0.05)), 0.35)

  bind_rows(
    tibble(
      index     = CpG_chr13$index,
      sort_pop   = "GFP_bright",
      site       = CpG_chr13$start,
      beta       = CpG_chr13$beta.bright,
      dm         = CpG_chr13$dm.bright,
      meth_prob  = meth_active,
      meth_reads = rbetabinom(N_CpG, COVERAGE_RRBS, meth_active, rho = 0.04),
      total_reads = COVERAGE_RRBS,
      meth_pct   = meth_reads / total_reads,
      rank       = CpG_chr13$Rank,
      sig        = CpG_chr13$sig
    ),
    tibble(
      index     = CpG_chr13$index,
      sort_pop   = "GFP_dim",
      site       = CpG_chr13$start,
      beta       = CpG_chr13$beta.dim,
      dm         = CpG_chr13$dm.dim,
      meth_prob  = meth_latent,
      meth_reads = rbetabinom(N_CpG, COVERAGE_RRBS, meth_latent, rho = 0.06),
      total_reads = COVERAGE_RRBS,
      meth_pct   = meth_reads / total_reads,
      rank       = CpG_chr13$Rank,
      sig        = CpG_chr13$sig
    )
  )
}

rrbs <- bind_rows(
  pmap_dfr(list(subpopulation_meta$subpopulation_id, xi_subpopulation[subpopulation_meta$subpopulation_id]),
           simulate_rrbs_subpopulation)
) %>%
  mutate(
    sort_pop = factor(sort_pop, levels = c("GFP_bright","GFP_dim"),
                      labels = c("GFP_bright","GFP_dim"))
  )

# 3 ▸ MODEL –  BETA-BINOMIAL DIFFERENTIAL METHYLATION  (GFP-Bright versus GFP-Dim)
# For each CpG site: test whether methylation differs between GFP-Bright and GFP-Dim
# Using beta-regression on the proportion

site_summary <- rrbs %>%
  group_by(index, site, rank, sig) %>%
  summarise(
    mean_meth = mean(meth_pct),
    sd_meth   = sd(meth_pct),
    n         = n(),
    .groups   = "drop"
  )

# Beta-regression per CpG site
beta_tests <- rrbs %>%
  # Ensure proportions strictly within (0,1) for beta regression
  mutate(meth_pct_adj = pmin(pmax(meth_pct, 0.001), 0.999)) %>%
  group_by(index, site, rank, sig) %>%
  group_modify(~{
    tryCatch({
      fit <- betareg(meth_pct_adj ~ sort_pop,
                     data = .x, link = "logit")
      cf  <- coef(summary(fit))$mean
      tibble(
        estimate = cf["sort_popGFP_dim", "Estimate"],
        std_err  = cf["sort_popGFP_dim", "Std. Error"],
        z_val    = cf["sort_popGFP_dim", "z value"],
        p_val    = cf["sort_popGFP_dim", "Pr(>|z|)"]
      )
    }, error = function(e) {
      tibble(estimate = NA, std_err = NA, z_val = NA, p_val = NA)
    })
  }) %>%
  ungroup() %>%
  mutate(
    p_adj       = p.adjust(p_val, method = "BH"),
    significant = p_adj < 0.05,
    direction   = ifelse(estimate < 0, "Hyper in GFP_dim", "Hyper in GFP_bright")
  )
