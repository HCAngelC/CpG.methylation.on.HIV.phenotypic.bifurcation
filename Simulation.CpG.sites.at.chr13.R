################################################################################
##                                                                            ##
##   Experimental design:                                                     ##
##     • Jurkat / primary CD4+ T-cell clones, each harbouring ONE HIV-GFP    ##
##       provirus at a defined genomic integration site                       ##
##     • FACS-sorted GFP+ (active) and GFP- (latent) sub-populations         ##
##       from the same clone → matched RRBS, RT-qPCR                         ##
##     • Readouts: RRBS (LTR CpG methylation), RT-qPCR (HIV mRNA),           ##
##       FACS (GFP MFI, %GFP+)                                                ##
##                                                                            ##
##   Models:                                                                  ##
##     1.  Beta-binomial RRBS differential methylation (GFP+ vs GFP-)        ##
##     2.  CpG methylation index (CMI) computation per clone/sort             ##
##     3.  Hill-function: CMI → HIV transcription (multi-readout NLS)         ##
##     4.  ODE bistability: single-provirus epigenetic switch                 ##
##     5.  Markov switching model: latency ↔ reactivation rates               ##
##     6.  Linear mixed-effects regression (clone + CpG random effects)       ##
##     7.  Correlation structure: RRBS–RT-qPCR–FACS integration              ##
##     8.  Survival analysis: time-to-reactivation from epigenetic state      ##
##     9.  Stochastic Gillespie SSA: cell-to-cell heterogeneity               ##
##    10.  Publication-quality multi-panel figure (PDF)                       ##
##                                                                            ##
##   References:                                                              ##
##     Blazkova et al. (2009) PLoS Pathog – LTR methylation in latency       ##
##     Kauder et al. (2009) Nat Med – DNMT3a promotes HIV latency            ##
##     Weinberger et al. (2005) Cell – stochastic HIV transcription          ##
##     Jordan et al. (2001) EMBO J – clonal HIV latency model                ##
##     Trejbalova et al. (2016) Nucleic Acids Res – RRBS of LTR              ##
##                                                                            ##
##   Author : Heng-Chang Chen                                                 ##
##   Date   : 2026-05-25                                                      ##
##                                                                            ##
################################################################################

```{r 0 ▸ DEPENDENCIES}
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

set.seed(2024)

cat(rep("═", 66), "\n", sep = "")
cat("  HIV CLONAL EPIGENETIC MEMORY MODEL  –  Initialising\n")
cat(rep("═", 66), "\n\n", sep = "")
```
```{r 1 ▸ EXPERIMENTAL DESIGN & KNOWN BIOLOGY}
# HIV-1 CpG sites at chr13
# Based on HXB2 sequence; sites in core promoter & nucleosome-0 region
# References: Blazkova 2009; Kauder 2009; Trejbalova 2016
CPG_chr13 <- df.f814.all.CpG.input.chr13.loci

# Clonal model parameters
N_CLONES     <- 2   # independent clonal cell lines
N_CELLS_SORT <- 10000  # cells per sorted fraction per clone
N_CPG        <- nrow(CPG_chr13)
COVERAGE_RRBS <- 28.76  # mean RRBS read depth per CpG in f814 (unsorted)

# Clone metadata: integration site chromatin context
clone_meta <- tibble(
  clone_id        = c("Bright", "Dim"),
  lambda = c(log(2)/0.5, log(2)/0.5),
  amplitude_fractions = c(0.22, 0.78), # Biphasic amplitude fractions (A + B = 1)
  # Baseline %GFP+ (prior to sorting) – reflects epigenetic history
  pct_GFP_bulk    = NA_real_   # filled below after simulation
)

# Per-clone epigenetic "memory strength" parameter (ξ)
# Higher ξ → stronger epigenetic silencing memory → higher CpG methylation
xi_clone <- setNames(
  runif(N_CLONES, 0.3, 0.9),
  clone_meta$clone_id
)
```
```{r 2 ▸ SIMULATE RRBS DATA  (GFP+ and GFP- sorted populations per clone)}
# Model: methylation at CpG site j in clone c, population p
#   θ_{c,j,p} ~ Beta(α_{c,j,p}, β_{c,j,p})
#   observed reads ~ BetaBinomial(n = COVERAGE, θ_{c,j,p}, ρ)
#
# GFP- (latent)  : high methylation, especially at Sp1 sites
# GFP+ (active)  : low methylation at all sites

rbetabinom <- function(n, size, prob, rho = 0.05) {
  # Beta-binomial: overdispersed binomial (rho = overdispersion)
  a <- prob       * (1 - rho) / rho
  b <- (1 - prob) * (1 - rho) / rho
  p <- rbeta(n, shape1 = pmax(a, 0.01), shape2 = pmax(b, 0.01))
  rbinom(n, as.integer(size), p) # size should be an integer
}

simulate_rrbs_clone <- function(clone_id, xi = 0.5) {

  # Per-CpG baseline methylation for this subpopulation
  cpg_sensitivity <- ifelse(CPG_chr13$sig == "sig", 1.4, 1.0)
  
  meth_latent <- pmin(0.05 + xi * 0.85 * cpg_sensitivity *
                        (1 + rnorm(N_CPG, 0, 0.07)), 0.98)
  meth_active <- pmin(0.02 + (1 - xi) * 0.25 * cpg_sensitivity *
                        (1 + rnorm(N_CPG, 0, 0.05)), 0.35)

  bind_rows(
    tibble(
      index     = CPG_chr13$index,
      sort_pop   = "GFP_bright",
      site       = CPG_chr13$start,
      beta       = CPG_chr13$beta.bright,
      dm         = CPG_chr13$dm.bright,
      meth_prob  = meth_active,
      meth_reads = rbetabinom(N_CPG, COVERAGE_RRBS, meth_active, rho = 0.04),
      total_reads = COVERAGE_RRBS,
      meth_pct   = meth_reads / total_reads,
      rank       = CPG_chr13$Rank,
      sig        = CPG_chr13$sig
    ),
    tibble(
      index     = CPG_chr13$index,
      sort_pop   = "GFP_dim",
      site       = CPG_chr13$start,
      beta       = CPG_chr13$beta.dim,
      dm         = CPG_chr13$dm.dim,
      meth_prob  = meth_latent,
      meth_reads = rbetabinom(N_CPG, COVERAGE_RRBS, meth_latent, rho = 0.06),
      total_reads = COVERAGE_RRBS,
      meth_pct   = meth_reads / total_reads,
      rank       = CPG_chr13$Rank,
      sig        = CPG_chr13$sig
    )
  )
}

rrbs <- bind_rows(
  pmap_dfr(list(clone_meta$clone_id, xi_clone[clone_meta$clone_id]),
           simulate_rrbs_clone)
) %>%
  mutate(
    sort_pop = factor(sort_pop, levels = c("GFP_bright","GFP_dim"),
                      labels = c("GFP_bright","GFP_dim"))
  )

cat("✓ RRBS simulated:", nrow(rrbs), "observations |",
    N_CLONES, "clones ×", N_CPG, "CpGs × 2 sorted populations\n")
```
```{r 3 ▸ CpG METHYLATION INDEX (CMI)  –  per clone per population}
# CMI_c = weighted mean methylation across CpGs at chr13
# Weights: coverage × significant sites (functionally significant sites weighted 2×)
CMI <- rrbs %>%
  mutate(weight = total_reads * ifelse(sig == "sig", 2, 1)) %>%
  group_by(index, sort_pop) %>%
  summarise(
    CMI        = sum(meth_reads * weight) / sum(total_reads * weight),
    CMI_unweighted = mean(meth_pct),
    n_cpg      = n(),
    .groups    = "drop"
    ) %>%
  left_join(CPG_chr13, by = "index")

cat("✓ CpG Methylation Index (CMI) computed per clone/population\n")
```
```{r 4 ▸ MODEL 1 –  BETA-BINOMIAL DIFFERENTIAL METHYLATION  (GFP+ vs GFP-)}
# For each CpG site: test whether methylation differs between GFP+ and GFP-
# Using beta-regression on the proportion (per-site, pooled across clones)

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

cat("── Beta-Regression Differential Methylation ───────────────────\n")
cat("  Significant CpGs (FDR < 5%):",
    sum(beta_tests$significant, na.rm = TRUE), "/", N_CPG, "\n")
cat("  Hypermethylated in GFP⁻ (latent):",
    sum(beta_tests$significant & beta_tests$direction == "Hyper in GFP⁻",
        na.rm = TRUE), "\n\n")
```
