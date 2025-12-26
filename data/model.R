# ==============================================================================
# CAUSAL MEDIATION (INTERVENTIONAL DIRECT / INDIRECT EFFECTS)
# Multi-valued treatment (Community) vs reference, binary outcomes
#
# Estimands (Risk Difference scale), for each community a != a0:
#   theta_aa = E[ Y^{a, G^a} ]   (equals E[Y^a] under consistency)
#   theta_a0 = E[ Y^{a, G^{a0}} ]
#   theta_00 = E[ Y^{a0, G^{a0}} ]
#   Total Effect (TE)  = theta_aa - theta_00
#   Interventional Direct Effect (IDE) = theta_a0 - theta_00
#   Interventional Indirect Effect (IIE) = theta_aa - theta_a0
#
# Identification: assumes no unmeasured confounding for A->(M,Y) given L and for M->Y given (A,L),
# plus positivity; uses *interventional* mediator distributions G^a.
#
# Features:
#   - Filter only communities with n>=100, plus pre-justified exclusion of C3_1 (positivity analysis)
#   - Rare-category lumping for baseline categorical covariates (explicit "Missing")
#   - Cross-fitting (K folds) to reduce overfitting bias
#   - Monte Carlo integration for mediator draws (joint via sequential factorization)
#   - Mediator simulation:
#       * Auto-detect bounded/percentage-like mediators -> logit transform on [lo, hi]
#       * Highly skewed positive mediators -> log transform (+ shift if needed)
#       * Residual bootstrap on transformed scale (instead of Gaussian noise), then inverse-transform
#   - Cluster-robust uncertainty via community-level Bayesian (weighted) bootstrap (Exp(1) weights)
#
# INPUT:  communities.csv
# OUTPUT: Causal_Mediation_Output/
#   - Full results CSV
#   - Figure PDF/PNG
# ==============================================================================

# -----------------------------
# 1) SETUP
# -----------------------------
setwd(".")
pkgs <- c("dplyr","tibble","readr","forcats","splines","ggplot2","stringr","tidyr","patchwork","purrr")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

DATA_PATH     <- "data.csv"
REF_COMMUNITY <- "C1_2"
SEED_VAL      <- 1234

KFOLDS        <- 4
MC_DRAWS      <- 100     # Increase to 500-1000 for final
BOOT_REPS     <- 100     # Increase to 500-1000 for final

MIN_LEVEL_N   <- 30      # Rare category threshold for baseline factors

dir.create("Causal_Mediation_Output", showWarnings = FALSE)

# -----------------------------
# 2) CONFOUNDERS & MEDIATORS
# -----------------------------
#BASE_VARS <- c("GENDER", "RACE", "ETHNICITY", "Ancestry", "AGE_imp", "AGE_miss")
BASE_VARS <- c("GENDER", "AGE_imp", "AGE_miss", "RACE")

MEDIATOR_CONFIG <- list(
  "Outcome_Asthma"     = c("A27_PM25", "A1_PPOV"),
  "Outcome_Dermatitis" = c("A27_PM25", "A1_PPOV"),
  "Outcome_Epilepsy"   = c("A17_UNEMPLOY", "A10_CRD_OBESITY", "A1_PPOV"),
  "Outcome_Influenza"  = c("A10_CRD_OBESITY", "A12_CRD_CSMOKING", "A1_PPOV")
)

OUTCOME_MAP <- list(
  Outcome_Asthma     = c(col="Asthma",            yes="Yes"),
  Outcome_Dermatitis = c(col="Dermatitis",        yes="Yes"),
  Outcome_Epilepsy   = c(col="Epilepsy",          yes="Yes"),
  Outcome_Influenza  = c(col="Influenza virus",   yes="Yes")
)

AGE_COL <- "AGE_AT_LAST_HS_VISIT"

# -----------------------------
# 3) HELPERS
# -----------------------------
clean_factor <- function(x, min_n = MIN_LEVEL_N) {
  f <- factor(x)
  f <- forcats::fct_explicit_na(f, na_level = "Missing")
  f <- forcats::fct_lump_min(f, min = min_n, other_level = "Other")
  if (!("Missing" %in% levels(f)) && any(is.na(x))) {
    f <- forcats::fct_explicit_na(f, na_level = "Missing")
  }
  f
}

median_imp_flag <- function(x) {
  xr <- suppressWarnings(as.numeric(x))
  miss <- as.integer(is.na(xr))
  imp <- ifelse(is.na(xr), stats::median(xr, na.rm = TRUE), xr)
  list(imp = as.numeric(imp), miss = as.integer(miss))
}

make_stratified_folds <- function(A, K = 5, seed = 1) {
  set.seed(seed)
  A <- factor(A)
  folds <- rep(NA_integer_, length(A))
  for (lvl in levels(A)) {
    idx <- which(A == lvl)
    idx <- sample(idx)
    folds[idx] <- rep(1:K, length.out = length(idx))
  }
  folds
}

rhs_with_splines <- function(vars, data, spline_df = 3) {
  out <- character(0)
  for (v in vars) {
    if (!v %in% names(data)) next
    if (is.numeric(data[[v]])) {
      u <- length(unique(data[[v]]))
      if (is.finite(u) && u >= 6) out <- c(out, paste0("splines::ns(", v, ", df=", spline_df, ")"))
      else out <- c(out, v)
    } else {
      out <- c(out, v)
    }
  }
  out
}

align_levels <- function(d, level_map) {
  for (nm in names(level_map)) {
    if (!nm %in% names(d)) next
    if (is.factor(d[[nm]]) || is.character(d[[nm]])) {
      d[[nm]] <- factor(as.character(d[[nm]]), levels = level_map[[nm]])
    }
  }
  d
}

add_zero_weight_rows_for_factors <- function(d_tr, factor_vars, level_map) {
  d_aug <- d_tr
  if (!(".w" %in% names(d_aug))) d_aug$.w <- 1
  template <- d_aug[1, , drop = FALSE]
  template$.w <- 0
  
  for (v in factor_vars) {
    if (!v %in% names(d_aug)) next
    if (!(is.factor(d_aug[[v]]) || is.character(d_aug[[v]]))) next
    
    present <- levels(factor(as.character(d_aug[[v]])))
    all_lv  <- level_map[[v]]
    miss <- setdiff(all_lv, present)
    miss <- miss[!is.na(miss)]
    if (length(miss) == 0) next
    
    for (lv in miss) {
      row <- template
      row[[v]] <- factor(lv, levels = level_map[[v]])
      row <- align_levels(row, level_map)
      d_aug <- dplyr::bind_rows(d_aug, row)
    }
  }
  d_aug
}

# Skewness helper (base R)
skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 5) return(0)
  m <- mean(x); s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(0)
  mean((x - m)^3) / (s^3)
}

# Transform specs for mediators
# Returns a list with:
#   type: "logit" | "log" | "identity"
#   lo, hi: bounds (for logit / clipping)
#   shift: for log (to make positive)
#   fwd(x), inv(z): transform and inverse
make_mediator_spec <- function(x, name = "", pct_like = FALSE) {
  x0 <- x[is.finite(x)]
  lo <- min(x0); hi <- max(x0)
  rng <- hi - lo
  sk <- skewness(x0)
  
  # Heuristic: percent/bounded if within [0,100] and has finite range
  bounded_0_100 <- (lo >= 0 && hi <= 100 && rng > 0)
  
  # Decide type:
  # - Very skewed positive -> log
  # - Else if bounded -> logit on padded bounds (keeps within plausible range)
  # - Else identity
  type <- "identity"
  shift <- 0
  
  if (lo > 0 && sk > 1.2) {
    type <- "log"
  } else if (bounded_0_100) {
    type <- "logit"
  } else if (lo >= 0 && sk > 1.5) {
    type <- "log"
    shift <- 1e-3
  }
  
  # Bounds padding for logit (avoid exact 0/1)
  if (type == "logit") {
    pad <- 0.01 * rng
    lo2 <- max(0, lo - pad)
    hi2 <- min(100, hi + pad)
    if (hi2 <= lo2) { lo2 <- lo; hi2 <- hi }
    lo <- lo2; hi <- hi2
  } else {
    # For log/identity, still keep a soft clipping window to prevent explosion
    pad <- 0.02 * rng
    lo <- lo - pad
    hi <- hi + pad
  }
  
  eps <- 1e-6
  
  if (type == "logit") {
    fwd <- function(v) {
      u <- (v - lo) / (hi - lo)
      u <- pmin(pmax(u, eps), 1 - eps)
      stats::qlogis(u)
    }
    inv <- function(z) {
      u <- stats::plogis(z)
      lo + (hi - lo) * u
    }
  } else if (type == "log") {
    # ensure positivity
    minx <- min(x0)
    if (minx <= 0) shift <- abs(minx) + 1e-3
    fwd <- function(v) log(v + shift)
    inv <- function(z) exp(z) - shift
  } else {
    fwd <- function(v) v
    inv <- function(z) z
  }
  
  list(type = type, lo = lo, hi = hi, shift = shift, skew = sk,
       fwd = fwd, inv = inv)
}

# Residual bootstrap simulator on transformed scale
# Returns: list(draws) where each element is a data.frame with updated mediator columns in original scale
simulate_mediators_residboot <- function(base_dat, a_prime,
                                         treat, lvls,
                                         meds_imp, meds_miss,
                                         L_vars,
                                         med_models, med_specs, med_resids,
                                         mc_draws = 200) {
  
  out <- base_dat
  out[[treat]] <- factor(a_prime, levels = lvls)
  
  draws <- vector("list", mc_draws)
  
  # Pre-sample residual indices for each mediator and draw to reduce overhead
  res_idx <- lapply(meds_imp, function(m) sample.int(length(med_resids[[m]]), size = nrow(out) * mc_draws, replace = TRUE))
  names(res_idx) <- meds_imp
  
  for (s in 1:mc_draws) {
    cur <- out
    for (j in seq_along(meds_imp)) {
      m_imp <- meds_imp[j]
      fit_m <- med_models[[m_imp]]
      spec  <- med_specs[[m_imp]]
      rr    <- med_resids[[m_imp]]
      
      mu_t <- as.numeric(stats::predict(fit_m, newdata = cur, type = "response"))
      
      # residual bootstrap in transformed space
      idx <- res_idx[[m_imp]][((s-1)*nrow(cur) + 1):(s*nrow(cur))]
      r_s <- rr[idx]
      z_t <- mu_t + r_s
      
      m_new <- spec$inv(z_t)
      
      # hard clip to spec bounds (important esp. for log/identity)
      m_new <- pmin(pmax(m_new, spec$lo), spec$hi)
      
      cur[[m_imp]] <- as.numeric(m_new)
    }
    draws[[s]] <- cur
  }
  draws
}

# -----------------------------
# 4) DATA CLEANING
# -----------------------------
cat(">>> LOADING DATA...\n")
df <- readr::read_csv(DATA_PATH, show_col_types = FALSE)

df_clean <- df %>%
  dplyr::group_by(Community) %>%
  dplyr::filter(dplyr::n() >= 100) %>%
  dplyr::ungroup() %>%
  dplyr::filter(Community != "C3_1") %>%
  dplyr::mutate(
    Community = factor(Community),
    GENDER    = clean_factor(GENDER),
    RACE      = clean_factor(RACE),
    ETHNICITY = clean_factor(ETHNICITY),
    Ancestry  = clean_factor(Ancestry)
  )

age_h <- median_imp_flag(df_clean[[AGE_COL]])
df_clean$AGE_imp  <- age_h$imp
df_clean$AGE_miss <- age_h$miss

for (nm in names(OUTCOME_MAP)) {
  col <- OUTCOME_MAP[[nm]][["col"]]
  yes <- OUTCOME_MAP[[nm]][["yes"]]
  if (!col %in% names(df_clean)) stop(paste("Missing outcome source column:", col))
  df_clean[[nm]] <- dplyr::if_else(df_clean[[col]] == yes, 1L, 0L, missing = NA_integer_)
}

all_meds <- unique(unlist(MEDIATOR_CONFIG))
for (m in all_meds) {
  if (!m %in% names(df_clean)) stop(paste("Missing mediator column:", m))
  h <- median_imp_flag(df_clean[[m]])
  df_clean[[paste0(m, "_imp")]]  <- h$imp
  df_clean[[paste0(m, "_miss")]] <- h$miss
}

df_clean <- df_clean %>% droplevels()

cat(">>> Data Ready. N =", nrow(df_clean), "| Communities =", length(levels(df_clean$Community)), "\n")

# -----------------------------
# 5) CORE ESTIMATOR: CROSS-FIT INTERVENTIONAL EFFECTS
# -----------------------------
estimate_interventional_effects <- function(d, outcome, mediators,
                                            L_vars = BASE_VARS,
                                            treat = "Community",
                                            ref = REF_COMMUNITY,
                                            K = KFOLDS, seed = SEED_VAL,
                                            mc_draws = MC_DRAWS,
                                            w = NULL,
                                            verbose_specs = FALSE) {
  
  set.seed(seed)
  
  meds_imp  <- paste0(mediators, "_imp")
  meds_miss <- paste0(mediators, "_miss")
  
  needed <- unique(c(outcome, treat, L_vars, meds_imp, meds_miss))
  dd <- d %>%
    dplyr::select(dplyr::all_of(needed)) %>%
    dplyr::filter(!is.na(.data[[outcome]])) %>%
    droplevels()
  
  if (is.null(w)) w <- rep(1, nrow(dd))
  dd$.w <- as.numeric(w)
  
  dd[[treat]] <- factor(dd[[treat]])
  if (ref %in% levels(dd[[treat]])) dd[[treat]] <- stats::relevel(dd[[treat]], ref = ref)
  lvls <- levels(dd[[treat]])
  if (!(ref %in% lvls)) stop("Reference community not present after filtering.")
  cmp_lvls <- setdiff(lvls, ref)
  
  # stable factor levels
  level_map <- list()
  level_map[[treat]] <- lvls
  factor_vars <- c(treat)
  
  for (v in needed) {
    if (!v %in% names(dd)) next
    if (is.factor(dd[[v]]) || is.character(dd[[v]])) {
      dd[[v]] <- factor(dd[[v]])
      level_map[[v]] <- levels(dd[[v]])
      factor_vars <- unique(c(factor_vars, v))
    }
  }
  dd <- align_levels(dd, level_map)
  
  n <- nrow(dd)
  folds <- make_stratified_folds(dd[[treat]], K = K, seed = seed)
  
  pred_00 <- rep(NA_real_, n)
  pred_aa <- matrix(NA_real_, n, length(cmp_lvls), dimnames = list(NULL, cmp_lvls))
  pred_a0 <- matrix(NA_real_, n, length(cmp_lvls), dimnames = list(NULL, cmp_lvls))
  
  # Outcome model formula
  f_q <- stats::as.formula(paste0(
    outcome, " ~ ",
    paste(rhs_with_splines(c(treat, L_vars, meds_imp, meds_miss), dd), collapse = " + ")
  ))
  
  # fold loop
  for (k in 1:K) {
    te <- which(folds == k)
    tr <- which(folds != k)
    if (length(te) == 0 || length(tr) == 0) next
    
    d_tr0 <- dd[tr, , drop = FALSE]
    d_te  <- dd[te, , drop = FALSE]
    
    d_tr0 <- align_levels(d_tr0, level_map)
    d_te  <- align_levels(d_te,  level_map)
    
    d_tr <- add_zero_weight_rows_for_factors(d_tr0, factor_vars = factor_vars, level_map = level_map)
    
    # ---- Mediator specs + models (sequential factorization) ----
    med_specs  <- list()
    med_models <- list()
    med_resids <- list()
    
    prev <- character(0)
    for (j in seq_along(meds_imp)) {
      m_imp <- meds_imp[j]
      
      # Decide transform spec from training data distribution
      spec <- make_mediator_spec(d_tr[[m_imp]], name = m_imp)
      med_specs[[m_imp]] <- spec
      
      # Create transformed response in training fold
      d_tr[[paste0(m_imp, "_t")]] <- spec$fwd(d_tr[[m_imp]])
      
      # Build conditional formula: T(Mj) ~ A + L + previous M_imp + miss flags
      rhs <- rhs_with_splines(c(treat, L_vars, prev, meds_miss), d_tr)
      fm  <- stats::as.formula(paste0(paste0(m_imp, "_t"), " ~ ", paste(rhs, collapse = " + ")))
      
      # Fit model (Gaussian)
      fit_m <- stats::glm(fm, data = d_tr, family = stats::gaussian(), weights = .w)
      med_models[[m_imp]] <- fit_m
      
      # Residuals in transformed space for bootstrap simulation
      mu_hat <- as.numeric(stats::predict(fit_m, newdata = d_tr, type = "response"))
      res <- d_tr[[paste0(m_imp, "_t")]] - mu_hat
      # Keep only non-zero weights residuals to avoid pseudo rows
      res <- res[d_tr$.w > 0]
      if (length(res) < 10) res <- res[is.finite(res)]
      med_resids[[m_imp]] <- res
      
      if (verbose_specs && k == 1) {
        cat(sprintf("Mediator %s: type=%s, skew=%.2f, lo=%.3f, hi=%.3f\n",
                    m_imp, spec$type, spec$skew, spec$lo, spec$hi))
      }
      
      prev <- c(prev, m_imp) # conditioning uses original-scale mediators
    }
    
    # ---- Outcome model ----
    fit_q <- stats::glm(f_q, data = d_tr, family = stats::binomial(), weights = .w)
    
    # ---- Simulate mediator draws under ref and each a ----
    draws_ref <- simulate_mediators_residboot(
      base_dat   = d_te, a_prime = ref,
      treat = treat, lvls = lvls,
      meds_imp = meds_imp, meds_miss = meds_miss,
      L_vars = L_vars,
      med_models = med_models, med_specs = med_specs, med_resids = med_resids,
      mc_draws = mc_draws
    )
    
    draws_by_a <- lapply(cmp_lvls, function(a) simulate_mediators_residboot(
      base_dat   = d_te, a_prime = a,
      treat = treat, lvls = lvls,
      meds_imp = meds_imp, meds_miss = meds_miss,
      L_vars = L_vars,
      med_models = med_models, med_specs = med_specs, med_resids = med_resids,
      mc_draws = mc_draws
    ))
    names(draws_by_a) <- cmp_lvls
    
    # Average Q(a, draws)
    avg_Q <- function(a_for_Y, draws_list) {
      qsum <- rep(0, nrow(d_te))
      for (s in 1:length(draws_list)) {
        nd <- draws_list[[s]]
        nd[[treat]] <- factor(a_for_Y, levels = lvls)
        qsum <- qsum + as.numeric(stats::predict(fit_q, newdata = nd, type = "response"))
      }
      qsum / length(draws_list)
    }
    
    pred_00[te] <- avg_Q(ref, draws_ref)
    
    for (a in cmp_lvls) {
      pred_a0[te, a] <- avg_Q(a, draws_ref)
      pred_aa[te, a] <- avg_Q(a, draws_by_a[[a]])
    }
  }
  
  # keep only rows with predictions
  ok <- which(!is.na(pred_00))
  if (length(ok) < n) {
    dd <- dd[ok, , drop = FALSE]
    pred_00 <- pred_00[ok]
    pred_aa <- pred_aa[ok, , drop = FALSE]
    pred_a0 <- pred_a0[ok, , drop = FALSE]
  }
  
  W <- dd$.w
  Wtot <- sum(W)
  
  theta_00 <- sum(W * pred_00) / Wtot
  
  tibble::tibble(
    treat_level = cmp_lvls,
    theta_00 = theta_00,
    theta_a0 = vapply(cmp_lvls, function(a) sum(W * pred_a0[, a]) / Wtot, numeric(1)),
    theta_aa = vapply(cmp_lvls, function(a) sum(W * pred_aa[, a]) / Wtot, numeric(1))
  ) %>%
    dplyr::mutate(
      TE  = theta_aa - theta_00,
      IDE = theta_a0 - theta_00,
      IIE = theta_aa - theta_a0
    )
}

# -----------------------------
# 6) CLUSTER-ROBUST WEIGHTED BOOTSTRAP (Exp(1) BY COMMUNITY)
# -----------------------------

bootstrap_cluster_weight <- function(d, outcome, mediators, B = BOOT_REPS, seed = SEED_VAL) {
  set.seed(seed)
  
  point <- estimate_interventional_effects(d, outcome, mediators, seed = seed, verbose_specs = TRUE)
  
  cmp <- point$treat_level
  k <- length(cmp)
  
  TE_b  <- matrix(NA_real_, B, k, dimnames = list(NULL, cmp))
  IDE_b <- matrix(NA_real_, B, k, dimnames = list(NULL, cmp))
  IIE_b <- matrix(NA_real_, B, k, dimnames = list(NULL, cmp))
  
  cl <- factor(d$Community)
  cl_lvls <- levels(cl)
  
  # -----------------------------
  # PARALLEL BOOTSTRAP (PSOCK, fault-tolerant)
  # -----------------------------
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
  
  n_cores <- 20L               
  chunk_size <- max(1L, 2L * n_cores)
  max_restarts <- 5L
  
  boot_one <- function(b) {
    tryCatch({
      set.seed(seed + 100000L + b)
      
      w_cl <- stats::rexp(length(cl_lvls), rate = 1)
      names(w_cl) <- cl_lvls
      w_row <- w_cl[as.character(cl)]
      
      est_b <- estimate_interventional_effects(d, outcome, mediators, seed = seed + b, w = w_row)
      est_b <- dplyr::right_join(tibble::tibble(treat_level = cmp), est_b, by = "treat_level")
      
      list(TE = est_b$TE, IDE = est_b$IDE, IIE = est_b$IIE)
    }, error = function(e) {
      list(TE = rep(NA_real_, k), IDE = rep(NA_real_, k), IIE = rep(NA_real_, k), .err = conditionMessage(e))
    })
  }
  
  init_cluster <- function(nc) {
    clust <- parallel::makeCluster(nc, outfile = "")  
    parallel::clusterEvalQ(clust, {
      Sys.setenv(
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        VECLIB_MAXIMUM_THREADS = "1"
      )
      library(dplyr)
      library(tibble)
      library(forcats)
      library(splines)
      invisible(NULL)
    })
    
    parallel::clusterExport(
      clust,
      varlist = c(
        "estimate_interventional_effects",
        "simulate_mediators_residboot",
        "make_mediator_spec",
        "skewness",
        "add_zero_weight_rows_for_factors",
        "align_levels",
        "rhs_with_splines",
        "make_stratified_folds",
        "BASE_VARS",
        "REF_COMMUNITY",
        "KFOLDS",
        "MC_DRAWS",
        "MIN_LEVEL_N"
      ),
      envir = .GlobalEnv
    )
    
    parallel::clusterExport(
      clust,
      varlist = c("d", "outcome", "mediators", "cmp", "k", "cl", "cl_lvls", "seed", "boot_one"),
      envir = environment()
    )
    
    clust
  }
  
  clust <- init_cluster(n_cores)
  cleanup <- function() {
    if (!is.null(clust)) {
      try(parallel::stopCluster(clust), silent = TRUE)
      clust <<- NULL
    }
  }
  on.exit(cleanup(), add = TRUE)
  
  res_list <- vector("list", B)
  pending <- seq_len(B)
  restarts <- 0L
  
  while (length(pending) > 0) {
    idx <- pending[seq_len(min(chunk_size, length(pending)))]
    
    tmp <- tryCatch(
      parallel::parLapplyLB(clust, X = idx, fun = boot_one),
      error = function(e) e
    )
    
    if (inherits(tmp, "error")) {
      restarts <- restarts + 1L
      message("\n[BOOT] worker down: ", conditionMessage(tmp),
              "\n[BOOT] restarting (restart #", restarts, ")...")
      
      cleanup()
      
      if (restarts > max_restarts) {
        warning("\n[BOOT] Demasiados crashes. Relleno los bootstraps restantes con NA para terminar.")
        for (bb in pending) {
          res_list[[bb]] <- list(TE = rep(NA_real_, k), IDE = rep(NA_real_, k), IIE = rep(NA_real_, k))
        }
        break
      }
      
      n_cores <- max(1L, floor(n_cores / 2L))
      chunk_size <- max(1L, 2L * n_cores)
      clust <- init_cluster(n_cores)
      next
    }
    
    for (i in seq_along(idx)) res_list[[idx[i]]] <- tmp[[i]]
    pending <- setdiff(pending, idx)
    cat(sprintf("\rBootstrap %d/%d (cores=%d)", B - length(pending), B, n_cores))
  }
  cat("\n")
  
  for (b in 1:B) {
    TE_b[b, ]  <- res_list[[b]]$TE
    IDE_b[b, ] <- res_list[[b]]$IDE
    IIE_b[b, ] <- res_list[[b]]$IIE
  }
  
  ci <- function(mat) {
    lo <- apply(mat, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
    hi <- apply(mat, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
    se <- apply(mat, 2, stats::sd, na.rm = TRUE)
    list(lo = lo, hi = hi, se = se)
  }
  
  te  <- ci(TE_b)
  ide <- ci(IDE_b)
  iie <- ci(IIE_b)
  
  point %>%
    dplyr::mutate(
      TE_se  = te$se,  TE_lo  = te$lo,  TE_hi  = te$hi,
      IDE_se = ide$se, IDE_lo = ide$lo, IDE_hi = ide$hi,
      IIE_se = iie$se, IIE_lo = iie$lo, IIE_hi = iie$hi,
      ci_method = "Cluster-weight bootstrap (Exp(1) by Community)",
      mc_draws = MC_DRAWS,
      kfolds = KFOLDS
    )
}


# -----------------------------
# 7) RUN ANALYSES
# -----------------------------
cat(">>> STARTING INTERVENTIONAL MEDIATION ANALYSES...\n")
final_res <- list()

for (out_name in names(MEDIATOR_CONFIG)) {
  cat("\n=== Outcome:", out_name, "===\n")
  meds <- MEDIATOR_CONFIG[[out_name]]
  
  res <- bootstrap_cluster_weight(df_clean, out_name, meds, B = BOOT_REPS, seed = SEED_VAL)
  res <- res %>% dplyr::mutate(Outcome = out_name)
  final_res[[out_name]] <- res
  
  readr::write_csv(res, file.path("Causal_Mediation_Output", paste0("Interventional_", out_name, ".csv")))
}

all_data <- dplyr::bind_rows(final_res)
readr::write_csv(all_data, "Causal_Mediation_Output/Full_Interventional_Mediation_Results.csv")
cat("\n>>> Saved results CSVs.\n")

# -----------------------------
# 8) FIGURE
# -----------------------------
cat(">>> GENERATING FIGURE...\n")



library(dplyr)
library(readr)

files <- list.files("Causal_Mediation_Output",
                    pattern = "^Interventional_Outcome_.*\\.csv$",
                    full.names = TRUE)

stopifnot(length(files) > 0)

all_data_rec <- bind_rows(lapply(files, readr::read_csv, show_col_types = FALSE))

out <- file.path("Causal_Mediation_Output",
                 paste0("Full_Interventional_Mediation_Results_RECOVERED_",
                        format(Sys.time(), "%Y%m%d_%H%M%S"),
                        ".csv"))

readr::write_csv(all_data_rec, out)
cat("OK, reconstruido en:\n", out, "\n")




df_viz <- all_data_rec %>%
  dplyr::mutate(
    Outcome_Label = stringr::str_remove(Outcome, "Outcome_"),
    Is_C2 = stringr::str_detect(treat_level, "^C2")
  ) %>%
  dplyr::rename(
    TE_value  = TE,
    IDE_value = IDE,
    IIE_value = IIE
  ) %>%
  tidyr::pivot_longer(
    cols = c(TE_value, TE_lo, TE_hi,
             IDE_value, IDE_lo, IDE_hi,
             IIE_value, IIE_lo, IIE_hi),
    names_to = c("Effect", ".value"),
    names_pattern = "^(TE|IDE|IIE)_(value|lo|hi)$"
  ) %>%
  dplyr::mutate(
    Effect = dplyr::recode(Effect,
                           TE="Total Effect",
                           IDE="Interventional Direct",
                           IIE="Interventional Indirect")
  )

theme_nature <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.5),
      axis.text = ggplot2::element_text(color = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 12),
      legend.position = "bottom",
      panel.grid.major.x = ggplot2::element_line(color="grey92"),
      panel.grid.major.y = ggplot2::element_line(color="grey95")
    )
}

p_forest <- ggplot2::ggplot(df_viz, ggplot2::aes(x = value, y = treat_level, color = Effect)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi),
                          height = 0.28, position = ggplot2::position_dodge(width = 0.65),
                          linewidth = 0.55) +
  ggplot2::geom_point(size = 2.3, position = ggplot2::position_dodge(width = 0.65)) +
  ggplot2::facet_wrap(~Outcome_Label, nrow = 1, scales = "free_x") +
  ggplot2::scale_color_manual(values = c(
    "Total Effect" = "black",
    "Interventional Direct" = "#0072B2",
    "Interventional Indirect" = "#D55E00"
  )) +
  ggplot2::labs(
    title = "Interventional Mediation Decomposition (Cluster C2)",
    subtitle = paste0("RD vs Ref ", REF_COMMUNITY,
                      " | Cross-fit K=", KFOLDS, " | MC draws=", MC_DRAWS,
                      " | Cluster bootstrap B=", BOOT_REPS,
                      " | Residual bootstrap for mediators"),
    x = "Risk Difference (RD) vs reference",
    y = "Community",
    color = "Component"
  ) +
  theme_nature()

ggplot2::ggsave("Causal_Mediation_Output/Figure_Interventional_Mediation.pdf", p_forest, width = 12, height = 6)
ggplot2::ggsave("Causal_Mediation_Output/Figure_Interventional_Mediation.png", p_forest, width = 12, height = 6, dpi = 300)

cat("\n>>> DONE. Outputs saved in 'Causal_Mediation_Output/'.\n")
