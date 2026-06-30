library(readxl)
library(plm)
library(car)
library(glmmTMB)
library(MASS)
library(fixest)
library(splines)
library(sandwich)
library(lmtest)
library(clubSandwich)
library(fwildclusterboot)
library(dqrng)
library(parallel)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(scales)
library(ggrepel)
library(extrafont)
library(dplyr)

options(scipen = 999)
sf_use_s2(FALSE)

OUTDIR   <- "."
DATA_DIR <- "data"
fig_dir  <- "figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
data_file <- file.path(DATA_DIR, "ssa_food_unrest_panel.xlsx")
if (!file.exists(data_file)) stop("Data file not found: data/ssa_food_unrest_panel.xlsx (run 01_build_panel.py first, and run R from the repo root)")

df <- read_excel(data_file)
df <- df %>% filter(iso3 != "SDN") %>% arrange(iso3, Year)
df <- as.data.frame(df)
pdata <- pdata.frame(df, index = c("iso3", "Year"))

gini_mean <- mean(df$gini, na.rm = TRUE)
und_mean  <- mean(df$undernourish_pct, na.rm = TRUE)

m <- df %>%
  arrange(iso3, Year) %>%
  group_by(iso3) %>%
  mutate(
    food_infl_l1    = dplyr::lag(food_infl),
    undernourish_l1 = dplyr::lag(undernourish_pct),
    gov_effect_l1   = dplyr::lag(gov_effect),
    oop_health_l1   = dplyr::lag(oop_health_pct),
    d_gini          = gini - dplyr::lag(gini)
  ) %>%
  ungroup() %>%
  mutate(
    iso3         = factor(iso3),
    yearf        = factor(Year),
    logP         = log(population),
    log_gdp      = log(gdp_pc),
    G_c          = gini - gini_mean,
    food_infl2   = food_infl^2,
    fp_Gc        = food_infl * G_c,
    fp_Gc_l1     = food_infl_l1 * G_c,
    u_rate       = log1p((unrest_n / population) * 1e6),
    log_unrest   = log1p(unrest_n),
    UND_c        = undernourish_pct - und_mean,
    UND_l1_c     = undernourish_l1 - mean(undernourish_l1, na.rm = TRUE),
    Gobs_c       = gini_wb - mean(gini_wb, na.rm = TRUE),
    fp_Gobs      = food_infl * (gini_wb - mean(gini_wb, na.rm = TRUE)),
    food_excess  = food_infl - inflation_pct,
    food_excess2 = (food_infl - inflation_pct)^2,
    fe_Gc        = (food_infl - inflation_pct) * G_c,
    infl_Gc      = inflation_pct * G_c,
    inflation2   = inflation_pct^2,
    fs_z         = as.numeric(scale(food_shock)),
    fs_z2        = as.numeric(scale(food_shock))^2,
    fsz_Gc       = as.numeric(scale(food_shock)) * G_c,
    shock2       = food_shock^2,
    shock_Gc     = food_shock * G_c
  )
m <- as.data.frame(m)

fw_lo <- quantile(m$food_infl, 0.01, na.rm = TRUE)
fw_hi <- quantile(m$food_infl, 0.99, na.rm = TRUE)
m$food_infl_w <- pmin(pmax(m$food_infl, fw_lo), fw_hi)

und_med <- median(m$undernourish_pct, na.rm = TRUE)
m$UND_high <- ifelse(m$undernourish_pct >= und_med, 1L, 0L)
m$und_tercile <- cut(m$undernourish_pct,
                     breaks = quantile(m$undernourish_pct, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                     include.lowest = TRUE, labels = c("low", "mid", "high"))

imp_share    <- tapply(m$food_import_pct, m$iso3, mean, na.rm = TRUE)
m$imp_base   <- as.numeric(imp_share[as.character(m$iso3)])
m$imp_base_c <- m$imp_base - mean(m$imp_base, na.rm = TRUE)
m$Z_bartik   <- m$global_food_infl * (m$imp_base / 100)
m$Z_bartikG  <- m$Z_bartik * m$G_c
m$Zg         <- m$global_food_infl
m$ZgG        <- m$Zg * m$G_c

f_count <- unrest_n ~ food_infl + food_infl2 + G_c + fp_Gc +
  log_gdp + urban_pct + food_import_pct + undernourish_pct +
  inflation_pct + democracy | iso3 + yearf

extract_core <- function(model, terms, mname) {
  s <- summary(model)$coefficients$cond
  cat(sprintf("\n--- %s | N = %d | logLik = %.1f | AIC = %.1f ---\n",
              mname, nobs(model), as.numeric(logLik(model)), AIC(model)))
  for (tm in terms) {
    if (tm %in% rownames(s)) {
      est <- s[tm, "Estimate"]; se <- s[tm, "Std. Error"]; p <- s[tm, "Pr(>|z|)"]
      cat(sprintf("  %-22s b = %8.4f | SE = %7.4f | IRR = %6.3f | p = %.4g\n",
                  tm, est, se, exp(est), p))
    }
  }
}

get_int <- function(mod, term = "food_infl:G_c") {
  s <- tryCatch(summary(mod)$coefficients$cond, error = function(e) NULL)
  if (is.null(s) || !(term %in% rownames(s)))
    return(setNames(c(NA, NA, NA), c("est", "se", "p")))
  setNames(as.numeric(s[term, c("Estimate", "Std. Error", "Pr(>|z|)")]), c("est", "se", "p"))
}

report <- function(lab, mod, term = "food_infl:G_c") {
  v <- get_int(mod, term)
  cat(sprintf("  %-40s b=%9.5f  SE=%8.5f  p=%.4f  %s\n",
              lab, v["est"], v["se"], v["p"],
              ifelse(!is.na(v["p"]) && v["p"] < 0.05, "reject H0", "--")))
}

wald_terms <- function(model, terms) {
  b <- fixef(model)$cond; V <- vcov(model)$cond
  terms <- terms[terms %in% names(b)]
  bs <- b[terms]; Vs <- V[terms, terms, drop = FALSE]
  W <- as.numeric(t(bs) %*% solve(Vs) %*% bs)
  data.frame(terms = paste(terms, collapse = " + "), chisq = W,
             df = length(terms), p = pchisq(W, length(terms), lower.tail = FALSE))
}

marg_fp <- function(model, FP, GC, intname = "food_infl:G_c") {
  b <- fixef(model)$cond; V <- vcov(model)$cond
  out <- list()
  for (fp in FP) for (gc in GC) {
    g <- setNames(rep(0, length(b)), names(b))
    g["food_infl"] <- 1; g["food_infl2"] <- 2 * fp; g[intname] <- gc
    est <- b["food_infl"] + 2 * b["food_infl2"] * fp + b[intname] * gc
    se  <- sqrt(as.numeric(t(g) %*% V %*% g))
    out[[length(out) + 1]] <- data.frame(FP = fp, G_c = gc, slope = as.numeric(est),
                                         SE = se, IRR = exp(as.numeric(est)),
                                         p = 2 * pnorm(abs(est / se), lower.tail = FALSE))
  }
  do.call(rbind, out)
}

ek <- function(mod, term, lab) {
  s <- summary(mod)$coefficients$cond
  data.frame(spec = lab, term = term, N = nobs(mod), AIC = round(AIC(mod), 1),
             est = s[term, "Estimate"], SE = s[term, "Std. Error"], IRR = exp(s[term, "Estimate"]),
             p = s[term, "Pr(>|z|)"], decision = ifelse(s[term, "Pr(>|z|)"] < .05, "Reject H0", "Fail"))
}

extr <- function(mod, nm) {
  s <- summary(mod)$coefficients$cond
  data.frame(model = nm, N = nobs(mod), AIC = round(AIC(mod), 1),
             est = s["food_infl:G_c", "Estimate"], SE = s["food_infl:G_c", "Std. Error"],
             IRR = exp(s["food_infl:G_c", "Estimate"]), p = s["food_infl:G_c", "Pr(>|z|)"])
}

base_fit <- function(dat) {
  glmmTMB(unrest_n ~ food_infl + food_infl2 + G_c + food_infl:G_c +
            log_gdp + urban_pct + food_import_pct + undernourish_pct +
            inflation_pct + democracy + (1 | iso3) + (1 | yearf),
          offset = logP, family = nbinom2, data = dat)
}

cd_test <- function(varname, data) {
  d <- droplevels(as.data.frame(data[!is.na(data[[varname]]), c("iso3", "Year", varname)]))
  pd <- pdata.frame(d, index = c("iso3", "Year"))
  f <- as.formula(paste(varname, "~ 1"))
  r <- tryCatch(pcdtest(f, data = pd, test = "cd"), error = function(e) NULL)
  if (!is.null(r))
    cat(sprintf("  %-16s CD = %8.3f | p = %.4g | %s\n", varname,
                as.numeric(r$statistic), r$p.value,
                ifelse(r$p.value < 0.05, "reject H0 (dependence)", "fail to reject")))
  else
    cat(sprintf("  %-16s CD testi calismadi\n", varname))
}

cips_run <- function(varname, data, drop_iso = character(0), start_year = NULL,
                     diff = FALSE, type = "drift", lags = 1, min_n = 5) {
  d <- as.data.frame(data[!is.na(data[[varname]]) & !(data$iso3 %in% drop_iso),
                          c("iso3", "Year", varname)])
  if (!is.null(start_year)) d <- d[d$Year >= start_year, ]
  d <- d[order(d$iso3, d$Year), ]
  ok <- names(which(table(droplevels(factor(d$iso3))) >= min_n))
  d <- droplevels(d[as.character(d$iso3) %in% ok, ])
  pd <- pdata.frame(d, index = c("iso3", "Year"))
  x <- pd[[varname]]
  if (diff) x <- diff(x)
  r <- tryCatch(cipstest(x, lags = lags, type = type, model = "cmg", truncated = TRUE),
                error = function(e) NULL)
  tag <- paste0(varname, if (diff) " (diff)" else " (levels)")
  if (!is.null(r)) {
    pv <- ifelse(r$p.value <= 0.01, "<= 0.01",
                 ifelse(r$p.value >= 0.10, ">= 0.10", sprintf("= %.3f", r$p.value)))
    cat(sprintf("  %-22s type=%-5s lag=%d | CIPS = %7.3f | p %s | %s\n",
                tag, type, lags, as.numeric(r$statistic), pv,
                ifelse(r$p.value < 0.05, "stationary", "unit root")))
  } else {
    cat(sprintf("  %-22s type=%-5s lag=%d | calismadi\n", tag, type, lags))
  }
}

dh_run <- function(data, yvar, xvar, lags, drop_iso = character(0),
                   start_year = NULL, label = "") {
  d <- data %>% filter(!iso3 %in% drop_iso) %>%
    dplyr::select(iso3, Year, all_of(c(yvar, xvar))) %>%
    filter(!is.na(.data[[yvar]]) & !is.na(.data[[xvar]]))
  if (!is.null(start_year)) d <- d %>% filter(Year >= start_year)
  d <- d %>% arrange(iso3, Year)
  ok <- d %>% group_by(iso3) %>% summarise(n = n(), .groups = "drop") %>%
    filter(n >= (2 * lags + 3))
  d <- d %>% filter(iso3 %in% ok$iso3)
  d <- droplevels(as.data.frame(d))
  pd <- pdata.frame(d, index = c("iso3", "Year"))
  f <- as.formula(paste(yvar, "~", xvar))
  r <- tryCatch(pgrangertest(f, data = pd, order = lags, test = "Ztilde"),
                error = function(e) { cat("   ERR:", conditionMessage(e), "\n"); NULL })
  cat(sprintf("  %-26s lag %d | panel: %d countries | ", label, lags, length(unique(d$iso3))))
  if (!is.null(r)) {
    cat(sprintf("Ztilde = %7.3f | p = %.4g | %s\n",
                as.numeric(r$statistic), r$p.value,
                ifelse(r$p.value < 0.05, "reject H0 (Granger-causes)", "fail to reject")))
  } else cat("\n")
  invisible(r)
}

cat("===== PANEL STRUCTURE =====\n")
cat("Countries:", length(unique(df$iso3)), "| Years:", min(df$Year), "-", max(df$Year), "\n")
cat("Total obs:", nrow(df), "\n")
print(pdim(pdata))

cat("\n===== MISSINGNESS (core vars) =====\n")
core_vars <- c("unrest_n", "food_infl", "gini")
for (v in core_vars) {
  cat(sprintf("%-12s non-missing: %d / %d (%.1f%%)\n",
              v, sum(!is.na(df[[v]])), nrow(df), 100 * mean(!is.na(df[[v]]))))
}

cat("\n===== NA STRUCTURE: food_infl & gini =====\n")
for (v in c("food_infl", "gini")) {
  cat(sprintf("\n### %s\n", v))
  na_by_country <- df %>% group_by(iso3) %>%
    summarise(n = n(), n_na = sum(is.na(.data[[v]])), .groups = "drop") %>%
    filter(n_na > 0) %>% arrange(desc(n_na))
  cat("Countries with any NA (count of NA years):\n")
  if (nrow(na_by_country) == 0) cat("  none\n") else print(as.data.frame(na_by_country), row.names = FALSE)
  na_by_year <- df %>% group_by(Year) %>%
    summarise(n_na = sum(is.na(.data[[v]])), .groups = "drop") %>% filter(n_na > 0)
  cat("\nYears with any NA (count of NA countries):\n")
  print(as.data.frame(na_by_year), row.names = FALSE)
  full_countries <- df %>% group_by(iso3) %>%
    summarise(n_na = sum(is.na(.data[[v]])), .groups = "drop") %>% filter(n_na == 0) %>% nrow()
  cat(sprintf("\nFully-observed countries for %s: %d / 47\n", v, full_countries))
}

cat("\n===== DESCRIPTIVES =====\n")
dvars <- c("unrest_n", "food_infl", "gini", "undernourish_pct", "gdp_pc", "log_gdp",
           "urban_pct", "food_import_pct", "inflation_pct", "democracy", "gov_effect", "oop_health_pct")
desc <- do.call(rbind, lapply(dvars, function(v) {
  x <- m[[v]]
  data.frame(var = v, N = sum(!is.na(x)),
             mean = round(mean(x, na.rm = TRUE), 2), sd = round(sd(x, na.rm = TRUE), 2),
             med = round(median(x, na.rm = TRUE), 2),
             q1 = round(quantile(x, .25, na.rm = TRUE), 2), q3 = round(quantile(x, .75, na.rm = TRUE), 2),
             min = round(min(x, na.rm = TRUE), 2), max = round(max(x, na.rm = TRUE), 2))
}))
print(desc, row.names = FALSE)

cat("\n===== CROSS-SECTIONAL DEPENDENCE (Pesaran CD; H0: independence) =====\n")
for (v in core_vars) cd_test(v, m)

cat("\n===== CIPS PANEL UNIT-ROOT (Pesaran, 2nd generation; H0: unit root) =====\n")
cat("--- core series, levels (drift + trend) ---\n")
cips_run("unrest_n",  m, lags = 1, type = "drift")
cips_run("unrest_n",  m, lags = 1, type = "trend")
cips_run("food_infl", m, drop_iso = c("ERI", "CAF"), start_year = 2001, lags = 1, type = "drift")
cips_run("food_infl", m, drop_iso = c("ERI", "CAF"), start_year = 2001, lags = 1, type = "trend")
cips_run("gini",      m, drop_iso = c("ERI", "SOM"), lags = 1, type = "drift")
cips_run("gini",      m, drop_iso = c("ERI", "SOM"), lags = 1, type = "trend")

cat("--- first differences ---\n")
cips_run("food_infl", m, drop_iso = c("ERI", "CAF"), start_year = 2001, diff = TRUE, lags = 1, type = "drift", min_n = 6)
cips_run("gini",      m, drop_iso = c("ERI", "SOM"), diff = TRUE, lags = 1, type = "drift", min_n = 6)

cat("--- log event rate u_rate (exposure-consistent outcome) ---\n")
cips_run("u_rate", m, lags = 1, type = "drift")
cips_run("u_rate", m, lags = 2, type = "drift")
cips_run("u_rate", m, lags = 1, type = "trend")
cips_run("u_rate", m, lags = 2, type = "trend")

cat("--- table version (drift, 2 lags) ---\n")
cips_run("unrest_n",  m, lags = 2, type = "drift")
cips_run("food_infl", m, drop_iso = c("ERI", "CAF"), start_year = 2001, lags = 2, type = "drift")
cips_run("gini",      m, drop_iso = c("ERI", "SOM"), lags = 2, type = "drift")
cips_run("gini",      m, drop_iso = c("ERI", "SOM"), diff = TRUE, lags = 2, type = "drift", min_n = 6)

cat("\n===== GINI VARIANCE DIAGNOSTICS =====\n")
gini_diag <- m %>% filter(iso3 != "ERI") %>%
  group_by(iso3) %>%
  summarise(n_obs = sum(!is.na(gini)), sd_gini = sd(gini, na.rm = TRUE), .groups = "drop") %>%
  arrange(sd_gini)
cat("Lowest-variance 15 countries (gini):\n")
print(as.data.frame(head(gini_diag, 15)), row.names = FALSE)
cat(sprintf("\nConstant (sd=0 or NA) countries: %d\n", sum(gini_diag$sd_gini == 0 | is.na(gini_diag$sd_gini))))
cat(sprintf("Countries with sd < 0.5: %d\n", sum(gini_diag$sd_gini < 0.5, na.rm = TRUE)))
cat(sprintf("gini_wb (raw, no interpolation) coverage: %d / %d\n", sum(!is.na(m$gini_wb)), nrow(m)))

cat("\n===== SWIID INTERPOLATION CHECK =====\n")
cat("gini variables present:", paste(grep("gini", names(m), value = TRUE, ignore.case = TRUE), collapse = ", "), "\n")
for (v in grep("gini", names(m), value = TRUE, ignore.case = TRUE)) {
  cat(sprintf("  %-14s NA=%4d  range=[%.1f, %.1f]\n", v, sum(is.na(m[[v]])),
              min(m[[v]], na.rm = TRUE), max(m[[v]], na.rm = TRUE)))
}
chk <- m[order(m$iso3, m$Year), c("iso3", "Year", "gini")]
chk$same_as_prev <- c(FALSE, chk$gini[-1] == chk$gini[-nrow(chk)] &
                        as.character(chk$iso3)[-1] == as.character(chk$iso3)[-nrow(chk)])
cat("Share of consecutive identical gini (interpolation signal):",
    round(mean(chk$same_as_prev, na.rm = TRUE), 3), "\n")
cat("food_infl  SD:", round(sd(m$food_infl, na.rm = TRUE), 3),
    "| food_shock SD:", round(sd(m$food_shock, na.rm = TRUE), 3),
    "| ratio ~", round(sd(m$food_infl, na.rm = TRUE) / sd(m$food_shock, na.rm = TRUE), 0), "\n")

cat("\n===== DUMITRESCU-HURLIN PANEL GRANGER NON-CAUSALITY (Ztilde) =====\n")
cat("H0: x does NOT homogeneously Granger-cause y\n")
cat("--- food_infl <-> unrest (levels, both I(0)) ---\n")
for (L in 1:3) {
  dh_run(m, "unrest_n", "food_infl", L, drop_iso = c("ERI", "CAF"), start_year = 2001,
         label = sprintf("food_infl =/=> unrest (L=%d)", L))
  dh_run(m, "food_infl", "unrest_n", L, drop_iso = c("ERI", "CAF"), start_year = 2001,
         label = sprintf("unrest =/=> food_infl (L=%d)", L))
}
cat("--- d_gini -> unrest (gini I(1), robustness direction) ---\n")
for (L in 1:2) {
  dh_run(m, "unrest_n", "d_gini", L, drop_iso = c("ERI", "SOM"),
         label = sprintf("d_gini =/=> unrest (L=%d)", L))
}

cat("\n===== MULTICOLLINEARITY (VIF) =====\n")
vif_data <- m %>%
  dplyr::select(food_infl, gini, gdp_pc, urban_pct, food_import_pct,
                undernourish_pct, inflation_pct, democracy) %>% na.omit()
cat("VIF sample size:", nrow(vif_data), "\n")
print(round(vif(lm(food_infl ~ gini + gdp_pc + urban_pct + food_import_pct +
                     undernourish_pct + inflation_pct + democracy, data = vif_data)), 3))
cat("\n--- with gov_effect added ---\n")
vif_data2 <- m %>%
  dplyr::select(food_infl, gini, gdp_pc, urban_pct, food_import_pct,
                undernourish_pct, inflation_pct, democracy, gov_effect) %>% na.omit()
cat("sample size:", nrow(vif_data2), "\n")
print(round(vif(lm(food_infl ~ gini + gdp_pc + urban_pct + food_import_pct +
                     undernourish_pct + inflation_pct + democracy + gov_effect, data = vif_data2)), 3))

cat("\n--- food_infl & inflation_pct together (collinearity) ---\n")
cat("Pearson r:", round(cor(m$food_infl, m$inflation_pct, use = "complete.obs"), 3), "\n")
lm_chk <- lm(unrest_n ~ food_infl + inflation_pct + G_c + log_gdp + urban_pct +
               food_import_pct + undernourish_pct + democracy, data = m)
print(round(vif(lm_chk), 2))

cat("\n===== FOOD INFLATION DISTRIBUTION =====\n")
print(summary(m$food_infl))
cat("quantiles:\n")
print(round(quantile(m$food_infl, c(0, .01, .05, .25, .5, .75, .95, .99, 1), na.rm = TRUE), 2))
cat("obs with food_infl > 50 :", sum(m$food_infl > 50, na.rm = TRUE), "\n")
cat("obs with food_infl > 100:", sum(m$food_infl > 100, na.rm = TRUE), "\n")
cat("obs with food_infl > 200:", sum(m$food_infl > 200, na.rm = TRUE), "\n")

cat("\n===== FOOD INFLATION VARIANCE DECOMPOSITION =====\n")
av_y <- anova(lm(food_infl ~ factor(Year), data = m))
av_c <- anova(lm(food_infl ~ factor(iso3), data = m))
cat(sprintf("  share explained by YEAR    : %.3f\n", av_y[1, "Sum Sq"] / sum(av_y[, "Sum Sq"])))
cat(sprintf("  share explained by COUNTRY : %.3f\n", av_c[1, "Sum Sq"] / sum(av_c[, "Sum Sq"])))
cat(sprintf("  food_infl winsor range: [%.1f, %.1f] (was [%.1f, %.1f])\n",
            min(m$food_infl_w, na.rm = TRUE), max(m$food_infl_w, na.rm = TRUE),
            min(m$food_infl, na.rm = TRUE),   max(m$food_infl, na.rm = TRUE)))

cat("\n===== HYPERINFLATION OBSERVATIONS (food_infl > 50) =====\n")
hyper <- m %>% filter(food_infl > 50) %>%
  dplyr::select(iso3, Year, food_infl, gini, unrest_n) %>% arrange(desc(food_infl))
print(as.data.frame(hyper), row.names = FALSE)
cat("\ncountries:\n")
print(table(droplevels(hyper$iso3)))
cat(sprintf("\n  mean gini at food_infl>50 : %.1f\n", mean(hyper$gini, na.rm = TRUE)))
cat(sprintf("  mean gini overall         : %.1f\n", mean(m$gini, na.rm = TRUE)))

cat("\n===== MODEL-2 VARIANTS: interaction food_infl:G_c =====\n")
m2a <- glmmTMB(unrest_n ~ food_infl * G_c + log_gdp + urban_pct + food_import_pct +
                 undernourish_pct + inflation_pct + democracy +
                 iso3 + offset(logP), family = nbinom2, data = m)
m2b <- glmmTMB(unrest_n ~ food_infl_w * G_c + log_gdp + urban_pct + food_import_pct +
                 undernourish_pct + inflation_pct + democracy +
                 iso3 + yearf + offset(logP), family = nbinom2, data = m)
m2c <- glmmTMB(unrest_n ~ food_infl_w * G_c + log_gdp + urban_pct + food_import_pct +
                 undernourish_pct + inflation_pct + democracy +
                 iso3 + offset(logP), family = nbinom2, data = m)
sa <- summary(m2a)$coefficients$cond
sb <- summary(m2b)$coefficients$cond
sc <- summary(m2c)$coefficients$cond
cat(sprintf("(a) NO year FE      : int b = %.4f | p = %.4g\n", sa["food_infl:G_c", "Estimate"], sa["food_infl:G_c", "Pr(>|z|)"]))
cat(sprintf("(b) winsor + yearFE : int b = %.4f | p = %.4g\n", sb["food_infl_w:G_c", "Estimate"], sb["food_infl_w:G_c", "Pr(>|z|)"]))
cat(sprintf("(c) winsor, NO yrFE : int b = %.4f | p = %.4g\n", sc["food_infl_w:G_c", "Estimate"], sc["food_infl_w:G_c", "Pr(>|z|)"]))

cat("\n===== NON-LINEARITY: linear vs quadratic price response =====\n")
m2 <- glmmTMB(unrest_n ~ food_infl * G_c + log_gdp + urban_pct + food_import_pct +
                undernourish_pct + inflation_pct + democracy +
                iso3 + yearf + offset(logP), family = nbinom2, data = m)
m2q <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 +
                 log_gdp + urban_pct + food_import_pct +
                 undernourish_pct + inflation_pct + democracy +
                 iso3 + yearf + offset(logP), family = nbinom2, data = m)
sq <- summary(m2q)$coefficients$cond
cat("quadratic model key terms:\n")
for (tm in c("food_infl", "food_infl2", "G_c", "food_infl:G_c")) {
  if (tm %in% rownames(sq)) cat(sprintf("  %-16s b = %10.5f | p = %.4g\n", tm, sq[tm, "Estimate"], sq[tm, "Pr(>|z|)"]))
}
b1 <- sq["food_infl", "Estimate"]; b2 <- sq["food_infl2", "Estimate"]
cat(sprintf("  inverted-U turning point at food_infl = %.1f%%\n", -b1 / (2 * b2)))
cat(sprintf("\n  Model 2 linear    : AIC = %.1f\n", AIC(m2)))
cat(sprintf("  Model 2 quadratic : AIC = %.1f\n", AIC(m2q)))
cat(sprintf("  improvement       : %.1f\n", AIC(m2) - AIC(m2q)))
cat("LR test linear vs quadratic:\n")
print(anova(m2, m2q))

cat("\n===== QUADRATIC ROBUSTNESS TO HYPERINFLATION =====\n")
m_norm <- m %>% filter(food_infl < 50)
m2n <- glmmTMB(unrest_n ~ food_infl * G_c + log_gdp + urban_pct + food_import_pct +
                 undernourish_pct + inflation_pct + democracy +
                 iso3 + yearf + offset(logP), family = nbinom2, data = m_norm)
sn <- summary(m2n)$coefficients$cond
cat(sprintf("Exclude food_infl>50 (N=%d): food_infl b = %.4f | p = %.4g ; food_infl:G_c b = %.4f | p = %.4g\n",
            nobs(m2n), sn["food_infl", "Estimate"], sn["food_infl", "Pr(>|z|)"],
            sn["food_infl:G_c", "Estimate"], sn["food_infl:G_c", "Pr(>|z|)"]))

m_no601 <- m %>% filter(food_infl < 200)
m2_no <- glmmTMB(unrest_n ~ food_infl * G_c + log_gdp + urban_pct + food_import_pct +
                   undernourish_pct + inflation_pct + democracy +
                   iso3 + yearf + offset(logP), family = nbinom2, data = m_no601)
sno <- summary(m2_no)$coefficients$cond
cat(sprintf("Exclude food_infl>200 (N=%d): food_infl:G_c b = %.4f | p = %.4g\n",
            nobs(m2_no), sno["food_infl:G_c", "Estimate"], sno["food_infl:G_c", "Pr(>|z|)"]))

m_noZS <- m %>% filter(!iso3 %in% c("ZWE", "SSD"))
m2q_noZS <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 +
                      log_gdp + urban_pct + food_import_pct +
                      undernourish_pct + inflation_pct + democracy +
                      iso3 + yearf + offset(logP), family = nbinom2, data = m_noZS)
sq2 <- summary(m2q_noZS)$coefficients$cond
cat(sprintf("Quadratic excl ZWE+SSD (N=%d):\n", nobs(m2q_noZS)))
for (tm in c("food_infl", "food_infl2", "G_c", "food_infl:G_c")) {
  if (tm %in% rownames(sq2)) cat(sprintf("  %-16s b = %11.6f | p = %.4g\n", tm, sq2[tm, "Estimate"], sq2[tm, "Pr(>|z|)"]))
}

cat("\n===== CENTRAL FE-NB MODELS (glmmTMB, nbinom2, dummy country+year FE, log offset) =====\n")
M1 <- glmmTMB(unrest_n ~ food_infl + log_gdp + urban_pct + food_import_pct +
                undernourish_pct + inflation_pct + democracy +
                iso3 + yearf + offset(logP), family = nbinom2, data = m)
M2 <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 +
                log_gdp + urban_pct + food_import_pct +
                undernourish_pct + inflation_pct + democracy +
                iso3 + yearf + offset(logP), family = nbinom2, data = m)
M3 <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 + food_infl_l1 * G_c +
                log_gdp + urban_pct + food_import_pct +
                undernourish_pct + inflation_pct + democracy +
                iso3 + yearf + offset(logP), family = nbinom2, data = m)

extract_core(M1, c("food_infl"), "MODEL 1 (linear price)")
extract_core(M2, c("food_infl", "food_infl2", "G_c", "food_infl:G_c"),
             "MODEL 2 (quadratic + interaction) *** CENTRAL ***")
extract_core(M3, c("food_infl", "food_infl2", "food_infl_l1", "G_c",
                   "food_infl:G_c", "food_infl_l1:G_c"), "MODEL 3 (+ lag)")
cat("\n=== dispersion (theta) ===\n")
cat(sprintf("  M1 theta = %.3f | M2 theta = %.3f | M3 theta = %.3f\n",
            sigma(M1), sigma(M2), sigma(M3)))
cat("\n=== full M2 coefficient table (non-FE terms) ===\n")
s2 <- summary(M2)$coefficients$cond
print(round(s2[!grepl("^iso3|^yearf|Intercept", rownames(s2)), ], 5))

cat("\n===== DYNAMICS: lagged price interaction (Model 3) =====\n")
s3 <- summary(M3)$coefficients$cond
lag_int <- if ("food_infl_l1:G_c" %in% rownames(s3)) "food_infl_l1:G_c" else "G_c:food_infl_l1"
for (tm in c("food_infl_l1", lag_int)) {
  if (tm %in% rownames(s3))
    cat(sprintf("  %-18s b = %10.5f | SE = %8.5f | p = %.4g\n",
                tm, s3[tm, "Estimate"], s3[tm, "Std. Error"], s3[tm, "Pr(>|z|)"]))
}
cat("\nH3 joint Wald (lagged price + lagged price x Gini = 0):\n")
print(wald_terms(M3, c("food_infl_l1", lag_int)))

cat("\n===== SAME-SAMPLE COMPARISON M2 vs M3 =====\n")
m3_data <- m %>% filter(!is.na(food_infl_l1))
M2_ss <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 +
                   log_gdp + urban_pct + food_import_pct +
                   undernourish_pct + inflation_pct + democracy +
                   iso3 + yearf + offset(logP), family = nbinom2, data = m3_data)
M3_ss <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 + food_infl_l1 * G_c +
                   log_gdp + urban_pct + food_import_pct +
                   undernourish_pct + inflation_pct + democracy +
                   iso3 + yearf + offset(logP), family = nbinom2, data = m3_data)
cat(sprintf("  same-sample N = %d\n", nobs(M2_ss)))
cat(sprintf("  M2 AIC = %.1f | M3 AIC = %.1f | delta = %.2f\n",
            AIC(M2_ss), AIC(M3_ss), AIC(M2_ss) - AIC(M3_ss)))
cat("  LR test M2 vs M3 (same sample):\n")
print(anova(M2_ss, M3_ss))
sc2 <- summary(M2_ss)$coefficients$cond
cat(sprintf("  contemporaneous food_infl:G_c (same sample) b = %.5f | p = %.4g\n",
            sc2["food_infl:G_c", "Estimate"], sc2["food_infl:G_c", "Pr(>|z|)"]))

cat("\n===== MARGINAL PRICE SLOPE BY GINI x INFLATION (Model 2) =====\n")
g_q  <- quantile(m$gini, c(.10, .25, .50, .75, .90), na.rm = TRUE)
fp_q <- quantile(m$food_infl, c(.25, .50, .75), na.rm = TRUE)
GCq  <- as.numeric(g_q) - gini_mean
ME_gini <- marg_fp(M2, FP = as.numeric(fp_q), GC = GCq)
ME_gini$gini <- gini_mean + ME_gini$G_c
print(round(ME_gini[, c("gini", "FP", "slope", "SE", "IRR", "p")], 5), row.names = FALSE)

cat("\n===== GINI THRESHOLD: where the price slope turns positive =====\n")
b <- fixef(M2)$cond
fp_med <- median(m$food_infl, na.rm = TRUE)
slope_at <- function(gc) b["food_infl"] + 2 * b["food_infl2"] * fp_med + b["food_infl:G_c"] * gc
gc_grid <- seq(min(m$G_c, na.rm = TRUE), max(m$G_c, na.rm = TRUE), length.out = 400)
thr_gc <- gc_grid[which.min(abs(sapply(gc_grid, slope_at)))]
cat(sprintf("  at median food_infl (%.1f%%), slope crosses zero at Gini = %.1f (G_c = %.2f)\n",
            fp_med, gini_mean + thr_gc, thr_gc))
cat(sprintf("  observed Gini range [%.1f, %.1f]; share of obs above threshold = %.1f%%\n",
            min(m$gini, na.rm = TRUE), max(m$gini, na.rm = TRUE),
            100 * mean(m$gini > (gini_mean + thr_gc), na.rm = TRUE)))

cat("\n===== ESTIMATOR CHECKS =====\n")
cat("--- Poisson vs NB (overdispersion) ---\n")
P_M2 <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 +
                  log_gdp + urban_pct + food_import_pct +
                  undernourish_pct + inflation_pct + democracy +
                  iso3 + yearf + offset(logP), family = poisson, data = m)
rp <- residuals(P_M2, type = "pearson")
cat(sprintf("  Poisson AIC = %.1f | NB AIC = %.1f | NB theta = %.3f\n", AIC(P_M2), AIC(M2), sigma(M2)))
cat(sprintf("  Poisson Pearson dispersion = %.2f (>>1 => NB justified)\n", sum(rp^2) / df.residual(P_M2)))
sp <- summary(P_M2)$coefficients$cond
cat(sprintf("  Poisson food_infl:G_c IRR = %.4f | p = %.4g\n",
            exp(sp["food_infl:G_c", "Estimate"]), sp["food_infl:G_c", "Pr(>|z|)"]))

cat("\n--- fixest within-FE estimators (country + year) ---\n")
d_pois <- droplevels(m[complete.cases(m[, c("unrest_n", "food_infl", "food_infl2", "G_c", "fp_Gc",
  "log_gdp", "urban_pct", "food_import_pct", "undernourish_pct", "inflation_pct",
  "democracy", "logP")]), ])
M2_fe   <- fenegbin(f_count, data = m, offset = ~logP)
M2_ppml <- fepois(f_count, data = d_pois, offset = ~logP)
cat("  fixest FE-NB (fenegbin), interaction fp_Gc:\n")
print(round(M2_fe$coeftable["fp_Gc", , drop = FALSE], 5))
cat("  fixest PPML (fepois), interaction fp_Gc:\n")
print(round(M2_ppml$coeftable["fp_Gc", , drop = FALSE], 5))
cat("\n  3-way comparison of interaction (food price x Gini):\n")
cat(sprintf("    glmmTMB NB  : b = %.5f\n", fixef(M2)$cond["food_infl:G_c"]))
cat(sprintf("    fixest NB   : b = %.5f\n", M2_fe$coeftable["fp_Gc", "Estimate"]))
cat(sprintf("    fixest PPML : b = %.5f\n", M2_ppml$coeftable["fp_Gc", "Estimate"]))

cat("\n===== INFLUENCE / NON-LINEARITY TABLE (interaction across specs) =====\n")
infl_tab <- rbind(
  ek(M2,       "food_infl:G_c",   "Central M2 (quadratic)"),
  ek(m2a,      "food_infl:G_c",   "No year FE"),
  ek(m2b,      "food_infl_w:G_c", "Winsorized + year FE"),
  ek(m2c,      "food_infl_w:G_c", "Winsorized, no year FE"),
  ek(m2n,      "food_infl:G_c",   "Exclude food_infl > 50"),
  ek(m2_no,    "food_infl:G_c",   "Exclude food_infl > 200"),
  ek(m2q_noZS, "food_infl:G_c",   "Exclude ZWE + SSD")
)
rownames(infl_tab) <- NULL
print(infl_tab, row.names = FALSE)

cat("\n===== RQ1: NUTRITIONAL VULNERABILITY x FOOD PRICE =====\n")
cat("coefficient of interest: food price x undernourishment (eta2)\n")
ctrl_noUND <- "log_gdp + urban_pct + food_import_pct + inflation_pct + democracy"

M_RQ1 <- glmmTMB(as.formula(paste(
  "unrest_n ~ food_infl * G_c + food_infl2 + food_infl * UND_l1_c +",
  ctrl_noUND, "+ iso3 + yearf + offset(logP)")),
  family = nbinom2, data = m)
r1 <- summary(M_RQ1)$coefficients$cond
cat(sprintf("(lagged, continuous) eta2 = food_infl:UND_l1_c  b = %.5f | SE = %.5f | p = %.4g\n",
            r1["food_infl:UND_l1_c", "Estimate"], r1["food_infl:UND_l1_c", "Std. Error"],
            r1["food_infl:UND_l1_c", "Pr(>|z|)"]))

M_RQ1b <- glmmTMB(as.formula(paste(
  "unrest_n ~ food_infl * G_c + food_infl2 + food_infl * UND_c +",
  ctrl_noUND, "+ iso3 + yearf + offset(logP)")),
  family = nbinom2, data = m)
r1b <- summary(M_RQ1b)$coefficients$cond
cat(sprintf("(contemporaneous)    eta2 = food_infl:UND_c     b = %.5f | p = %.4g\n",
            r1b["food_infl:UND_c", "Estimate"], r1b["food_infl:UND_c", "Pr(>|z|)"]))

M_RQ1c <- glmmTMB(as.formula(paste(
  "unrest_n ~ food_infl * G_c * UND_c + food_infl2 +",
  ctrl_noUND, "+ iso3 + yearf + offset(logP)")),
  family = nbinom2, data = m)
r1c <- summary(M_RQ1c)$coefficients$cond
tt <- if ("food_infl:G_c:UND_c" %in% rownames(r1c)) "food_infl:G_c:UND_c" else grep("UND_c", grep(":G_c", rownames(r1c), value = TRUE), value = TRUE)[1]
cat(sprintf("(triple interaction) %-22s b = %.5f | p = %.4g\n", tt, r1c[tt, "Estimate"], r1c[tt, "Pr(>|z|)"]))

M_RQ1d <- glmmTMB(as.formula(paste(
  "unrest_n ~ food_infl * G_c + food_infl2 + food_infl * UND_high +",
  ctrl_noUND, "+ iso3 + yearf + offset(logP)")),
  family = nbinom2, data = m)
r1d <- summary(M_RQ1d)$coefficients$cond
cat(sprintf("(high/low split)     food_infl:UND_high       b = %.5f | p = %.4g\n",
            r1d["food_infl:UND_high", "Estimate"], r1d["food_infl:UND_high", "Pr(>|z|)"]))

cat("\n--- tercile & spline forms (coefficients) ---\n")
M_RQ1e <- glmmTMB(as.formula(paste(
  "unrest_n ~ food_infl * G_c + food_infl2 + food_infl * und_tercile +",
  ctrl_noUND, "+ iso3 + yearf + offset(logP)")),
  family = nbinom2, data = m)
r1e <- summary(M_RQ1e)$coefficients$cond
print(round(r1e[grep("food_infl:und_tercile", rownames(r1e)), , drop = FALSE], 5))
M_RQ1f <- glmmTMB(as.formula(paste(
  "unrest_n ~ food_infl * G_c + food_infl2 + food_infl * ns(undernourish_pct, 3) +",
  ctrl_noUND, "+ iso3 + yearf + offset(logP)")),
  family = nbinom2, data = m)
r1f <- summary(M_RQ1f)$coefficients$cond
cat("central food_infl:G_c still present in RQ1e:",
    round(r1e["food_infl:G_c", "Estimate"], 5), "| p =", round(r1e["food_infl:G_c", "Pr(>|z|)"], 4), "\n")

cat("\n--- LR tests (complete-case, dummy FE): does food x UND improve fit? ---\n")
dd <- droplevels(m[complete.cases(m[, c("unrest_n", "food_infl", "food_infl2", "G_c",
  "undernourish_pct", "log_gdp", "urban_pct", "food_import_pct", "inflation_pct",
  "democracy", "logP", "und_tercile")]), ])
M_base_t <- glmmTMB(as.formula(paste(
  "unrest_n ~ food_infl * G_c + food_infl2 + und_tercile +",
  ctrl_noUND, "+ iso3 + yearf + offset(logP)")),
  family = nbinom2, data = dd)
M_terc <- glmmTMB(as.formula(paste(
  "unrest_n ~ food_infl * G_c + food_infl2 + food_infl * und_tercile +",
  ctrl_noUND, "+ iso3 + yearf + offset(logP)")),
  family = nbinom2, data = dd)
cat("Tercile interaction LR test:\n")
print(anova(M_base_t, M_terc))
M_base_s <- glmmTMB(as.formula(paste(
  "unrest_n ~ food_infl * G_c + food_infl2 + ns(undernourish_pct, 3) +",
  ctrl_noUND, "+ iso3 + yearf + offset(logP)")),
  family = nbinom2, data = dd)
M_spl <- glmmTMB(as.formula(paste(
  "unrest_n ~ food_infl * G_c + food_infl2 + food_infl * ns(undernourish_pct, 3) +",
  ctrl_noUND, "+ iso3 + yearf + offset(logP)")),
  family = nbinom2, data = dd)
cat("Spline interaction LR test:\n")
print(anova(M_base_s, M_spl))

cat("\n--- marginal price slope by undernourishment (from lagged RQ1 model) ---\n")
und_mean_l1 <- mean(m$undernourish_l1, na.rm = TRUE)
und_q <- quantile(m$undernourish_l1, c(.10, .25, .50, .75, .90), na.rm = TRUE)
ME_und <- marg_fp(M_RQ1, FP = median(m$food_infl, na.rm = TRUE),
                  GC = as.numeric(und_q) - und_mean_l1, intname = "food_infl:UND_l1_c")
ME_und$undernourish <- und_mean_l1 + ME_und$G_c
print(round(ME_und[, c("undernourish", "FP", "slope", "SE", "IRR", "p")], 5), row.names = FALSE)

cat("\n===== UNDERNOURISHMENT MAIN EFFECT ACROSS M1 / M2 / M3 =====\n")
for (mm in list(c("M1", "M1"), c("M2", "M2"), c("M3", "M3"))) {
  s <- summary(get(mm[1]))$coefficients$cond
  if ("undernourish_pct" %in% rownames(s))
    cat(sprintf("  %-4s undernourish_pct b = %.5f | SE = %.5f | p = %.4g\n",
                mm[2], s["undernourish_pct", "Estimate"], s["undernourish_pct", "Std. Error"],
                s["undernourish_pct", "Pr(>|z|)"]))
}

cat("\n===== SENSITIVITY: GOVERNANCE & HEALTH FINANCING =====\n")
mk_sens <- function(extra) {
  f <- as.formula(paste("unrest_n ~ food_infl * G_c + food_infl2 +",
    "log_gdp + urban_pct + food_import_pct + undernourish_pct + inflation_pct + democracy +",
    extra, "+ iso3 + yearf + offset(logP)"))
  glmmTMB(f, family = nbinom2, data = m)
}
M_gov     <- mk_sens("gov_effect_l1")
M_oop     <- mk_sens("oop_health_l1")
M_gov_oop <- mk_sens("gov_effect_l1 + oop_health_l1")
gov_tab <- rbind(extr(M2, "Central M2"),
                 extr(M_gov, "+ gov_effect (l1)"),
                 extr(M_oop, "+ oop_health (l1)"),
                 extr(M_gov_oop, "+ both"))
rownames(gov_tab) <- NULL
print(gov_tab, row.names = FALSE)
cat("\nVIF with out-of-pocket health added:\n")
vif_oop <- m %>% dplyr::select(food_infl, gini, gdp_pc, urban_pct, food_import_pct,
                               undernourish_pct, inflation_pct, democracy, oop_health_l1) %>% na.omit()
print(round(vif(lm(food_infl ~ gini + gdp_pc + urban_pct + food_import_pct +
                     undernourish_pct + inflation_pct + democracy + oop_health_l1, data = vif_oop)), 2))

cat("\n===== ENDOGENEITY: INSTRUMENTAL VARIABLES / CONTROL FUNCTION =====\n")
ctrl_iv <- "log_gdp + urban_pct + food_import_pct + undernourish_pct + inflation_pct + democracy"
has_temp <- "temp_anomaly" %in% names(m)

base_imp <- m %>% filter(Year <= 2004) %>% group_by(iso3) %>%
  summarise(import_base = mean(food_import_pct, na.rm = TRUE), .groups = "drop")
m <- m %>% left_join(base_imp, by = "iso3")
m$Z1  <- m$global_food_infl * m$import_base
m$Z1G <- m$Z1 * m$G_c

cat("--- (a) shift-share / Bartik instruments (relevance check) ---\n")
fs_bartik <- lm(as.formula(paste("food_infl ~ Z_bartik + Z_bartikG +", ctrl_iv,
                                 "+ factor(iso3) + factor(yearf)")), data = m)
F_bartik <- tryCatch(linearHypothesis(fs_bartik, c("Z_bartik = 0", "Z_bartikG = 0"))$F[2],
                     error = function(e) NA)
cat(sprintf("  Bartik first-stage joint F (Z_bartik, Z_bartikG) = %.2f %s\n",
            F_bartik, ifelse(!is.na(F_bartik) && F_bartik < 10, "(WEAK; discarded)", "")))
ss_inst <- if (has_temp) "Z1 + Z1G + temp_anomaly + I(temp_anomaly * G_c)" else "Z1 + Z1G"
IV_ss <- tryCatch(feols(as.formula(paste("u_rate ~", ctrl_iv, "| iso3 + yearf |",
                                         "food_infl + fp_Gc ~", ss_inst)),
                        data = m, cluster = ~iso3), error = function(e) NULL)
if (!is.null(IV_ss)) {
  cat("  shift-share 2SLS (food_infl, fp_Gc instrumented):\n")
  print(round(coef(IV_ss)[c("fit_food_infl", "fit_fp_Gc")], 5))
  cat("  first-stage diagnostics:\n")
  print(tryCatch(fitstat(IV_ss, c("ivf1", "ivwald1"), verbose = FALSE), error = function(e) "n/a"))
}

cat("\n--- (b) direct global-price exposure instrument (principal) ---\n")
fs_dir <- lm(as.formula(paste("food_infl ~ Zg + ZgG +", ctrl_iv,
                              "+ factor(iso3) + factor(yearf)")), data = m)
F_dir <- tryCatch(linearHypothesis(fs_dir, c("Zg = 0", "ZgG = 0"))$F[2],
                  error = function(e) NA)
cat(sprintf("  direct first-stage joint F (Zg, ZgG) = %.2f %s\n",
            F_dir, ifelse(!is.na(F_dir) && F_dir >= 10, "(strong)", "")))

cc_cf <- droplevels(m[complete.cases(m[, c("unrest_n", "food_infl", "G_c", "Zg", "ZgG",
  "log_gdp", "urban_pct", "food_import_pct", "undernourish_pct", "inflation_pct",
  "democracy", "logP")]), ])
fs_cf <- lm(as.formula(paste("food_infl ~ Zg + ZgG +", ctrl_iv,
                             "+ factor(iso3) + factor(yearf)")), data = cc_cf)
cc_cf$nu_hat <- resid(fs_cf)
CF2 <- glmmTMB(as.formula(paste(
  "unrest_n ~ food_infl * G_c + food_infl2 + nu_hat + nu_hat:G_c +",
  ctrl_iv, "+ iso3 + yearf + offset(logP)")),
  family = nbinom2, data = cc_cf)
scf <- summary(CF2)$coefficients$cond
cat("  control-function NB (CF2): structural interaction + endogeneity terms\n")
print(round(scf[intersect(c("food_infl", "food_infl2", "food_infl:G_c", "nu_hat", "nu_hat:G_c"),
                          rownames(scf)), , drop = FALSE], 5))
cat("  (joint significance of nu_hat terms = Wu-Hausman-style endogeneity test)\n")

IV2 <- tryCatch(feols(as.formula(paste("u_rate ~", ctrl_iv, "| iso3 + yearf |",
                                       "food_infl + fp_Gc ~ Zg + ZgG")),
                      data = m, cluster = ~iso3), error = function(e) NULL)
if (!is.null(IV2)) {
  cat("\n  direct-instrument 2SLS (IV2):\n")
  print(summary(IV2))
  cat("\n  first stage:\n")
  print(summary(IV2, stage = 1))
  cat("\n  IV diagnostics (KP first-stage F / Wu-Hausman):\n")
  print(tryCatch(fitstat(IV2, c("ivf1", "ivwald1", "wh"), verbose = FALSE), error = function(e) "n/a"))
}

cat("\n===== ROBUSTNESS BATTERY (interaction food price x inequality) =====\n")
cat("Central reference:\n")
report("Central M2 (glmmTMB NB)", M2, "food_infl:G_c")

cat("\n--- alternative price measures ---\n")
M_detr <- glmmTMB(unrest_n ~ food_shock * G_c + shock2 +
                    log_gdp + urban_pct + food_import_pct + undernourish_pct +
                    inflation_pct + democracy + (1 | iso3) + (1 | yearf),
                  offset = logP, family = nbinom2, data = m)
report("Detrended food shock (raw)", M_detr, "food_shock:G_c")
M_detr_z <- glmmTMB(unrest_n ~ fs_z * G_c + fs_z2 +
                      log_gdp + urban_pct + food_import_pct + undernourish_pct +
                      inflation_pct + democracy + (1 | iso3) + (1 | yearf),
                    offset = logP, family = nbinom2, data = m)
report("Detrended food shock (standardized)", M_detr_z, "fs_z:G_c")
M_excess <- glmmTMB(unrest_n ~ food_excess * G_c + food_excess2 +
                      log_gdp + urban_pct + food_import_pct + undernourish_pct +
                      inflation_pct + democracy + (1 | iso3) + (1 | yearf),
                    offset = logP, family = nbinom2, data = m)
report("Excess food inflation (food - general)", M_excess, "food_excess:G_c")
m_cc <- droplevels(m[complete.cases(m[, c("food_infl", "inflation_pct", "unrest_n", "G_c",
  "log_gdp", "urban_pct", "food_import_pct", "undernourish_pct", "democracy", "logP")]), ])
m_cc$food_orth <- resid(lm(food_infl ~ inflation_pct, data = m_cc))
M_orth <- glmmTMB(unrest_n ~ food_orth * G_c + I(food_orth^2) +
                    log_gdp + urban_pct + food_import_pct + undernourish_pct +
                    inflation_pct + democracy + (1 | iso3) + (1 | yearf),
                  offset = logP, family = nbinom2, data = m_cc)
report("Food inflation orthogonalized to general", M_orth, "food_orth:G_c")

cat("\n--- general inflation channel ---\n")
M_noinf <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 +
                     log_gdp + urban_pct + food_import_pct + undernourish_pct +
                     democracy + (1 | iso3) + (1 | yearf),
                   offset = logP, family = nbinom2, data = m)
report("Drop general inflation control", M_noinf, "food_infl:G_c")
M_gen <- glmmTMB(unrest_n ~ food_infl + food_infl2 + inflation_pct * G_c +
                   log_gdp + urban_pct + food_import_pct + undernourish_pct +
                   democracy + (1 | iso3) + (1 | yearf),
                 offset = logP, family = nbinom2, data = m)
report("General-inflation x Gini (placebo)", M_gen, "inflation_pct:G_c")
M_both <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 + inflation_pct * G_c +
                    log_gdp + urban_pct + food_import_pct + undernourish_pct +
                    democracy + (1 | iso3) + (1 | yearf),
                  offset = logP, family = nbinom2, data = m)
report("Both interactions: food x Gini", M_both, "food_infl:G_c")
report("Both interactions: general x Gini", M_both, "inflation_pct:G_c")

cat("\n--- alternative inequality measure ---\n")
M_obs <- glmmTMB(unrest_n ~ food_infl * Gobs_c + food_infl2 +
                   log_gdp + urban_pct + food_import_pct + undernourish_pct +
                   inflation_pct + democracy + (1 | iso3) + (1 | yearf),
                 offset = logP, family = nbinom2, data = m)
report("Observed World Bank Gini", M_obs, "food_infl:Gobs_c")

cat("\n--- sample / period / estimator variations ---\n")
M_y05 <- base_fit(m %>% filter(Year >= 2005))
report("Start year 2005", M_y05, "food_infl:G_c")
top3 <- names(sort(tapply(m$unrest_n, m$iso3, sum, na.rm = TRUE), decreasing = TRUE))[1:3]
cat("  highest-count countries excluded:", paste(top3, collapse = ", "), "\n")
M_notop <- base_fit(droplevels(m %>% filter(!iso3 %in% top3)))
report("Exclude 3 highest-count countries", M_notop, "food_infl:G_c")
M_pois <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 +
                    log_gdp + urban_pct + food_import_pct + undernourish_pct +
                    inflation_pct + democracy + (1 | iso3) + (1 | yearf),
                  offset = logP, family = poisson, data = m)
report("Poisson (random effects)", M_pois, "food_infl:G_c")
M_loglin <- lm(log_unrest ~ food_infl * G_c + food_infl2 +
                 log_gdp + urban_pct + food_import_pct + undernourish_pct +
                 inflation_pct + democracy + factor(iso3) + factor(yearf), data = m)
sll <- summary(M_loglin)$coefficients
cat(sprintf("  %-40s b=%9.5f  SE=%8.5f  p=%.4f\n", "log(unrest+1) OLS (linear FE)",
            sll["food_infl:G_c", "Estimate"], sll["food_infl:G_c", "Std. Error"],
            sll["food_infl:G_c", "Pr(>|t|)"]))

cat("\n--- Southern-Africa influence ---\n")
M_noZAF <- base_fit(droplevels(m %>% filter(iso3 != "ZAF")))
report("Exclude South Africa (ZAF)", M_noZAF, "food_infl:G_c")
south4 <- c("ZAF", "NAM", "SWZ", "BWA")
M_noS4 <- fenegbin(f_count, data = droplevels(m %>% filter(!iso3 %in% south4)), offset = ~logP)
cat(sprintf("  %-40s b=%9.5f  SE=%8.5f  p=%.4f  (fixest NB, Table-5 row)\n",
            "Exclude Southern-4 (ZAF NAM SWZ BWA)",
            M_noS4$coeftable["fp_Gc", "Estimate"], M_noS4$coeftable["fp_Gc", "Std. Error"],
            M_noS4$coeftable["fp_Gc", "Pr(>|z|)"]))

cat("\n--- leave-one-country-out jackknife (fixest NB interaction fp_Gc) ---\n")
ctrys <- levels(droplevels(m$iso3))
jack <- sapply(ctrys, function(cc) {
  fit <- tryCatch(fenegbin(f_count, data = droplevels(m[m$iso3 != cc, ]), offset = ~logP), error = function(e) NULL)
  if (is.null(fit) || !("fp_Gc" %in% rownames(fit$coeftable))) NA else fit$coeftable["fp_Gc", "Estimate"]
})
full_b <- M2_fe$coeftable["fp_Gc", "Estimate"]
cat(sprintf("  full-sample fp_Gc = %.5f\n", full_b))
cat(sprintf("  jackknife range  = [%.5f, %.5f]\n", min(jack, na.rm = TRUE), max(jack, na.rm = TRUE)))
cat("  most influential (largest drop when removed):\n")
infl_rank <- sort(full_b - jack, decreasing = TRUE)
print(round(head(infl_rank, 6), 5))
cat(sprintf("  sign flips to negative when removing: %s\n",
            paste(names(jack)[which(jack < 0)], collapse = ", ")))

cat("\n===== INFERENCE ON THE INTERACTION (food price x Gini) =====\n")
cat("Term reported throughout: fp_Gc (fixest) / food_infl:G_c (glmmTMB)\n")

cat("\n--- (1) linear two-way FE benchmark: classic / cluster / Driscoll-Kraay ---\n")
d_lin <- droplevels(m[complete.cases(m[, c("u_rate", "food_infl", "food_infl2", "G_c", "fp_Gc",
  "log_gdp", "urban_pct", "food_import_pct", "undernourish_pct", "inflation_pct", "democracy")]), ])
pd_lin <- pdata.frame(d_lin, index = c("iso3", "Year"))
fe_plm <- plm(u_rate ~ food_infl + food_infl2 + G_c + fp_Gc + log_gdp + urban_pct +
                food_import_pct + undernourish_pct + inflation_pct + democracy,
              data = pd_lin, model = "within", effect = "twoways")
se_classic <- coeftest(fe_plm)["fp_Gc", ]
se_clust   <- coeftest(fe_plm, vcov = function(x) vcovHC(x, method = "arellano", type = "HC1", cluster = "group"))["fp_Gc", ]
se_dk      <- coeftest(fe_plm, vcov = function(x) vcovSCC(x, type = "HC1"))["fp_Gc", ]
cat(sprintf("  classic        b=%9.5f  SE=%8.5f  p=%.4f\n", se_classic[1], se_classic[2], se_classic[4]))
cat(sprintf("  cluster(iso3)  b=%9.5f  SE=%8.5f  p=%.4f\n", se_clust[1], se_clust[2], se_clust[4]))
cat(sprintf("  Driscoll-Kraay b=%9.5f  SE=%8.5f  p=%.4f\n", se_dk[1], se_dk[2], se_dk[4]))

cat("\n--- (2) fixest FE-NB (fenegbin): iid / cluster / two-way ---\n")
nb_iid <- summary(M2_fe, vcov = "iid")$coeftable["fp_Gc", ]
nb_cl  <- summary(M2_fe, cluster = ~iso3)$coeftable["fp_Gc", ]
nb_tw  <- summary(M2_fe, cluster = ~iso3 + yearf)$coeftable["fp_Gc", ]
cat(sprintf("  iid            b=%9.5f  SE=%8.5f  p=%.4f\n", nb_iid[1], nb_iid[2], nb_iid[4]))
cat(sprintf("  cluster(iso3)  b=%9.5f  SE=%8.5f  p=%.4f\n", nb_cl[1], nb_cl[2], nb_cl[4]))
cat(sprintf("  two-way        b=%9.5f  SE=%8.5f  p=%.4f\n", nb_tw[1], nb_tw[2], nb_tw[4]))

cat("\n--- (3) PPML (fepois): iid / cluster / two-way ---\n")
pp_iid <- summary(M2_ppml, vcov = "iid")$coeftable["fp_Gc", ]
pp_cl  <- summary(M2_ppml, cluster = ~iso3)$coeftable["fp_Gc", ]
pp_tw  <- summary(M2_ppml, cluster = ~iso3 + yearf)$coeftable["fp_Gc", ]
cat(sprintf("  iid            b=%9.5f  SE=%8.5f  p=%.4f\n", pp_iid[1], pp_iid[2], pp_iid[4]))
cat(sprintf("  cluster(iso3)  b=%9.5f  SE=%8.5f  p=%.4f\n", pp_cl[1], pp_cl[2], pp_cl[4]))
cat(sprintf("  two-way        b=%9.5f  SE=%8.5f  p=%.4f\n", pp_tw[1], pp_tw[2], pp_tw[4]))

cat("\n--- (4) CR2 bias-reduced cluster SEs + Satterthwaite df (clubSandwich) ---\n")
M_pois_glm <- glm(unrest_n ~ food_infl + food_infl2 + G_c + fp_Gc + log_gdp + urban_pct +
                    food_import_pct + undernourish_pct + inflation_pct + democracy +
                    factor(iso3) + factor(yearf) + offset(logP), family = poisson, data = d_pois)
M_glmnb <- glm.nb(unrest_n ~ food_infl + food_infl2 + G_c + fp_Gc + log_gdp + urban_pct +
                    food_import_pct + undernourish_pct + inflation_pct + democracy +
                    factor(iso3) + factor(yearf) + offset(logP), data = d_pois)
ct_pois <- tryCatch(coef_test(M_pois_glm, vcov = "CR2", cluster = d_pois$iso3,
                              test = "Satterthwaite", coefs = "fp_Gc"), error = function(e) NULL)
ct_nb <- tryCatch(coef_test(M_glmnb, vcov = "CR2", cluster = d_pois$iso3,
                            test = "Satterthwaite", coefs = "fp_Gc"), error = function(e) NULL)
gv <- function(x, nm) { for (n in nm) if (n %in% names(x)) return(x[[n]]); NA }
if (!is.null(ct_pois)) cat(sprintf("  PPML  CR2  b=%9.5f  SE=%8.5f  df=%6.2f  p=%.4f\n",
                                   gv(ct_pois, c("beta", "Estimate")), gv(ct_pois, c("SE", "Std. Error")),
                                   gv(ct_pois, c("df_Satt", "df")), gv(ct_pois, c("p_Satt", "p_val"))))
if (!is.null(ct_nb))   cat(sprintf("  NB    CR2  b=%9.5f  SE=%8.5f  df=%6.2f  p=%.4f\n",
                                   gv(ct_nb, c("beta", "Estimate")), gv(ct_nb, c("SE", "Std. Error")),
                                   gv(ct_nb, c("df_Satt", "df")), gv(ct_nb, c("p_Satt", "p_val"))))
cat("  (very small Satterthwaite df signals unreliable large-sample cluster inference)\n")

cat("\n--- (5) linear FE wild cluster bootstrap (fwildclusterboot) ---\n")
set.seed(42)
lm_wcb <- lm(u_rate ~ food_infl + food_infl2 + G_c + fp_Gc + log_gdp + urban_pct +
               food_import_pct + undernourish_pct + inflation_pct + democracy +
               factor(iso3) + factor(yearf), data = d_lin)
wcb_webb <- tryCatch(boottest(lm_wcb, clustid = ~iso3, param = "fp_Gc", B = 9999, type = "webb"),
                     error = function(e) NULL)
wcb_rade <- tryCatch(boottest(lm_wcb, clustid = ~iso3, param = "fp_Gc", B = 9999, type = "rademacher"),
                     error = function(e) NULL)
if (!is.null(wcb_webb)) cat(sprintf("  Webb 6-point   p=%.4f  95%% CI = [%.5f, %.5f]\n",
                                    wcb_webb$p_val, wcb_webb$conf_int[1], wcb_webb$conf_int[2]))
if (!is.null(wcb_rade)) cat(sprintf("  Rademacher     p=%.4f  95%% CI = [%.5f, %.5f]\n",
                                    wcb_rade$p_val, wcb_rade$conf_int[1], wcb_rade$conf_int[2]))

cat("\n--- (6) PPML score wild cluster bootstrap (Webb weights, CGM unrestricted) ---\n")
set.seed(42)
sw <- tryCatch({
  U  <- estfun(M_pois_glm)
  V0 <- vcov(M_pois_glm)
  IF <- U %*% V0
  j  <- which(colnames(U) == "fp_Gc")
  gidx <- split(seq_len(nrow(d_pois)), d_pois$iso3)
  IFg <- t(sapply(gidx, function(ix) colSums(IF[ix, , drop = FALSE])))
  b_j  <- coef(M_pois_glm)[j]
  se_j <- sqrt(sum(IFg[, j]^2))
  t_obs <- b_j / se_j
  webb <- c(-sqrt(1.5), -1, -sqrt(0.5), sqrt(0.5), 1, sqrt(1.5))
  G <- length(gidx); Bsw <- 9999
  tb <- replicate(Bsw, {
    w <- sample(webb, G, replace = TRUE)
    num <- sum(w * IFg[, j]); den <- sqrt(sum((w * IFg[, j])^2))
    num / den
  })
  list(b = b_j, se = se_j, t = t_obs, p = mean(abs(tb) >= abs(t_obs)))
}, error = function(e) NULL)
if (!is.null(sw)) cat(sprintf("  score WCB  b=%9.5f  CR-SE=%8.5f  t=%6.3f  p=%.4f\n",
                              sw$b, sw$se, sw$t, sw$p))

cat("\n--- (7) PPML pairs-cluster bootstrap (parallel, B = 999) ---\n")
B <- 999
ctry_list <- split(seq_len(nrow(d_pois)), d_pois$iso3)
ncl <- length(ctry_list)
boot_one <- function(b) {
  pick <- sample(seq_len(ncl), ncl, replace = TRUE)
  idx  <- unlist(ctry_list[pick], use.names = FALSE)
  db   <- d_pois[idx, ]
  db$iso3 <- factor(rep(seq_len(ncl), times = lengths(ctry_list[pick])))
  fit <- tryCatch(fepois(f_count, data = db, offset = ~logP), error = function(e) NULL)
  if (is.null(fit) || !("fp_Gc" %in% rownames(fit$coeftable))) return(NA_real_)
  fit$coeftable["fp_Gc", "Estimate"]
}
boots <- tryCatch({
  ncores <- max(1L, detectCores() - 1L)
  cl <- makeCluster(ncores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, library(fixest))
  clusterExport(cl, c("d_pois", "ctry_list", "ncl", "f_count"), envir = environment())
  clusterSetRNGStream(cl, 42)
  parSapply(cl, seq_len(B), boot_one)
}, error = function(e) { set.seed(42); sapply(seq_len(B), boot_one) })
boots <- boots[is.finite(boots)]
if (length(boots) > 50) {
  b_hat <- M2_ppml$coeftable["fp_Gc", "Estimate"]
  ci <- quantile(boots, c(.025, .975))
  p_boot <- 2 * min(mean(boots <= 0), mean(boots >= 0))
  cat(sprintf("  valid resamples = %d\n", length(boots)))
  cat(sprintf("  point b = %.5f | boot SE = %.5f | 95%% CI = [%.5f, %.5f] | p = %.4f\n",
              b_hat, sd(boots), ci[1], ci[2], p_boot))
}

cat("\n===== INFERENCE SUMMARY (interaction food price x Gini) =====\n")
cat("  model-based (glmmTMB NB) p :", round(get_int(M2)["p"], 4), "\n")
cat("  cluster-robust (PPML)    p :", round(pp_cl[4], 4), "\n")
cat("  Driscoll-Kraay (linear)  p :", round(se_dk[4], 4), "\n")
if (!is.null(wcb_webb)) cat("  wild bootstrap (Webb)    p :", round(wcb_webb$p_val, 4), "\n")
if (!is.null(sw))       cat("  PPML score WCB           p :", round(sw$p, 4), "\n")
if (length(boots) > 50) cat("  PPML pairs bootstrap     p :", round(2 * min(mean(boots <= 0), mean(boots >= 0)), 4), "\n")

cat("\n===== TABLE: MODEL COEFFICIENTS (M1 / M2 / M3, non-FE terms) =====\n")
dump_coefs <- function(model, name) {
  s <- summary(model)$coefficients$cond
  s <- s[!grepl("^iso3|^yearf|Intercept", rownames(s)), , drop = FALSE]
  cat(sprintf("\n--- %s (N=%d, AIC=%.1f, theta=%.3f) ---\n", name, nobs(model), AIC(model), sigma(model)))
  print(round(s, 5))
}
dump_coefs(M1, "Model 1 (linear price)")
dump_coefs(M2, "Model 2 (quadratic + interaction, CENTRAL)")
dump_coefs(M3, "Model 3 (+ lagged price)")

cat("\n===== TABLE: CONTROL-FUNCTION / IV RESULTS =====\n")
scf2 <- summary(CF2)$coefficients$cond
print(round(scf2[!grepl("^iso3|^yearf|Intercept", rownames(scf2)), , drop = FALSE], 5))
if (exists("IV2") && !is.null(IV2)) {
  cat("\nIV2 (2SLS) instrumented terms:\n")
  print(round(coef(IV2)[grep("fit_", names(coef(IV2)))], 5))
}

cat("\n===== EXPORTS =====\n")
b2 <- fixef(M2)$cond
v2 <- as.matrix(vcov(M2)$cond)
write.csv(data.frame(term = names(b2), estimate = as.numeric(b2)),
          file.path(OUTDIR, "m2_coef.csv"), row.names = FALSE)
write.csv(as.data.frame(v2), file.path(OUTDIR, "m2_vcov.csv"))
ctrl_means <- sapply(c("log_gdp", "urban_pct", "food_import_pct", "undernourish_pct",
                       "inflation_pct", "democracy"), function(v) mean(m[[v]], na.rm = TRUE))
meta <- data.frame(
  model              = "M2_glmmTMB_NB_quadratic_interaction",
  n_obs              = nobs(M2),
  n_countries        = length(unique(d_pois$iso3)),
  year_min           = min(d_pois$Year), year_max = max(d_pois$Year),
  gini_mean          = gini_mean,
  und_mean           = und_mean,
  food_infl_median   = fp_med,
  zero_crossing_gini = gini_mean + thr_gc,
  theta              = sigma(M2),
  t(ctrl_means)
)
write.csv(meta, file.path(OUTDIR, "m2_meta.csv"), row.names = FALSE)
cat("wrote m2_coef.csv, m2_vcov.csv, m2_meta.csv to", OUTDIR, "\n")

cat("\n--- forest-plot data (interaction across specifications) ---\n")
specs <- list(
  list(lab = "Central M2",            mod = M2,        term = "food_infl:G_c"),
  list(lab = "Winsorized + year FE",  mod = m2b,       term = "food_infl_w:G_c"),
  list(lab = "Exclude food_infl>50",  mod = m2n,       term = "food_infl:G_c"),
  list(lab = "Exclude food_infl>200", mod = m2_no,     term = "food_infl:G_c"),
  list(lab = "Exclude ZWE+SSD",       mod = m2q_noZS,  term = "food_infl:G_c"),
  list(lab = "+ gov_effect (l1)",     mod = M_gov,     term = "food_infl:G_c"),
  list(lab = "+ oop_health (l1)",     mod = M_oop,     term = "food_infl:G_c"),
  list(lab = "+ both",                mod = M_gov_oop, term = "food_infl:G_c")
)
fp <- do.call(rbind, lapply(specs, function(s) {
  v <- get_int(s$mod, s$term)
  data.frame(spec = s$lab, est = v["est"], se = v["se"], p = v["p"],
             lo = v["est"] - 1.96 * v["se"], hi = v["est"] + 1.96 * v["se"])
}))
rownames(fp) <- NULL
fp$irr <- exp(fp$est)
write.csv(fp, file.path(OUTDIR, "fig2_forest.csv"), row.names = FALSE)
print(round(fp[, c("est", "se", "p", "lo", "hi", "irr")], 5))

cat("\n--- Figure 2 data: marginal price slope across Gini (fixest NB) ---\n")
bfe <- coef(M2_fe)
Vm  <- vcov(M2_fe, vcov = "iid")
Vc  <- vcov(M2_fe, cluster = ~iso3)
gini_grid <- seq(quantile(m$gini, .05, na.rm = TRUE), quantile(m$gini, .95, na.rm = TRUE), length.out = 60)
gc_grid2  <- gini_grid - gini_mean
gvec <- function(gc) { g <- setNames(rep(0, length(bfe)), names(bfe)); g["food_infl"] <- 1; g["food_infl2"] <- 2 * fp_med; g["fp_Gc"] <- gc; g }
slope <- bfe["food_infl"] + 2 * bfe["food_infl2"] * fp_med + bfe["fp_Gc"] * gc_grid2
se_m  <- sapply(gc_grid2, function(gc) { g <- gvec(gc); sqrt(as.numeric(t(g) %*% Vm %*% g)) })
se_c  <- sapply(gc_grid2, function(gc) { g <- gvec(gc); sqrt(as.numeric(t(g) %*% Vc %*% g)) })
fig2 <- data.frame(gini = gini_grid, slope = slope,
                   lo_model = slope - 1.96 * se_m, hi_model = slope + 1.96 * se_m,
                   lo_clust = slope - 1.96 * se_c, hi_clust = slope + 1.96 * se_c)
write.csv(fig2, file.path(OUTDIR, "fig2_data.csv"), row.names = FALSE)
cat(sprintf("  marginal slope crosses zero at Gini = %.1f (model band excludes 0 above ~%.1f; cluster band wider)\n",
            gini_grid[which.min(abs(slope))], gini_grid[which.min(abs(slope))]))

cat("\n===== FIGURE 3: MAP OF COUNTRY-SPECIFIC PRICE SLOPES =====\n")
M_full <- base_fit(m)
bf <- fixef(M_full)$cond
ctry_g <- m %>% group_by(iso3) %>% summarise(gini = mean(gini, na.rm = TRUE), .groups = "drop")
ctry_g$slope <- bf["food_infl"] + 2 * bf["food_infl2"] * fp_med + bf["food_infl:G_c"] * (ctry_g$gini - gini_mean)
ctry_g$iso3 <- as.character(ctry_g$iso3)
ft <- tryCatch({
  if (!("Times New Roman" %in% fonts())) { font_import(pattern = "[Tt]imes", prompt = FALSE); loadfonts(quiet = TRUE) }
  if ("Times New Roman" %in% fonts()) "Times New Roman" else ""
}, error = function(e) "")
fig3 <- tryCatch({
  world  <- ne_countries(scale = "medium", returnclass = "sf")
  africa <- world[!is.na(world$continent) & world$continent == "Africa", ]
  map_df <- merge(africa, as.data.frame(ctry_g), by.x = "iso_a3", by.y = "iso3", all.x = TRUE)
  south4 <- c("ZAF", "NAM", "SWZ", "BWA")
  lab_df <- map_df[map_df$iso_a3 %in% south4, ]
  cent <- suppressWarnings(st_coordinates(st_centroid(lab_df)))
  lab_df$lon <- cent[, 1]; lab_df$lat <- cent[, 2]
  g <- ggplot(map_df) +
    geom_sf(aes(fill = slope), color = "grey40", linewidth = 0.1) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0, na.value = "grey90", name = "Price slope\n(at median inflation)") +
    geom_text_repel(data = lab_df, aes(x = lon, y = lat, label = iso_a3),
                    size = 3, family = ft, min.segment.length = 0, box.padding = 0.6,
                    segment.color = "grey20") +
    coord_sf(xlim = c(-20, 52), ylim = c(-37, 20), expand = FALSE) +
    labs(title = "Country-specific food-price slope on social unrest",
         subtitle = "Positive (red) effects concentrated in high-inequality Southern Africa",
         x = NULL, y = NULL) +
    theme_minimal(base_family = ft) +
    theme(legend.position = "right", panel.grid = element_line(color = "grey95"))
  ggsave(file.path(fig_dir, "fig3_map.png"), g, width = 8, height = 8, dpi = 300)
  ggsave(file.path(fig_dir, "fig3_map.pdf"), g, width = 8, height = 8, device = cairo_pdf)
  "ok"
}, error = function(e) { cat("  map step skipped:", conditionMessage(e), "\n"); NULL })
if (!is.null(fig3)) cat("  wrote fig3_map.png and fig3_map.pdf to", fig_dir, "\n")

cat("\n===== SAVE WORKSPACE =====\n")
key_objects <- c("M1", "M2", "M3", "M2_ss", "M3_ss", "M2_fe", "M2_ppml", "M_glmnb",
                 "P_M2", "CF2", "IV2", "IV_ss", "M_RQ1", "M_RQ1b", "M_RQ1c", "M_RQ1d",
                 "M_terc", "M_spl", "M_gov", "M_oop", "M_gov_oop", "M_full",
                 "ME_gini", "ME_und", "infl_tab", "gov_tab", "fp", "fig2", "meta",
                 "thr_gc", "jack", "boots", "desc")
key_objects <- key_objects[sapply(key_objects, exists)]
save(list = key_objects, file = file.path(OUTDIR, "ssa_workspace.RData"))
cat("saved objects:", paste(key_objects, collapse = ", "), "\n")

cat("\n===== DONE =====\n")
