library(plm)
library(readxl)
library(dplyr)
library(car)
library(glmmTMB)
library(splines)


OUTDIR <- "."
df <- read_excel("data/ssa_food_unrest_panel.xlsx")

df <- df %>% filter(iso3 != "SDN")

df <- df %>% arrange(iso3, Year)
pdata <- pdata.frame(df, index = c("iso3", "Year"))

core_vars <- c("unrest_n", "food_infl", "gini")

cat("=== PANEL STRUCTURE ===\n")
cat("Countries:", length(unique(df$iso3)), "| Years:", min(df$Year), "-", max(df$Year), "\n")
cat("Total obs:", nrow(df), "\n")
pdim(pdata)

cat("\n=== MISSINGNESS (core vars) ===\n")
for (v in core_vars) {
  cat(sprintf("%-12s non-missing: %d / %d (%.1f%%)\n",
              v, sum(!is.na(df[[v]])), nrow(df), 100*mean(!is.na(df[[v]]))))
}

cat("\n=== CROSS-SECTIONAL DEPENDENCE (Pesaran CD) ===\n")
cat("H0: cross-sectional independence\n\n")
for (v in core_vars) {
  f <- as.formula(paste(v, "~ 1"))
  cd <- tryCatch(pcdtest(f, data = pdata, test = "cd"),
                 error = function(e) NULL)
  if (!is.null(cd)) {
    cat(sprintf("%-12s CD = %8.3f | p = %.4g\n", v, cd$statistic, cd$p.value))
  } else {
    cat(sprintf("%-12s CD testi calismadi (eksik veri olabilir)\n", v))
  }
}


library(plm)

# H0: series has a unit root (non-stationary), against stationarity
# Pesaran CIPS — second-generation, robust to cross-sectional dependence

core_vars <- c("unrest_n", "food_infl", "gini")

run_cips <- function(varname, data, lags = 1) {
  x <- data[[varname]]
  res_level <- tryCatch(
    cipstest(x, lags = lags, type = "trend", model = "cmg", truncated = TRUE),
    error = function(e) NULL)
  res_drift <- tryCatch(
    cipstest(x, lags = lags, type = "drift", model = "cmg", truncated = TRUE),
    error = function(e) NULL)
  list(trend = res_level, drift = res_drift)
}

cat("=== CIPS PANEL UNIT-ROOT TEST (Pesaran, 2nd generation) ===\n")
cat("H0: unit root (non-stationary)  |  lags = 1\n")
cat("------------------------------------------------------------\n")

for (v in core_vars) {
  cat(sprintf("\n### %s\n", v))
  out <- run_cips(v, pdata, lags = 1)
  if (!is.null(out$drift)) {
    cat(sprintf("  LEVEL  (drift): CIPS = %7.3f | p %s\n",
                out$drift$statistic,
                ifelse(out$drift$p.value <= 0.01, "<= 0.01",
                       ifelse(out$drift$p.value >= 0.10, ">= 0.10",
                              sprintf("= %.3f", out$drift$p.value)))))
  }
  if (!is.null(out$trend)) {
    cat(sprintf("  LEVEL  (trend): CIPS = %7.3f | p %s\n",
                out$trend$statistic,
                ifelse(out$trend$p.value <= 0.01, "<= 0.01",
                       ifelse(out$trend$p.value >= 0.10, ">= 0.10",
                              sprintf("= %.3f", out$trend$p.value)))))
  }
}

cat("\n\n=== FIRST DIFFERENCES ===\n")
cat("H0: unit root (non-stationary)  |  lags = 1\n")
cat("------------------------------------------------------------\n")

for (v in core_vars) {
  dx <- diff(pdata[[v]])
  cat(sprintf("\n### d.%s\n", v))
  rd <- tryCatch(cipstest(dx, lags = 1, type = "drift", model = "cmg", truncated = TRUE),
                 error = function(e) NULL)
  if (!is.null(rd)) {
    cat(sprintf("  DIFF   (drift): CIPS = %7.3f | p %s\n",
                rd$statistic,
                ifelse(rd$p.value <= 0.01, "<= 0.01",
                       ifelse(rd$p.value >= 0.10, ">= 0.10",
                              sprintf("= %.3f", rd$p.value)))))
  }
}


library(dplyr)

cat("=== NA STRUCTURE: food_infl & gini ===\n\n")

for (v in c("food_infl", "gini")) {
  cat(sprintf("### %s\n", v))
  
  # NA per country
  na_by_country <- df %>%
    group_by(iso3) %>%
    summarise(n = n(), n_na = sum(is.na(.data[[v]])), .groups = "drop") %>%
    filter(n_na > 0) %>%
    arrange(desc(n_na))
  
  cat("Countries with any NA (count of NA years):\n")
  if (nrow(na_by_country) == 0) {
    cat("  none\n")
  } else {
    print(as.data.frame(na_by_country), row.names = FALSE)
  }
  
  # NA per year
  na_by_year <- df %>%
    group_by(Year) %>%
    summarise(n_na = sum(is.na(.data[[v]])), .groups = "drop") %>%
    filter(n_na > 0)
  cat("\nYears with any NA (count of NA countries):\n")
  print(as.data.frame(na_by_year), row.names = FALSE)
  
  # how many countries are FULLY observed (no NA at all)?
  full_countries <- df %>%
    group_by(iso3) %>%
    summarise(n_na = sum(is.na(.data[[v]])), .groups = "drop") %>%
    filter(n_na == 0) %>%
    nrow()
  cat(sprintf("\nFully-observed countries for %s: %d / 47\n", v, full_countries))
  cat("------------------------------------------------------------\n\n")
}




run_cips_clean <- function(varname, drop_iso, start_year = NULL) {
  d <- df %>% filter(!iso3 %in% drop_iso)
  if (!is.null(start_year)) d <- d %>% filter(Year >= start_year)
  d <- d %>% filter(!is.na(.data[[varname]])) %>% arrange(iso3, Year)
  
  # keep only countries with a usable run of observations
  ok <- d %>% group_by(iso3) %>% summarise(n = n(), .groups="drop") %>% filter(n >= 5)
  d <- d %>% filter(iso3 %in% ok$iso3)
  
  pd <- pdata.frame(d, index = c("iso3", "Year"))
  cat(sprintf("  [panel: %d countries, %d obs]\n",
              length(unique(d$iso3)), nrow(d)))
  
  for (ty in c("drift", "trend")) {
    r <- tryCatch(cipstest(pd[[varname]], lags = 1, type = ty, model = "cmg", truncated = TRUE),
                  error = function(e) NULL)
    if (!is.null(r)) {
      pv <- ifelse(r$p.value <= 0.01, "<= 0.01",
                   ifelse(r$p.value >= 0.10, ">= 0.10", sprintf("= %.3f", r$p.value)))
      cat(sprintf("  LEVEL (%-5s): CIPS = %7.3f | p %s\n", ty, r$statistic, pv))
    } else {
      cat(sprintf("  LEVEL (%-5s): calismadi\n", ty))
    }
  }
}

run_cips_clean_diff <- function(varname, drop_iso, start_year = NULL) {
  d <- df %>% filter(!iso3 %in% drop_iso)
  if (!is.null(start_year)) d <- d %>% filter(Year >= start_year)
  d <- d %>% filter(!is.na(.data[[varname]])) %>% arrange(iso3, Year)
  ok <- d %>% group_by(iso3) %>% summarise(n = n(), .groups="drop") %>% filter(n >= 6)
  d <- d %>% filter(iso3 %in% ok$iso3)
  pd <- pdata.frame(d, index = c("iso3", "Year"))
  dx <- diff(pd[[varname]])
  r <- tryCatch(cipstest(dx, lags = 1, type = "drift", model = "cmg", truncated = TRUE),
                error = function(e) NULL)
  if (!is.null(r)) {
    pv <- ifelse(r$p.value <= 0.01, "<= 0.01",
                 ifelse(r$p.value >= 0.10, ">= 0.10", sprintf("= %.3f", r$p.value)))
    cat(sprintf("  DIFF  (drift): CIPS = %7.3f | p %s\n", r$statistic, pv))
  }
}

cat("=== CIPS (clean sub-panels) — H0: unit root ===\n\n")

cat("### unrest_n (all 47 countries)\n")
run_cips_clean("unrest_n", drop_iso = character(0))
cat("\n")

cat("### food_infl (drop ERI, CAF; start 2001)\n")
run_cips_clean("food_infl", drop_iso = c("ERI","CAF"), start_year = 2001)
cat("\n")

cat("### gini (drop ERI)\n")
run_cips_clean("gini", drop_iso = c("ERI"))
cat("\n")

cat("=== FIRST DIFFERENCES ===\n\n")
cat("### d.food_infl\n")
run_cips_clean_diff("food_infl", drop_iso = c("ERI","CAF"), start_year = 2001)
cat("\n### d.gini\n")
run_cips_clean_diff("gini", drop_iso = c("ERI"))





gini_diag <- df %>%
  filter(iso3 != "ERI") %>%
  group_by(iso3) %>%
  summarise(
    n_obs   = sum(!is.na(gini)),
    sd_gini = sd(gini, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(sd_gini)

cat("=== gini: en dusuk varyansli 15 ulke ===\n")
print(as.data.frame(head(gini_diag, 15)), row.names = FALSE)

cat(sprintf("\nSabit (sd=0 veya NA) ulke sayisi: %d\n",
            sum(gini_diag$sd_gini == 0 | is.na(gini_diag$sd_gini))))
cat(sprintf("sd < 0.5 olan ulke sayisi: %d\n",
            sum(gini_diag$sd_gini < 0.5, na.rm = TRUE)))

cat(sprintf("\ngini_wb (ham, interpolasyonsuz) doluluk: %d / %d\n",
            sum(!is.na(df$gini_wb)), nrow(df)))









run_gini_cips <- function(drop_iso) {
  d <- df %>% filter(!iso3 %in% drop_iso) %>%
    filter(!is.na(gini)) %>% arrange(iso3, Year)
  pd <- pdata.frame(d, index = c("iso3", "Year"))
  cat(sprintf("  [panel: %d countries, %d obs]\n", length(unique(d$iso3)), nrow(d)))
  for (ty in c("drift", "trend")) {
    r <- tryCatch(cipstest(pd[["gini"]], lags = 1, type = ty, model = "cmg", truncated = TRUE),
                  error = function(e) NULL)
    if (!is.null(r)) {
      pv <- ifelse(r$p.value <= 0.01, "<= 0.01",
                   ifelse(r$p.value >= 0.10, ">= 0.10", sprintf("= %.3f", r$p.value)))
      cat(sprintf("  LEVEL (%-5s): CIPS = %7.3f | p %s\n", ty, r$statistic, pv))
    } else {
      cat(sprintf("  LEVEL (%-5s): calismadi\n", ty))
    }
  }
  dx <- diff(pd[["gini"]])
  rd <- tryCatch(cipstest(dx, lags = 1, type = "drift", model = "cmg", truncated = TRUE),
                 error = function(e) NULL)
  if (!is.null(rd)) {
    pv <- ifelse(rd$p.value <= 0.01, "<= 0.01",
                 ifelse(rd$p.value >= 0.10, ">= 0.10", sprintf("= %.3f", rd$p.value)))
    cat(sprintf("  DIFF  (drift): CIPS = %7.3f | p %s\n", rd$statistic, pv))
  }
}

cat("### gini (drop ERI, SOM)\n")
run_gini_cips(drop_iso = c("ERI", "SOM"))











run_dh <- function(data, yvar, xvar, lags, drop_iso = character(0),
                   start_year = NULL, label = "") {
  d <- data %>% filter(!iso3 %in% drop_iso) %>%
    select(iso3, Year, all_of(c(yvar, xvar))) %>%
    filter(!is.na(.data[[yvar]]) & !is.na(.data[[xvar]]))
  if (!is.null(start_year)) d <- d %>% filter(Year >= start_year)
  d <- d %>% arrange(iso3, Year)
  ok <- d %>% group_by(iso3) %>% summarise(n = n(), .groups = "drop") %>%
    filter(n >= (2 * lags + 3))
  d <- d %>% filter(iso3 %in% ok$iso3)
  pd <- pdata.frame(d, index = c("iso3", "Year"))
  f <- as.formula(paste(yvar, "~", xvar))
  r <- tryCatch(pgrangertest(f, data = pd, order = lags, test = "Ztilde"),
                error = function(e) {cat("   ERR:", conditionMessage(e), "\n"); NULL})
  cat(sprintf("\n### %s  [ %s  =/=>  %s ]\n", label, xvar, yvar))
  cat(sprintf("    panel: %d countries | lag = %d\n", length(unique(d$iso3)), lags))
  if (!is.null(r)) {
    stat <- as.numeric(r$statistic)
    cat(sprintf("    Ztilde = %7.3f | p = %.4g  -> %s\n",
                stat, r$p.value,
                ifelse(r$p.value < 0.05, "REJECT H0 (Granger-causes)",
                       "fail to reject (no causality)")))
  }
  invisible(r)
}

cat("==============================================================\n")
cat("DUMITRESCU-HURLIN (Ztilde, unbalanced-panel appropriate)\n")
cat("H0: x does NOT homogeneously Granger-cause y\n")
cat("==============================================================\n")

cat("\n--- MAIN: food_infl <-> unrest_n (levels, both I(0)) ---\n")
for (L in 1:3) {
  cat(sprintf("\n[ lag = %d ]", L))
  run_dh(df, "unrest_n", "food_infl", L, drop_iso = c("ERI","CAF"),
         start_year = 2001, label = sprintf("H5 (L=%d)", L))
  run_dh(df, "food_infl", "unrest_n", L, drop_iso = c("ERI","CAF"),
         start_year = 2001, label = sprintf("H6 (L=%d)", L))
}

df <- df %>%
  arrange(iso3, Year) %>%
  group_by(iso3) %>%
  mutate(d_gini = gini - dplyr::lag(gini)) %>%
  ungroup()

cat("d_gini eklendi. non-missing:", sum(!is.na(df$d_gini)), "/", nrow(df), "\n\n")

cat("--- ROBUSTNESS: d_gini -> unrest_n (gini I(1)) ---\n")
for (L in 1:2) {
  run_dh(df, "unrest_n", "d_gini", L, drop_iso = c("ERI","SOM"),
         label = sprintf("H7-rob (L=%d)", L))
}


vif_data <- df %>%
  filter(iso3 != "SDN") %>%
  select(food_infl, gini, gdp_pc, urban_pct, food_import_pct,
         undernourish_pct, inflation_pct, democracy) %>%
  na.omit()

cat("VIF sample size:", nrow(vif_data), "\n\n")

lm_vif <- lm(food_infl ~ gini + gdp_pc + urban_pct + food_import_pct +
               undernourish_pct + inflation_pct + democracy,
             data = vif_data)

cat("=== VIF (main-model covariates) ===\n")
print(round(vif(lm_vif), 3))

cat("\n=== with gov_effect added (check overlap with democracy) ===\n")
vif_data2 <- df %>%
  filter(iso3 != "SDN") %>%
  select(food_infl, gini, gdp_pc, urban_pct, food_import_pct,
         undernourish_pct, inflation_pct, democracy, gov_effect) %>%
  na.omit()
cat("sample size:", nrow(vif_data2), "\n")
lm_vif2 <- lm(food_infl ~ gini + gdp_pc + urban_pct + food_import_pct +
                undernourish_pct + inflation_pct + democracy + gov_effect,
              data = vif_data2)
print(round(vif(lm_vif2), 3))


m <- df %>%
  filter(iso3 != "SDN") %>%
  arrange(iso3, Year) %>%
  mutate(
    iso3   = factor(iso3),
    yearf  = factor(Year),
    G_c    = gini - mean(gini, na.rm = TRUE),          # mean-centred Gini
    logP   = log(population)
  ) %>%
  group_by(iso3) %>%
  mutate(food_infl_l1 = dplyr::lag(food_infl)) %>%
  ungroup()

cat("=== FE NEGATIVE BINOMIAL (glmmTMB, nbinom2) ===\n")
cat("offset = log(population); FE = country + year\n\n")

## Model 1: price only
m1 <- glmmTMB(unrest_n ~ food_infl + gdp_pc + urban_pct + food_import_pct +
                undernourish_pct + inflation_pct + democracy +
                iso3 + yearf + offset(logP),
              family = nbinom2, data = m)

## Model 2: + inequality + interaction (central)
m2 <- glmmTMB(unrest_n ~ food_infl * G_c + gdp_pc + urban_pct + food_import_pct +
                undernourish_pct + inflation_pct + democracy +
                iso3 + yearf + offset(logP),
              family = nbinom2, data = m)

## Model 3: + lagged price + lagged interaction
m3 <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl_l1 * G_c +
                gdp_pc + urban_pct + food_import_pct +
                undernourish_pct + inflation_pct + democracy +
                iso3 + yearf + offset(logP),
              family = nbinom2, data = m)

extract_core <- function(model, terms, mname) {
  s <- summary(model)$coefficients$cond
  cat(sprintf("\n--- %s | N = %d | logLik = %.1f | AIC = %.1f ---\n",
              mname, nobs(model), as.numeric(logLik(model)), AIC(model)))
  for (tm in terms) {
    if (tm %in% rownames(s)) {
      est <- s[tm, "Estimate"]; se <- s[tm, "Std. Error"]
      p   <- s[tm, "Pr(>|z|)"]
      cat(sprintf("  %-22s b = %8.4f | SE = %7.4f | IRR = %6.3f | p = %.4g\n",
                  tm, est, se, exp(est), p))
    }
  }
}

extract_core(m1, c("food_infl"), "MODEL 1 (price only)")
extract_core(m2, c("food_infl","G_c","food_infl:G_c"), "MODEL 2 (+ interaction)")
extract_core(m3, c("food_infl","food_infl_l1","G_c",
                   "food_infl:G_c","food_infl_l1:G_c"), "MODEL 3 (+ lag)")













cat("=== food_infl distribution ===\n")
print(summary(m$food_infl))
cat("quantiles:\n")
print(round(quantile(m$food_infl, c(0,.01,.05,.25,.5,.75,.95,.99,1), na.rm=TRUE), 2))
cat("\nobs with food_infl > 50:", sum(m$food_infl > 50, na.rm=TRUE), "\n")
cat("obs with food_infl > 100:", sum(m$food_infl > 100, na.rm=TRUE), "\n")

# how much food_infl variation is BETWEEN years vs within?
cat("\n=== variance decomposition of food_infl ===\n")
av <- anova(lm(food_infl ~ factor(Year), data = m))
cat("share of food_infl variance explained by YEAR dummies:",
    round(av[1,"Sum Sq"]/sum(av[,"Sum Sq"]), 3), "\n")
avc <- anova(lm(food_infl ~ factor(iso3), data = m))
cat("share explained by COUNTRY dummies:",
    round(avc[1,"Sum Sq"]/sum(avc[,"Sum Sq"]), 3), "\n")

# Model 2 variants: (a) WITHOUT year FE, (b) winsorized price, (c) food_shock instead
m$food_infl_w <- pmin(pmax(m$food_infl, quantile(m$food_infl,.01,na.rm=TRUE)),
                      quantile(m$food_infl,.99,na.rm=TRUE))
m$log_gdp <- log(m$gdp_pc)

cat("\n=== Model 2 variants: interaction food_infl:G_c ===\n")

# (a) drop year FE
m2a <- glmmTMB(unrest_n ~ food_infl * G_c + log_gdp + urban_pct + food_import_pct +
                 undernourish_pct + inflation_pct + democracy +
                 iso3 + offset(logP), family = nbinom2, data = m)
sa <- summary(m2a)$coefficients$cond
cat(sprintf("(a) NO year FE     : int b = %.4f | p = %.4g\n",
            sa["food_infl:G_c","Estimate"], sa["food_infl:G_c","Pr(>|z|)"]))

# (b) winsorized price, keep year FE
m2b <- glmmTMB(unrest_n ~ food_infl_w * G_c + log_gdp + urban_pct + food_import_pct +
                 undernourish_pct + inflation_pct + democracy +
                 iso3 + yearf + offset(logP), family = nbinom2, data = m)
sb <- summary(m2b)$coefficients$cond
cat(sprintf("(b) winsor + yearFE: int b = %.4f | p = %.4g\n",
            sb["food_infl_w:G_c","Estimate"], sb["food_infl_w:G_c","Pr(>|z|)"]))

# (c) winsorized price, NO year FE
m2c <- glmmTMB(unrest_n ~ food_infl_w * G_c + log_gdp + urban_pct + food_import_pct +
                 undernourish_pct + inflation_pct + democracy +
                 iso3 + offset(logP), family = nbinom2, data = m)
sc <- summary(m2c)$coefficients$cond
cat(sprintf("(c) winsor, NO yrFE: int b = %.4f | p = %.4g\n",
            sc["food_infl_w:G_c","Estimate"], sc["food_infl_w:G_c","Pr(>|z|)"]))


m <- m %>%
  mutate(
    log_gdp = log(gdp_pc),
    food_infl_w = pmin(pmax(food_infl,
                            quantile(food_infl, 0.01, na.rm = TRUE)),
                       quantile(food_infl, 0.99, na.rm = TRUE))
  )

cat("=== food_infl variance: between-year vs between-country ===\n")
av_y <- anova(lm(food_infl ~ factor(Year),  data = m))
av_c <- anova(lm(food_infl ~ factor(iso3),  data = m))
cat(sprintf("  share explained by YEAR    : %.3f\n", av_y[1,"Sum Sq"]/sum(av_y[,"Sum Sq"])))
cat(sprintf("  share explained by COUNTRY : %.3f\n", av_c[1,"Sum Sq"]/sum(av_c[,"Sum Sq"])))
cat(sprintf("  food_infl winsor range: [%.1f, %.1f] (was [%.1f, %.1f])\n",
            min(m$food_infl_w,na.rm=TRUE), max(m$food_infl_w,na.rm=TRUE),
            min(m$food_infl,na.rm=TRUE),   max(m$food_infl,na.rm=TRUE)))

cat("\n=== FE-NB with log(gdp_pc) ===\n")

m1 <- glmmTMB(unrest_n ~ food_infl + log_gdp + urban_pct + food_import_pct +
                undernourish_pct + inflation_pct + democracy +
                iso3 + yearf + offset(logP), family = nbinom2, data = m)

m2 <- glmmTMB(unrest_n ~ food_infl * G_c + log_gdp + urban_pct + food_import_pct +
                undernourish_pct + inflation_pct + democracy +
                iso3 + yearf + offset(logP), family = nbinom2, data = m)

m3 <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl_l1 * G_c +
                log_gdp + urban_pct + food_import_pct +
                undernourish_pct + inflation_pct + democracy +
                iso3 + yearf + offset(logP), family = nbinom2, data = m)

extract_core(m1, c("food_infl"), "MODEL 1 (price only, log gdp)")
extract_core(m2, c("food_infl","G_c","food_infl:G_c"), "MODEL 2 (+ interaction, log gdp)")
extract_core(m3, c("food_infl","food_infl_l1","G_c",
                   "food_infl:G_c","food_infl_l1:G_c"), "MODEL 3 (+ lag, log gdp)")


cat("=== WHO are the hyperinflation observations (food_infl > 50)? ===\n")
hyper <- m %>% filter(food_infl > 50) %>%
  select(iso3, Year, food_infl, gini, unrest_n) %>%
  arrange(desc(food_infl))
print(as.data.frame(hyper), row.names = FALSE)

cat("\n=== which countries do they belong to? ===\n")
print(table(droplevels(hyper$iso3)))

cat("\n=== gini at hyperinflation obs vs overall ===\n")
cat(sprintf("  mean gini at food_infl>50 : %.1f\n", mean(hyper$gini, na.rm=TRUE)))
cat(sprintf("  mean gini overall         : %.1f\n", mean(m$gini, na.rm=TRUE)))

cat("\n=== NON-LINEARITY: quadratic price term (Model 2 + food_infl^2) ===\n")
m$food_infl2 <- m$food_infl^2
m2q <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 +
                 log_gdp + urban_pct + food_import_pct +
                 undernourish_pct + inflation_pct + democracy +
                 iso3 + yearf + offset(logP), family = nbinom2, data = m)
sq <- summary(m2q)$coefficients$cond
for (tm in c("food_infl","food_infl2","G_c","food_infl:G_c")) {
  if (tm %in% rownames(sq))
    cat(sprintf("  %-16s b = %10.5f | p = %.4g\n", tm, sq[tm,"Estimate"], sq[tm,"Pr(>|z|)"]))
}

cat("\n=== Model 2 on NORMAL-inflation subsample (food_infl < 50) ===\n")
m_norm <- m %>% filter(food_infl < 50)
m2n <- glmmTMB(unrest_n ~ food_infl * G_c + log_gdp + urban_pct + food_import_pct +
                 undernourish_pct + inflation_pct + democracy +
                 iso3 + yearf + offset(logP), family = nbinom2, data = m_norm)
sn <- summary(m2n)$coefficients$cond
cat(sprintf("  N = %d (dropped %d hyperinflation obs)\n", nobs(m2n), nobs(m2)-nobs(m2n)))
cat(sprintf("  food_infl      b = %.4f | p = %.4g\n", sn["food_infl","Estimate"], sn["food_infl","Pr(>|z|)"]))
cat(sprintf("  food_infl:G_c  b = %.4f | p = %.4g\n", sn["food_infl:G_c","Estimate"], sn["food_infl:G_c","Pr(>|z|)"]))

cat("\n=== sensitivity: drop the single most extreme obs (601%) ===\n")
m_no601 <- m %>% filter(food_infl < 200)
m2_no <- glmmTMB(unrest_n ~ food_infl * G_c + log_gdp + urban_pct + food_import_pct +
                   undernourish_pct + inflation_pct + democracy +
                   iso3 + yearf + offset(logP), family = nbinom2, data = m_no601)
sno <- summary(m2_no)$coefficients$cond
cat(sprintf("  drop food_infl>200 (N=%d): food_infl:G_c b = %.4f | p = %.4g\n",
            nobs(m2_no), sno["food_infl:G_c","Estimate"], sno["food_infl:G_c","Pr(>|z|)"]))

# Is the quadratic robust, or just fitting Zimbabwe/South Sudan?
cat("=== quadratic WITHOUT ZWE & SSD (the two hyperinflation countries) ===\n")
m_noZS <- m %>% filter(!iso3 %in% c("ZWE","SSD"))
m2q_noZS <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 +
                      log_gdp + urban_pct + food_import_pct +
                      undernourish_pct + inflation_pct + democracy +
                      iso3 + yearf + offset(logP), family = nbinom2, data = m_noZS)
sq2 <- summary(m2q_noZS)$coefficients$cond
cat(sprintf("  N = %d (dropped ZWE+SSD)\n", nobs(m2q_noZS)))
for (tm in c("food_infl","food_infl2","G_c","food_infl:G_c")) {
  if (tm %in% rownames(sq2))
    cat(sprintf("  %-16s b = %11.6f | p = %.4g\n", tm, sq2[tm,"Estimate"], sq2[tm,"Pr(>|z|)"]))
}

# Where is the turning point of the inverted-U? (in food_infl units)
b1 <- summary(m2q)$coefficients$cond["food_infl","Estimate"]
b2 <- summary(m2q)$coefficients$cond["food_infl2","Estimate"]
cat(sprintf("\n  inverted-U turning point at food_infl = %.1f%%\n", -b1/(2*b2)))
cat("  (above this inflation level, marginal effect on unrest declines)\n")

# Compare model fit: linear vs quadratic (AIC)
cat("\n=== model comparison (AIC) ===\n")
cat(sprintf("  Model 2 linear    : AIC = %.1f\n", AIC(m2)))
cat(sprintf("  Model 2 quadratic : AIC = %.1f\n", AIC(m2q)))
cat(sprintf("  improvement       : %.1f\n", AIC(m2) - AIC(m2q)))

# LR test
cat("\n=== LR test: does quadratic significantly improve fit? ===\n")
print(anova(m2, m2q))



cat("=== FINAL MAIN MODELS (quadratic price, log gdp, country+year FE) ===\n")
cat("    raw food_infl (not winsorized); non-linearity via food_infl^2\n\n")

# Model 1: linear price only (baseline)
M1 <- glmmTMB(unrest_n ~ food_infl + log_gdp + urban_pct + food_import_pct +
                undernourish_pct + inflation_pct + democracy +
                iso3 + yearf + offset(logP), family = nbinom2, data = m)

# Model 2: quadratic price + inequality interaction (CENTRAL MODEL)
M2 <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 +
                log_gdp + urban_pct + food_import_pct +
                undernourish_pct + inflation_pct + democracy +
                iso3 + yearf + offset(logP), family = nbinom2, data = m)

# Model 3: + lagged price + lagged interaction
m <- m %>% mutate(food_infl2_l1 = food_infl_l1^2)
M3 <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 +
                food_infl_l1 * G_c + food_infl2_l1 +
                log_gdp + urban_pct + food_import_pct +
                undernourish_pct + inflation_pct + democracy +
                iso3 + yearf + offset(logP), family = nbinom2, data = m)

extract_core(M1, c("food_infl"), "MODEL 1 (linear price)")
extract_core(M2, c("food_infl","food_infl2","G_c","food_infl:G_c"),
             "MODEL 2 (quadratic + interaction) *** CENTRAL ***")
extract_core(M3, c("food_infl","food_infl2","food_infl_l1","food_infl2_l1",
                   "G_c","food_infl:G_c","food_infl_l1:G_c"),
             "MODEL 3 (+ lag)")

# dispersion parameter (theta) for each
cat("\n=== dispersion (theta) ===\n")
cat(sprintf("  M1 theta = %.3f | M2 theta = %.3f | M3 theta = %.3f\n",
            sigma(M1), sigma(M2), sigma(M3)))

# save full coefficient tables for later (Results table)
cat("\n=== full M2 coefficient table (all non-FE terms) ===\n")
s2 <- summary(M2)$coefficients$cond
keep <- !grepl("^iso3|^yearf|Intercept", rownames(s2))
print(round(s2[keep, ], 5))




# Model 3 fix: drop lagged quadratic (caused collinearity), keep lagged linear
M3 <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 +
                food_infl_l1 * G_c +
                log_gdp + urban_pct + food_import_pct +
                undernourish_pct + inflation_pct + democracy +
                iso3 + yearf + offset(logP), family = nbinom2, data = m)

cat("=== MODEL 3 (fixed: lagged linear only, no lagged quadratic) ===\n")
if (!is.null(M3$fit$convergence) && M3$fit$convergence == 0) cat("  converged OK\n")
extract_core(M3, c("food_infl","food_infl2","food_infl_l1",
                   "G_c","food_infl:G_c","food_infl_l1:G_c"), "MODEL 3 (fixed)")
cat(sprintf("  N = %d | AIC = %.1f | theta = %.3f\n", nobs(M3), AIC(M3), sigma(M3)))

# check Hessian / SE validity
s3 <- summary(M3)$coefficients$cond
cat("\n  any NaN in SE?", any(is.nan(s3[,"Std. Error"])), "\n")



library(fixest)

cat("===== 1. M3 lagged interaction =====\n")
print(rownames(summary(M3)$coefficients$cond)[!grepl("^iso3|^yearf", rownames(summary(M3)$coefficients$cond))])
lag_int <- grep("food_infl_l1:G_c|G_c:food_infl_l1", names(fixef(M3)$cond), value = TRUE)
cat("lagged interaction term:", lag_int, "\n")
s3 <- summary(M3)$coefficients$cond
print(round(s3[c("food_infl","food_infl2","food_infl_l1","G_c","food_infl:G_c",lag_int), ], 5))

cat("\n===== 2. H3 joint Wald =====\n")
wald_terms <- function(model, terms) {
  b <- fixef(model)$cond; V <- vcov(model)$cond
  terms <- terms[terms %in% names(b)]
  bs <- b[terms]; Vs <- V[terms, terms, drop = FALSE]
  W <- as.numeric(t(bs) %*% solve(Vs) %*% bs)
  data.frame(terms = paste(terms, collapse=" + "), chisq = W,
             df = length(terms), p = pchisq(W, length(terms), lower.tail = FALSE))
}
print(wald_terms(M3, c("food_infl_l1", lag_int)))

cat("\n===== 3. Same-sample M2 vs M3 =====\n")
m3_data <- m %>% filter(!is.na(food_infl_l1)) %>% droplevels()
M2_ss <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 +
                   log_gdp + urban_pct + food_import_pct + undernourish_pct + inflation_pct +
                   democracy + iso3 + yearf + offset(logP), family = nbinom2, data = m3_data)
M3_ss <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 + food_infl_l1 * G_c +
                   log_gdp + urban_pct + food_import_pct + undernourish_pct + inflation_pct +
                   democracy + iso3 + yearf + offset(logP), family = nbinom2, data = m3_data)
print(AIC(M2_ss, M3_ss))
print(anova(M2_ss, M3_ss))

gini_mean <- mean(m$gini, na.rm=TRUE)

cat("\n===== 4. Marginal food-price slope by Gini x inflation =====\n")
marg_fp <- function(model, FP, GC, intname = "food_infl:G_c") {
  b <- fixef(model)$cond; V <- vcov(model)$cond
  out <- list()
  for (fp in FP) for (gc in GC) {
    g <- setNames(rep(0, length(b)), names(b))
    g["food_infl"] <- 1; g["food_infl2"] <- 2*fp; g[intname] <- gc
    est <- b["food_infl"] + 2*b["food_infl2"]*fp + b[intname]*gc
    se  <- sqrt(as.numeric(t(g) %*% V %*% g))
    out[[length(out)+1]] <- data.frame(FP = fp, G_c = gc,
                                       slope = as.numeric(est), SE = se, IRR = exp(as.numeric(est)),
                                       p = 2*pnorm(abs(est/se), lower.tail = FALSE))
  }
  do.call(rbind, out)
}
g_q  <- quantile(m$gini, c(.25,.50,.75,.90), na.rm = TRUE)
fp_q <- quantile(m$food_infl, c(.25,.50,.75), na.rm = TRUE)
ME_gini <- marg_fp(M2, as.numeric(fp_q), as.numeric(g_q - gini_mean))
ME_gini$gini_pct <- rep(names(g_q), each = length(fp_q))
ME_gini$gini_val <- rep(as.numeric(g_q), each = length(fp_q))
nc <- sapply(ME_gini, is.numeric); ME_gini[nc] <- round(ME_gini[nc], 5)
print(ME_gini)



cat("\n===== 5. Gini threshold for positive food-price effect =====\n")
b2 <- fixef(M2)$cond
thr <- data.frame(FP = as.numeric(fp_q),
                  Gini_threshold = as.numeric(gini_mean - (b2["food_infl"] + 2*b2["food_infl2"]*as.numeric(fp_q)) / b2["food_infl:G_c"]))
print(round(thr, 3))

cat("\n===== 6. RQ1 nutritional vulnerability model =====\n")
m <- m %>% arrange(iso3, Year) %>% group_by(iso3) %>%
  mutate(undernourish_l1 = dplyr::lag(undernourish_pct)) %>% ungroup() %>%
  mutate(UND_l1_c = undernourish_l1 - mean(undernourish_l1, na.rm = TRUE))
M_RQ1 <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 + food_infl * UND_l1_c +
                   log_gdp + urban_pct + food_import_pct + inflation_pct + democracy +
                   iso3 + yearf + offset(logP), family = nbinom2, data = m)
s_rq1 <- summary(M_RQ1)$coefficients$cond
und_int <- grep("food_infl:UND_l1_c|UND_l1_c:food_infl", rownames(s_rq1), value = TRUE)
print(round(s_rq1[c("food_infl","food_infl2","G_c","food_infl:G_c","UND_l1_c",und_int), ], 5))
cat(sprintf("N = %d | AIC = %.1f | theta = %.3f\n", nobs(M_RQ1), AIC(M_RQ1), sigma(M_RQ1)))

cat("\n===== 7. Marginal food-price slope by undernourishment =====\n")
und_q <- quantile(m$undernourish_l1, c(.25,.50,.75,.90), na.rm = TRUE)
und_mean <- mean(m$undernourish_l1, na.rm = TRUE)
ME_und <- marg_fp(M_RQ1, as.numeric(fp_q), as.numeric(und_q - und_mean), intname = und_int)
ME_und$und_pct <- rep(names(und_q), each = length(fp_q))
nc <- sapply(ME_und, is.numeric); ME_und[nc] <- round(ME_und[nc], 5)
print(ME_und)

cat("\n===== 8. gov_effect / oop_health sensitivity =====\n")
m <- m %>% arrange(iso3, Year) %>% group_by(iso3) %>%
  mutate(gov_effect_l1 = dplyr::lag(gov_effect),
         oop_health_l1 = dplyr::lag(oop_health_pct)) %>% ungroup()
mk <- function(extra) glmmTMB(reformulate(c("food_infl * G_c","food_infl2","log_gdp",
                                            "urban_pct","food_import_pct","undernourish_pct","inflation_pct","democracy",
                                            extra,"iso3","yearf","offset(logP)"), response = "unrest_n"),
                              family = nbinom2, data = m)
M_gov     <- mk("gov_effect_l1")
M_oop     <- mk("oop_health_l1")
M_gov_oop <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 + log_gdp + urban_pct +
                       food_import_pct + undernourish_pct + inflation_pct + democracy +
                       gov_effect_l1 + oop_health_l1 + iso3 + yearf + offset(logP), family = nbinom2, data = m)
extr <- function(mod, nm) { s <- summary(mod)$coefficients$cond
data.frame(model = nm, N = nobs(mod), AIC = round(AIC(mod),1),
           est = s["food_infl:G_c","Estimate"], SE = s["food_infl:G_c","Std. Error"],
           IRR = exp(s["food_infl:G_c","Estimate"]), p = s["food_infl:G_c","Pr(>|z|)"]) }
print(rbind(extr(M2,"Central"), extr(M_gov,"+gov"), extr(M_oop,"+oop"),
            extr(M_gov_oop,"+gov+oop")) %>% mutate(across(where(is.numeric), ~round(.,5))))

cat("\n===== 9. VIF oop / gov+oop =====\n")
vd1 <- m %>% select(food_infl, gini, log_gdp, urban_pct, food_import_pct,
                    undernourish_pct, inflation_pct, democracy, oop_health_l1) %>% na.omit()
cat("oop only, N =", nrow(vd1), "\n")
print(round(car::vif(lm(food_infl ~ gini + log_gdp + urban_pct + food_import_pct +
                          undernourish_pct + inflation_pct + democracy + oop_health_l1, data = vd1)), 3))
vd2 <- m %>% select(food_infl, gini, log_gdp, urban_pct, food_import_pct,
                    undernourish_pct, inflation_pct, democracy, gov_effect_l1, oop_health_l1) %>% na.omit()
cat("gov+oop, N =", nrow(vd2), "\n")
print(round(car::vif(lm(food_infl ~ gini + log_gdp + urban_pct + food_import_pct +
                          undernourish_pct + inflation_pct + democracy + gov_effect_l1 + oop_health_l1, data = vd2)), 3))

cat("\n===== 10. PPML (fixest) =====\n")
PPML_M2 <- fepois(unrest_n ~ food_infl * G_c + food_infl2 + log_gdp + urban_pct +
                    food_import_pct + undernourish_pct + inflation_pct + democracy +
                    offset(logP) | iso3 + yearf, data = m, cluster = ~ iso3)
print(round(coeftable(PPML_M2)[c("food_infl","food_infl2","G_c","food_infl:G_c"), ], 5))

cat("\n===== 11. Poisson vs NB overdispersion =====\n")
P_M2 <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 + log_gdp + urban_pct +
                  food_import_pct + undernourish_pct + inflation_pct + democracy +
                  iso3 + yearf + offset(logP), family = poisson, data = m)
print(AIC(P_M2, M2))
print(anova(P_M2, M2))

cat("\n===== 12. Influence / nonlinearity H2 table =====\n")
ek <- function(mod, term, lab) { s <- summary(mod)$coefficients$cond
data.frame(spec = lab, term = term, N = nobs(mod), AIC = round(AIC(mod),1),
           est = s[term,"Estimate"], SE = s[term,"Std. Error"], IRR = exp(s[term,"Estimate"]),
           p = s[term,"Pr(>|z|)"], decision = ifelse(s[term,"Pr(>|z|)"]<.05,"Reject H0","Fail")) }
infl_tab <- rbind(
  ek(M2,"food_infl:G_c","Central nonlinear FE-NB"),
  ek(m2b,"food_infl_w:G_c","Winsor 1-99 + yearFE"),
  ek(m2c,"food_infl_w:G_c","Winsor 1-99 no yearFE"),
  ek(m2n,"food_infl:G_c","Exclude food_infl>50"),
  ek(m2_no,"food_infl:G_c","Exclude food_infl>200"),
  ek(m2q_noZS,"food_infl:G_c","Quadratic excl ZWE+SSD"))
print(infl_tab %>% mutate(across(where(is.numeric), ~round(.,5))))

cat("\n===== 13. Hyperinflation raw vs complete-case =====\n")
cc_vars <- c("unrest_n","food_infl","G_c","food_infl2","log_gdp","urban_pct",
             "food_import_pct","undernourish_pct","inflation_pct","democracy","logP")
mcc <- m %>% filter(complete.cases(across(all_of(cc_vars))))
print(data.frame(threshold = c(">50",">100",">200"),
                 raw_n = c(sum(m$food_infl>50,na.rm=T), sum(m$food_infl>100,na.rm=T), sum(m$food_infl>200,na.rm=T)),
                 complete_case_n = c(sum(mcc$food_infl>50,na.rm=T), sum(mcc$food_infl>100,na.rm=T), sum(mcc$food_infl>200,na.rm=T))))

cat("\n===== 14. Control-function NB (linear, H4) =====\n")
base_imp <- m %>% filter(Year >= 2000, Year <= 2004) %>% group_by(iso3) %>%
  summarise(import_base = mean(food_import_pct, na.rm = TRUE), .groups = "drop")
iv_data <- m %>%
  filter(!is.na(food_infl), !is.na(G_c), !is.na(global_food_infl), !is.na(temp_anomaly),
         !is.na(food_import_pct), !is.na(log_gdp), !is.na(urban_pct), !is.na(undernourish_pct),
         !is.na(inflation_pct), !is.na(democracy), !is.na(unrest_n), !is.na(logP)) %>%
  arrange(iso3, Year) %>% left_join(base_imp, by = "iso3") %>%
  mutate(Z1 = global_food_infl * import_base, Z2 = temp_anomaly,
         Z1G = Z1 * G_c, Z2G = Z2 * G_c, FPG = food_infl * G_c)
iv_data <- iv_data %>%
  filter(!is.na(Z1), !is.na(Z2), !is.na(Z1G), !is.na(Z2G), !is.na(import_base)) %>%
  droplevels()
cat("iv_data rows after instrument NA filter:", nrow(iv_data), "\n")

FS <- lm(food_infl ~ Z1 + Z2 + Z1G + Z2G + G_c + log_gdp + urban_pct + food_import_pct +
           undernourish_pct + inflation_pct + democracy + iso3 + yearf, data = iv_data)
iv_data$nu_hat <- resid(FS)
cat("first-stage F (instruments jointly):\n")
print(car::linearHypothesis(FS, c("Z1=0","Z2=0","Z1G=0","Z2G=0")))

CF_NB <- glmmTMB(unrest_n ~ food_infl * G_c + nu_hat + nu_hat:G_c +
                   log_gdp + urban_pct + food_import_pct + undernourish_pct + inflation_pct +
                   democracy + iso3 + yearf + offset(logP), family = nbinom2, data = iv_data)
s_cf <- summary(CF_NB)$coefficients$cond
nu_int <- grep("nu_hat:G_c|G_c:nu_hat", rownames(s_cf), value = TRUE)
print(round(s_cf[c("food_infl","G_c","food_infl:G_c","nu_hat",nu_int), ], 5))

cat("===== 15. 2SLS benchmark (hatasiz) =====\n")

iv_data$u_rate <- log1p((iv_data$unrest_n / iv_data$population) * 1e6)
cat("u_rate eklendi, range:", round(range(iv_data$u_rate, na.rm=TRUE),2), "\n")
cat("iv_data rows:", nrow(iv_data), "| u_rate NA:", sum(is.na(iv_data$u_rate)), "\n\n")

IV_2SLS <- tryCatch(
  feols(
    u_rate ~ log_gdp + urban_pct + food_import_pct + undernourish_pct +
      inflation_pct + democracy | iso3 + yearf |
      food_infl + FPG ~ Z1 + Z2 + Z1G + Z2G,
    data = iv_data, cluster = ~ iso3
  ),
  error = function(e) { cat("2SLS error:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(IV_2SLS)) {
  print(summary(IV_2SLS, stage = 2))
  cat("\n--- first-stage diagnostics ---\n")
  print(fitstat(IV_2SLS, type = c("ivf","ivwald","cd")))
} else {
  cat("\n2SLS kurulamadi; control-function (blok 14) zaten H4 icin yeterli.\n")
}

cat("\n===== DONE =====\n")

cat("===== IV-v2: global_food_infl as direct instrument =====\n")
iv2 <- m %>%
  filter(!is.na(food_infl), !is.na(G_c), !is.na(global_food_infl),
         !is.na(log_gdp), !is.na(urban_pct), !is.na(food_import_pct),
         !is.na(undernourish_pct), !is.na(inflation_pct), !is.na(democracy),
         !is.na(unrest_n), !is.na(logP)) %>%
  arrange(iso3, Year) %>%
  mutate(Zg = global_food_infl, ZgG = global_food_infl * G_c, FPG = food_infl * G_c) %>%
  droplevels()
cat("iv2 rows:", nrow(iv2), "\n")

FS2 <- lm(food_infl ~ Zg + ZgG + G_c + log_gdp + urban_pct + food_import_pct +
            undernourish_pct + inflation_pct + democracy + iso3 + yearf, data = iv2)
cat("first-stage F for food_infl (instruments Zg, ZgG jointly):\n")
print(car::linearHypothesis(FS2, c("Zg=0","ZgG=0")))
iv2$nu_hat <- resid(FS2)

CF2 <- glmmTMB(unrest_n ~ food_infl * G_c + nu_hat + nu_hat:G_c +
                 log_gdp + urban_pct + food_import_pct + undernourish_pct + inflation_pct +
                 democracy + iso3 + yearf + offset(logP), family = nbinom2, data = iv2)
s_cf2 <- summary(CF2)$coefficients$cond
nu_int2 <- grep("nu_hat:G_c|G_c:nu_hat", rownames(s_cf2), value = TRUE)
cat("\ncontrol-function (IV-v2):\n")
print(round(s_cf2[c("food_infl","G_c","food_infl:G_c","nu_hat",nu_int2), ], 5))

cat("\n2SLS (IV-v2):\n")
iv2$u_rate <- log1p((iv2$unrest_n / iv2$population) * 1e6)
IV2 <- tryCatch(feols(u_rate ~ log_gdp + urban_pct + food_import_pct + undernourish_pct +
                        inflation_pct + democracy | iso3 + yearf |
                        food_infl + FPG ~ Zg + ZgG, data = iv2, cluster = ~ iso3),
                error = function(e) {cat("err:", conditionMessage(e), "\n"); NULL})
if (!is.null(IV2)) { print(summary(IV2, stage = 2))
  print(fitstat(IV2, type = c("ivf","ivwald"))) }

cat("\n\n===== RQ1-v2: contemporaneous undernourishment =====\n")
m <- m %>% mutate(UND_c = undernourish_pct - mean(undernourish_pct, na.rm = TRUE))
M_RQ1b <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 + food_infl * UND_c +
                    log_gdp + urban_pct + food_import_pct + inflation_pct + democracy +
                    iso3 + yearf + offset(logP), family = nbinom2, data = m)
s_rq1b <- summary(M_RQ1b)$coefficients$cond
undc_int <- grep("food_infl:UND_c|UND_c:food_infl", rownames(s_rq1b), value = TRUE)
cat("contemporaneous UND:\n")
print(round(s_rq1b[c("food_infl","food_infl2","G_c","food_infl:G_c","UND_c",undc_int), ], 5))
cat(sprintf("N = %d | AIC = %.1f\n", nobs(M_RQ1b), AIC(M_RQ1b)))

cat("\n--- compare RQ1 variants (food_infl:UND interaction) ---\n")
cat(sprintf("  lagged UND   : b = %.5f | p = %.4f (AIC %.1f)\n",
            s_rq1["food_infl:UND_l1_c","Estimate"], s_rq1["food_infl:UND_l1_c","Pr(>|z|)"], AIC(M_RQ1)))
cat(sprintf("  contemp UND  : b = %.5f | p = %.4f (AIC %.1f)\n",
            s_rq1b[undc_int,"Estimate"], s_rq1b[undc_int,"Pr(>|z|)"], AIC(M_RQ1b)))

cat("\n===== DONE =====\n")









cat("===== RQ1-v3: triple interaction food_infl x G_c x UND =====\n")
m <- m %>%
  arrange(iso3, Year) %>% group_by(iso3) %>%
  mutate(undernourish_l1 = dplyr::lag(undernourish_pct)) %>% ungroup() %>%
  mutate(UND_l1_c = undernourish_l1 - mean(undernourish_l1, na.rm = TRUE),
         UND_c    = undernourish_pct - mean(undernourish_pct, na.rm = TRUE))

M_RQ1c <- glmmTMB(unrest_n ~ food_infl * G_c * UND_c + food_infl2 +
                    log_gdp + urban_pct + food_import_pct + inflation_pct + democracy +
                    iso3 + yearf + offset(logP), family = nbinom2, data = m)
s_c <- summary(M_RQ1c)$coefficients$cond
keep_c <- rownames(s_c)[grepl("food_infl|G_c|UND_c", rownames(s_c)) & !grepl("^iso3|^yearf", rownames(s_c))]
cat("triple-interaction model terms:\n")
print(round(s_c[keep_c, ], 5))
cat(sprintf("N = %d | AIC = %.1f\n", nobs(M_RQ1c), AIC(M_RQ1c)))
triple_term <- grep("food_infl:G_c:UND_c|food_infl:UND_c:G_c|G_c:food_infl:UND_c|UND_c:G_c:food_infl", rownames(s_c), value = TRUE)
cat("triple term:", triple_term, "| p =", round(s_c[triple_term,"Pr(>|z|)"],4), "\n")

cat("\n===== RQ1-v4: high/low undernourishment split =====\n")
und_med <- median(m$undernourish_pct, na.rm = TRUE)
m <- m %>% mutate(UND_high = ifelse(undernourish_pct >= und_med, 1L, 0L))
cat("undernourishment median:", round(und_med,1), "| high group n:", sum(m$UND_high==1, na.rm=TRUE), "\n")

M_RQ1d <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 + food_infl * UND_high +
                    log_gdp + urban_pct + food_import_pct + inflation_pct + democracy +
                    iso3 + yearf + offset(logP), family = nbinom2, data = m)
s_d <- summary(M_RQ1d)$coefficients$cond
undh_int <- grep("food_infl:UND_high|UND_high:food_infl", rownames(s_d), value = TRUE)
print(round(s_d[c("food_infl","food_infl2","G_c","food_infl:G_c","UND_high",undh_int), ], 5))

cat("\n--- ALL RQ1 variants summary (food_infl x UND interaction) ---\n")
cat(sprintf("  v1 lagged UND continuous : b = %+.5f | p = %.4f\n",
            s_rq1["food_infl:UND_l1_c","Estimate"], s_rq1["food_infl:UND_l1_c","Pr(>|z|)"]))
cat(sprintf("  v2 contemp UND continuous: b = %+.5f | p = %.4f\n",
            summary(M_RQ1b)$coefficients$cond[grep("food_infl:UND_c|UND_c:food_infl", rownames(summary(M_RQ1b)$coefficients$cond), value=TRUE),"Estimate"],
            summary(M_RQ1b)$coefficients$cond[grep("food_infl:UND_c|UND_c:food_infl", rownames(summary(M_RQ1b)$coefficients$cond), value=TRUE),"Pr(>|z|)"]))
cat(sprintf("  v3 triple interaction    : b = %+.5f | p = %.4f\n",
            s_c[triple_term,"Estimate"], s_c[triple_term,"Pr(>|z|)"]))
cat(sprintf("  v4 high/low UND split    : b = %+.5f | p = %.4f\n",
            s_d[undh_int,"Estimate"], s_d[undh_int,"Pr(>|z|)"]))
cat("\n  (food_infl:G_c stays significant throughout? check below)\n")
cat(sprintf("  food_infl:G_c in v3 = %+.5f (p=%.4f) | v4 = %+.5f (p=%.4f)\n",
            s_c["food_infl:G_c","Estimate"], s_c["food_infl:G_c","Pr(>|z|)"],
            s_d["food_infl:G_c","Estimate"], s_d["food_infl:G_c","Pr(>|z|)"]))

cat("\n===== DONE =====\n")


library(splines)

cat("===== RQ1-v5a: undernourishment terciles x food_infl =====\n")
und_t <- quantile(m$undernourish_pct, c(1/3, 2/3), na.rm = TRUE)
m <- m %>% mutate(
  UND_grp = cut(undernourish_pct, breaks = c(-Inf, und_t[1], und_t[2], Inf),
                labels = c("low","mid","high"))
)
cat("tercile cutoffs:", round(as.numeric(und_t),1), "| group sizes:\n")
print(table(m$UND_grp))

M_RQ1e <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 + food_infl * UND_grp +
                    log_gdp + urban_pct + food_import_pct + inflation_pct + democracy +
                    iso3 + yearf + offset(logP), family = nbinom2, data = m)
s_e <- summary(M_RQ1e)$coefficients$cond
ke <- rownames(s_e)[grepl("food_infl|UND_grp|G_c", rownames(s_e)) & !grepl("^iso3|^yearf", rownames(s_e))]
print(round(s_e[ke, ], 5))
cat(sprintf("N = %d | AIC = %.1f\n", nobs(M_RQ1e), AIC(M_RQ1e)))

cat("\n===== RQ1-v5b: natural spline interaction food_infl x ns(UND,3) =====\n")
m_sp <- m %>% filter(!is.na(undernourish_pct), !is.na(food_infl)) %>% droplevels()
M_RQ1f <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 +
                    food_infl * ns(undernourish_pct, 3) +
                    log_gdp + urban_pct + food_import_pct + inflation_pct + democracy +
                    iso3 + yearf + offset(logP), family = nbinom2, data = m_sp)
cat(sprintf("spline model: N = %d | AIC = %.1f\n", nobs(M_RQ1f), AIC(M_RQ1f)))

cat("\n--- LR test: does food_infl x UND interaction (spline) add anything? ---\n")
M_base_sp <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 +
                       ns(undernourish_pct, 3) +
                       log_gdp + urban_pct + food_import_pct + inflation_pct + democracy +
                       iso3 + yearf + offset(logP), family = nbinom2, data = m_sp)
print(anova(M_base_sp, M_RQ1f))

cat("\n--- LR test: does the tercile interaction add anything? ---\n")
M_base_grp <- glmmTMB(unrest_n ~ food_infl * G_c + food_infl2 + UND_grp +
                        log_gdp + urban_pct + food_import_pct + inflation_pct + democracy +
                        iso3 + yearf + offset(logP), family = nbinom2, data = m)
print(anova(M_base_grp, M_RQ1e))

cat("\n--- food_infl:G_c (main thesis) still significant in v5? ---\n")
cat(sprintf("  v5a tercile: food_infl:G_c = %+.5f (p=%.4f)\n",
            s_e["food_infl:G_c","Estimate"], s_e["food_infl:G_c","Pr(>|z|)"]))
sf <- summary(M_RQ1f)$coefficients$cond
cat(sprintf("  v5b spline : food_infl:G_c = %+.5f (p=%.4f)\n",
            sf["food_infl:G_c","Estimate"], sf["food_infl:G_c","Pr(>|z|)"]))

cat("\n===== DONE =====\n")

# kayit

DATA_DIR <- "data"
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

save_path <- file.path(DATA_DIR, "ssa_analysis_workspace.RData")
save.image(file = save_path)

cat("Tum workspace kaydedildi:\n", save_path, "\n\n")
cat("Dosya boyutu:", round(file.size(save_path)/1024^2, 1), "MB\n\n")

key_objects <- c("m", "M1", "M2", "M3", "M_RQ1", "M_RQ1b", "M_RQ1c", "M_RQ1d",
                 "M_RQ1e", "M_RQ1f", "M_gov", "M_oop", "M_gov_oop", "PPML_M2",
                 "CF2", "IV2", "FS2", "iv2", "ME_gini", "ME_und", "infl_tab",
                 "thr", "gini_mean", "fp_q", "g_q")
present <- key_objects[sapply(key_objects, exists)]
missing <- setdiff(key_objects, present)
cat("Kaydedilen kilit nesneler (", length(present), "):\n", sep="")
cat(paste(present, collapse=", "), "\n")
if (length(missing)) cat("\nBellekte olmayan (sorun degil):", paste(missing, collapse=", "), "\n")

DATA_DIR <- "data"
b <- fixef(M2)$cond
V <- vcov(M2)$cond
keep <- c("(Intercept)","food_infl","food_infl2","G_c","food_infl:G_c",
          "log_gdp","urban_pct","food_import_pct","undernourish_pct","inflation_pct","democracy")
write.csv(data.frame(term=keep, coef=b[keep]), file.path(DATA_DIR,"m2_coef.csv"), row.names=FALSE)
write.csv(as.data.frame(V[keep,keep]), file.path(DATA_DIR,"m2_vcov.csv"), row.names=TRUE)
ctrl <- with(m, data.frame(
  log_gdp=mean(log_gdp,na.rm=TRUE), urban_pct=mean(urban_pct,na.rm=TRUE),
  food_import_pct=mean(food_import_pct,na.rm=TRUE), undernourish_pct=mean(undernourish_pct,na.rm=TRUE),
  inflation_pct=mean(inflation_pct,na.rm=TRUE), democracy=mean(democracy,na.rm=TRUE),
  gini_mean=mean(gini,na.rm=TRUE), fp_median=median(food_infl,na.rm=TRUE),
  gini_p05=quantile(gini,.05,na.rm=TRUE), gini_p95=quantile(gini,.95,na.rm=TRUE),
  gini_p10=quantile(gini,.10,na.rm=TRUE), gini_p50=quantile(gini,.50,na.rm=TRUE),
  gini_p90=quantile(gini,.90,na.rm=TRUE)))
write.csv(ctrl, file.path(DATA_DIR,"m2_meta.csv"), row.names=FALSE)
cat("yazildi: m2_coef.csv, m2_vcov.csv, m2_meta.csv\n")


DATA_DIR <- "data"

get_int <- function(mod, term) {
  s <- summary(mod)$coefficients$cond
  data.frame(est = s[term,"Estimate"], se = s[term,"Std. Error"], p = s[term,"Pr(>|z|)"])
}

specs <- list(
  list(lab="Central (non-linear FE-NB)",        mod=M2,        term="food_infl:G_c"),
  list(lab="+ Governance (gov. effectiveness)", mod=M_gov,     term="food_infl:G_c"),
  list(lab="+ Health financing (OOP)",          mod=M_oop,     term="food_infl:G_c"),
  list(lab="+ Governance + OOP",                mod=M_gov_oop, term="food_infl:G_c"),
  list(lab="Winsorised 1-99%",                  mod=m2b,       term="food_infl_w:G_c"),
  list(lab="Exclude inflation > 50%",           mod=m2n,       term="food_infl:G_c"),
  list(lab="Exclude inflation > 200%",          mod=m2_no,     term="food_infl:G_c"),
  list(lab="Exclude Zimbabwe + S. Sudan",       mod=m2q_noZS,  term="food_infl:G_c")
)

fp <- do.call(rbind, lapply(specs, function(s){
  v <- get_int(s$mod, s$term)
  data.frame(spec=s$lab, est=v$est, se=v$se, p=v$p,
             lo=v$est-1.96*v$se, hi=v$est+1.96*v$se)
}))
fp$irr <- exp(fp$est); fp$irr_lo <- exp(fp$lo); fp$irr_hi <- exp(fp$hi)
write.csv(fp, file.path(DATA_DIR,"fig2_forest.csv"), row.names=FALSE)
print(fp[,c("spec","est","p","irr","irr_lo","irr_hi")])
cat("\nyazildi: fig2_forest.csv\n")

DATA_DIR <- "data"
meta <- read.csv(file.path(DATA_DIR,"m2_meta.csv"))
meta$model       <- "Two-way FE negative binomial"
meta$n_obs       <- nobs(M2)
meta$n_countries <- length(unique(m$iso3[!is.na(m$food_infl)]))
meta$n_years     <- length(unique(m$Year))
write.csv(meta, file.path(DATA_DIR,"m2_meta.csv"), row.names=FALSE)
cat("meta guncellendi: N =", nobs(M2),
    "| countries =", meta$n_countries, "| years =", meta$n_years, "\n")


DATA_DIR <- "data"

mvars <- c("unrest_n","food_infl","food_infl2","G_c","log_gdp","urban_pct",
           "food_import_pct","undernourish_pct","inflation_pct","democracy","logP")
mcc <- m[complete.cases(m[, mvars]), ]
n_ctry_model <- length(unique(mcc$iso3))
n_year_model <- length(unique(mcc$Year))

cat("M2 efektif: N =", nrow(mcc),
    "| ulke =", n_ctry_model, "| yil =", n_year_model, "\n")
cat("dropped countries (panelde olup modelde olmayan):\n")
print(setdiff(unique(m$iso3), unique(mcc$iso3)))

meta <- read.csv(file.path(DATA_DIR,"m2_meta.csv"))
meta$model       <- "Two-way FE negative binomial"
meta$n_obs       <- nobs(M2)
meta$n_countries <- n_ctry_model
meta$n_years     <- n_year_model
write.csv(meta, file.path(DATA_DIR,"m2_meta.csv"), row.names=FALSE)
cat("\nmeta guncellendi: N =", meta$n_obs, "| countries =", meta$n_countries, "\n")

# ek bilgiler

options(scipen=999)
library(plm)

cat("===== TABLE 2 PANEL A: Pesaran CD (H0: cross-sectional independence) =====\n")
cd_vars <- c("unrest_n","food_infl","gini")
for (v in cd_vars) {
  d <- m[!is.na(m[[v]]), c("iso3","Year",v)]
  pd <- pdata.frame(d, index=c("iso3","Year"))
  f <- as.formula(paste(v, "~ 1"))
  r <- tryCatch(pcdtest(f, data=pd, test="cd"),
                error=function(e){cat(v,"ERR:",conditionMessage(e),"\n");NULL})
  if(!is.null(r)) cat(sprintf("  %-16s CD = %8.3f | p = %.4g | %s\n", v,
                              as.numeric(r$statistic), r$p.value,
                              ifelse(r$p.value<0.05,"reject H0 (dependence)","fail to reject")))
}

cat("\n===== TABLE 2 PANEL B: CIPS unit root (H0: unit root / non-stationary) =====\n")
cips_one <- function(v, drop_iso, start_year, diff=FALSE) {
  d <- m[!is.na(m[[v]]) & !(m$iso3 %in% drop_iso), c("iso3","Year",v)]
  if(!is.null(start_year)) d <- d[d$Year>=start_year,]
  d <- d[order(d$iso3,d$Year),]
  pd <- pdata.frame(d, index=c("iso3","Year"))
  x <- pd[[v]]
  if(diff) x <- diff(x)
  r <- tryCatch(cipstest(x, type="drift", model="cmg", lags=2),
                error=function(e){cat(v,"ERR:",conditionMessage(e),"\n");NULL})
  if(!is.null(r)) cat(sprintf("  %-22s CIPS = %7.3f | p %s %.3f | %s\n",
                              paste0(v, if(diff)" (diff)" else " (levels)"),
                              as.numeric(r$statistic), ifelse(r$p.value<=0.01,"<","="),
                              ifelse(r$p.value<=0.01,0.01,r$p.value),
                              ifelse(r$p.value<0.05,"reject H0 (stationary)","fail to reject (unit root)")))
}
cips_one("unrest_n", character(0), NULL, FALSE)
cips_one("food_infl", c("ERI","CAF"), 2001, FALSE)
cips_one("gini", c("ERI","SOM"), NULL, FALSE)
cips_one("gini", c("ERI","SOM"), NULL, TRUE)

cat("\n===== TABLE 2 PANEL C: Dumitrescu-Hurlin (H0: homogeneous non-causality) =====\n")
dh_one <- function(yv, xv, lag, drop_iso, start_year, lab) {
  d <- m[!is.na(m[[yv]]) & !is.na(m[[xv]]) & !(m$iso3 %in% drop_iso), c("iso3","Year",yv,xv)]
  if(!is.null(start_year)) d <- d[d$Year>=start_year,]
  d <- d[order(d$iso3,d$Year),]
  ok <- names(which(table(d$iso3) >= (2*lag+3)))
  d <- d[d$iso3 %in% ok,]
  pd <- pdata.frame(d, index=c("iso3","Year"))
  f <- as.formula(paste(yv,"~",xv))
  r <- tryCatch(pgrangertest(f, data=pd, order=lag, test="Ztilde"),
                error=function(e){cat("ERR:",conditionMessage(e),"\n");NULL})
  if(!is.null(r)) cat(sprintf("  %-26s lag %d: Ztilde = %7.3f | p = %.4g | %s\n",
                              lab, lag, as.numeric(r$statistic), r$p.value,
                              ifelse(r$p.value<0.05,"reject H0","fail to reject")))
}
# d_gini for H7
m$d_gini <- ave(m$gini, m$iso3, FUN=function(x) c(NA,diff(x)))
for(L in 1:3) dh_one("unrest_n","food_infl",L,c("ERI","CAF"),2001,"food_infl =/=> unrest")
for(L in 1:3) dh_one("food_infl","unrest_n",L,c("ERI","CAF"),2001,"unrest =/=> food_infl")
for(L in 1:2) dh_one("unrest_n","d_gini",L,c("ERI","SOM"),NULL,"d_gini =/=> unrest")

cat("\n===== DONE TABLE 2 =====\n")












options(scipen=999)
library(plm)
m$u_rate <- log1p((m$unrest_n / m$population) * 1e6)

cips_test <- function(v, drop_iso=character(0), start_year=NULL, diff=FALSE, type="drift", lags=2) {
  d <- m[!is.na(m[[v]]) & !(m$iso3 %in% drop_iso), c("iso3","Year",v)]
  if(!is.null(start_year)) d <- d[d$Year>=start_year,]
  d <- d[order(d$iso3,d$Year),]
  pd <- pdata.frame(d, index=c("iso3","Year"))
  x <- pd[[v]]; if(diff) x <- diff(x)
  r <- tryCatch(cipstest(x, type=type, model="cmg", lags=lags),
                error=function(e)NULL)
  if(!is.null(r)) cat(sprintf("  %-20s type=%-5s lag=%d: CIPS=%7.3f | p=%s%.3f | %s\n",
                              paste0(v, if(diff)"(d)" else ""), type, lags, as.numeric(r$statistic),
                              ifelse(r$p.value<=0.01,"<",ifelse(r$p.value>=0.10,">","=")),
                              ifelse(r$p.value<=0.01,0.01,ifelse(r$p.value>=0.10,0.10,r$p.value)),
                              ifelse(r$p.value<0.05,"stationary","unit root")))
}

cat("=== CIPS on u_rate (log event rate) ===\n")
cips_test("u_rate", type="drift", lags=1)
cips_test("u_rate", type="drift", lags=2)
cips_test("u_rate", type="trend", lags=1)
cips_test("u_rate", type="trend", lags=2)

cat("\n=== food_infl (confirm robust) ===\n")
cips_test("food_infl", c("ERI","CAF"), 2001, type="drift", lags=1)
cips_test("food_infl", c("ERI","CAF"), 2001, type="trend", lags=1)

cat("\n=== gini levels vs diff ===\n")
cips_test("gini", c("ERI","SOM"), type="drift", lags=1)
cips_test("gini", c("ERI","SOM"), diff=TRUE, type="drift", lags=1)



options(scipen=999)

cat("===== DESCRIPTIVES =====\n")
dvars <- c("unrest_n","food_infl","gini","undernourish_pct","gdp_pc","log_gdp",
           "urban_pct","food_import_pct","inflation_pct","democracy","gov_effect","oop_health_pct")
desc <- do.call(rbind, lapply(dvars, function(v){
  x <- m[[v]]
  data.frame(var=v, N=sum(!is.na(x)),
             mean=round(mean(x,na.rm=TRUE),2), sd=round(sd(x,na.rm=TRUE),2),
             med=round(median(x,na.rm=TRUE),2),
             q1=round(quantile(x,.25,na.rm=TRUE),2), q3=round(quantile(x,.75,na.rm=TRUE),2),
             min=round(min(x,na.rm=TRUE),2), max=round(max(x,na.rm=TRUE),2))
}))
print(desc, row.names=FALSE)

cat("\n===== TABLE 3: ALL coefficients M1, M2, M3 =====\n")
for (nm in c("M1","M2","M3")) {
  s <- summary(get(nm))$coefficients$cond
  keep <- !grepl("^iso3|^yearf|Intercept", rownames(s))
  cat("\n---", nm, "(N =", nobs(get(nm)), ", theta =", round(sigma(get(nm)),3),
      ", AIC =", round(AIC(get(nm)),1), ") ---\n")
  print(round(s[keep,], 5))
}

cat("\n===== TABLE 4: control-function (CF2) full =====\n")
s_cf <- summary(CF2)$coefficients$cond
print(round(s_cf[!grepl("^iso3|^yearf|Intercept", rownames(s_cf)),], 5))










options(scipen=999)

cat("===== IV2 / 2SLS benchmark gercek degerler =====\n")
print(class(IV2))
cat("\n--- summary ---\n")
print(summary(IV2))

cat("\n===== sadece katsayilar (food_infl ve interaction) =====\n")
co <- coef(summary(IV2))
print(round(co, 5))

cat("\n===== first-stage F (KP) =====\n")
if (!is.null(IV2$iv_first_stage)) print(IV2$iv_first_stage)
cat("\n--- fitstat KP ---\n")
tryCatch(print(fitstat(IV2, ~ ivf1 + ivwald)), error=function(e) cat("fitstat:", conditionMessage(e), "\n"))

cat("\n===== Wu-Hausman / endogeneity =====\n")
tryCatch(print(fitstat(IV2, ~ wh)), error=function(e) cat("wh:", conditionMessage(e), "\n"))

options(scipen=999)

get_int <- function(mod, term="food_infl:G_c") {
  s <- tryCatch(summary(mod)$coefficients$cond, error=function(e) NULL)
  if (is.null(s) || !(term %in% rownames(s))) return(c(NA,NA,NA))
  s[term, c("Estimate","Std. Error","Pr(>|z|)")]
}
report <- function(lab, mod, term="food_infl:G_c") {
  v <- get_int(mod, term)
  cat(sprintf("  %-38s b=%9.5f  SE=%8.5f  p=%.4f  %s\n",
              lab, v[1], v[2], v[3], ifelse(!is.na(v[3]) && v[3]<0.05,"reject H0","--")))
}

m$logP <- log(m$population)
m$G_c  <- m$gini - mean(m$gini, na.rm=TRUE)

cat("===== BASELINE (teyit) =====\n")
report("Central M2", M2)

cat("\n===== 1. DETRENDED FOOD PRICE SHOCK =====\n")
# food_shock = detrended; onun karesi ve G_c etkilesimi
m$shock_c   <- m$food_shock
m$shock2    <- m$food_shock^2
m$shock_Gc  <- m$food_shock * m$G_c
M_detr <- tryCatch(glmmTMB(unrest_n ~ food_shock + shock2 + G_c + shock_Gc +
                             log_gdp + urban_pct + food_import_pct + undernourish_pct +
                             inflation_pct + democracy + (1|iso3) + (1|yearf),
                           offset=logP, family=nbinom2, data=m),
                   error=function(e){cat("  ERR:",conditionMessage(e),"\n");NULL})
if(!is.null(M_detr)) report("Detrended shock x inequality", M_detr, "shock_Gc")

cat("\n===== 2. OBSERVED INEQUALITY ONLY (gini_wb ham) =====\n")
m$Gobs_c    <- m$gini_wb - mean(m$gini_wb, na.rm=TRUE)
m$fp_Gobs   <- m$food_infl * m$Gobs_c
M_obs <- tryCatch(glmmTMB(unrest_n ~ food_infl + food_infl2 + Gobs_c + fp_Gobs +
                            log_gdp + urban_pct + food_import_pct + undernourish_pct +
                            inflation_pct + democracy + (1|iso3) + (1|yearf),
                          offset=logP, family=nbinom2, data=m),
                  error=function(e){cat("  ERR:",conditionMessage(e),"\n");NULL})
if(!is.null(M_obs)) report("Observed Gini (WB) x food infl", M_obs, "fp_Gobs")

cat("\n===== 3. EXCLUDING GENERAL INFLATION =====\n")
M_noinf <- tryCatch(glmmTMB(unrest_n ~ food_infl + food_infl2 + G_c + food_infl:G_c +
                              log_gdp + urban_pct + food_import_pct + undernourish_pct +
                              democracy + (1|iso3) + (1|yearf),
                            offset=logP, family=nbinom2, data=m),
                    error=function(e){cat("  ERR:",conditionMessage(e),"\n");NULL})
if(!is.null(M_noinf)) report("Excl. general inflation", M_noinf)

cat("\n===== 4. PANEL START-YEAR SENSITIVITY (>=2005) =====\n")
m05 <- m[m$Year >= 2005, ]
M_y05 <- tryCatch(glmmTMB(unrest_n ~ food_infl + food_infl2 + G_c + food_infl:G_c +
                            log_gdp + urban_pct + food_import_pct + undernourish_pct +
                            inflation_pct + democracy + (1|iso3) + (1|yearf),
                          offset=logP, family=nbinom2, data=m05),
                  error=function(e){cat("  ERR:",conditionMessage(e),"\n");NULL})
if(!is.null(M_y05)) report("Start year 2005", M_y05)

cat("\n===== kontrol: food_shock ve gini_wb mevcut mu, kac NA =====\n")
cat("food_shock NA:", sum(is.na(m$food_shock)), "/ var:", "food_shock" %in% names(m), "\n")
cat("gini_wb    NA:", sum(is.na(m$gini_wb)),    "/ var:", "gini_wb" %in% names(m), "\n")
cat("food_shock ozet:\n"); print(summary(m$food_shock))














options(scipen=999)
cat("===== HAM SWIID kontrol: hangi degiskenler var =====\n")
print(grep("gini", names(m), value=TRUE, ignore.case=TRUE))
cat("\nHer gini degiskeninin NA sayisi ve ozeti:\n")
for(v in grep("gini", names(m), value=TRUE, ignore.case=TRUE)) {
  cat(sprintf("  %-14s NA=%4d  range=[%.1f, %.1f]\n", v, sum(is.na(m[[v]])),
              min(m[[v]],na.rm=TRUE), max(m[[v]],na.rm=TRUE)))
}

cat("\n===== gini (kullanilan) interpole mi, kac essiz ardisik tekrar =====\n")
# eger interpolasyon varsa ayni deger ardisik tekrarlanir
chk <- m[order(m$iso3, m$Year), c("iso3","Year","gini")]
chk$same_as_prev <- c(FALSE, chk$gini[-1] == chk$gini[-nrow(chk)] &
                        chk$iso3[-1] == chk$iso3[-nrow(chk)])
cat("Ardisik ayni gini orani (interpolasyon isareti):",
    round(mean(chk$same_as_prev, na.rm=TRUE),3), "\n")

cat("\n===== detrended shock olcek: food_infl ile karsilastir =====\n")
cat("food_infl  SD:", round(sd(m$food_infl,na.rm=TRUE),3), "\n")
cat("food_shock SD:", round(sd(m$food_shock,na.rm=TRUE),3), "\n")
cat("=> detrended katsayisi food_infl'inkinden ~", 
    round(sd(m$food_infl,na.rm=TRUE)/sd(m$food_shock,na.rm=TRUE),0), "kat buyuk beklenir\n")





















options(scipen=999)
library(glmmTMB)

cat("===== food_infl vs inflation_pct korelasyon =====\n")
cat("Pearson r:", round(cor(m$food_infl, m$inflation_pct, use="complete.obs"),3), "\n")
cat("VIF benzeri: ikisi ayni modelde mi sorun?\n")

cat("\n===== excl-inflation modelinde food_infl ana etki ne oluyor =====\n")
m$logP <- log(m$population); m$G_c <- m$gini - mean(m$gini,na.rm=TRUE)
M_noinf <- glmmTMB(unrest_n ~ food_infl + food_infl2 + G_c + food_infl:G_c +
                     log_gdp + urban_pct + food_import_pct + undernourish_pct +
                     democracy + (1|iso3) + (1|yearf),
                   offset=logP, family=nbinom2, data=m)
cat("--- M_noinf tum food/gini terimleri ---\n")
s <- summary(M_noinf)$coefficients$cond
print(round(s[grep("food|G_c", rownames(s)), ], 5))

cat("\n===== KARSILASTIRMA: baseline M2 ayni terimler =====\n")
s2 <- summary(M2)$coefficients$cond
print(round(s2[grep("food|G_c", rownames(s2)), ], 5))

cat("\n===== ALTERNATIF: food-spesifik enflasyon (food_infl - genel enflasyon farki) =====\n")
# gida enflasyonunun genel enflasyondan ARINMIS kismi
m$food_excess <- m$food_infl - m$inflation_pct
m$food_excess2 <- m$food_excess^2
m$fe_Gc <- m$food_excess * m$G_c
M_excess <- tryCatch(glmmTMB(unrest_n ~ food_excess + food_excess2 + G_c + fe_Gc +
                               log_gdp + urban_pct + food_import_pct + undernourish_pct +
                               inflation_pct + democracy + (1|iso3) + (1|yearf),
                             offset=logP, family=nbinom2, data=m),
                     error=function(e){cat("ERR:",conditionMessage(e),"\n");NULL})
if(!is.null(M_excess)){
  cat("--- relative food inflation (food - general) x inequality ---\n")
  se <- summary(M_excess)$coefficients$cond
  print(round(se[grep("food|fe_Gc|G_c", rownames(se)), ], 5))
}

cat("\n===== food_infl ve inflation_pct ikisinin de interaction'i =====\n")
m$infl_Gc <- m$inflation_pct * m$G_c
M_both <- tryCatch(glmmTMB(unrest_n ~ food_infl + food_infl2 + G_c + food_infl:G_c +
                             inflation_pct + infl_Gc +
                             log_gdp + urban_pct + food_import_pct + undernourish_pct +
                             democracy + (1|iso3) + (1|yearf),
                           offset=logP, family=nbinom2, data=m),
                   error=function(e){cat("ERR:",conditionMessage(e),"\n");NULL})
if(!is.null(M_both)){
  cat("--- food x gini VE general-infl x gini ayni modelde ---\n")
  sb <- summary(M_both)$coefficients$cond
  print(round(sb[grep("food|infl_Gc|G_c", rownames(sb)), ], 5))
}















options(scipen=999)
library(glmmTMB)
m$logP <- log(m$population); m$G_c <- m$gini - mean(m$gini,na.rm=TRUE)

cat("===== 1. food_infl ve inflation_pct VIF (ayni modelde) =====\n")
library(car)
lm_chk <- lm(unrest_n ~ food_infl + inflation_pct + G_c + log_gdp + urban_pct +
               food_import_pct + undernourish_pct + democracy, data=m)
print(round(vif(lm_chk),2))

cat("\n===== 2. food CPI seviyesi vs genel CPI - bunlar gercekten farkli mi? =====\n")
cat("food_infl ile inflation_pct ozet farki:\n")
cat("  food_infl  mean/median:", round(mean(m$food_infl,na.rm=TRUE),2), "/", round(median(m$food_infl,na.rm=TRUE),2),"\n")
cat("  inflation  mean/median:", round(mean(m$inflation_pct,na.rm=TRUE),2), "/", round(median(m$inflation_pct,na.rm=TRUE),2),"\n")
cat("  fark (food - genel) mean/median:", round(mean(m$food_infl-m$inflation_pct,na.rm=TRUE),2),"/", round(median(m$food_infl-m$inflation_pct,na.rm=TRUE),2),"\n")
cat("  korelasyon zaten r=0.985\n")

cat("\n===== 3. SADECE genel enflasyon interaction (food yok) =====\n")
m$infl_Gc <- m$inflation_pct * m$G_c; m$inflation2 <- m$inflation_pct^2
M_genonly <- glmmTMB(unrest_n ~ inflation_pct + inflation2 + G_c + infl_Gc +
                       log_gdp + urban_pct + food_import_pct + undernourish_pct +
                       democracy + (1|iso3) + (1|yearf),
                     offset=logP, family=nbinom2, data=m)
cat("--- general inflation x gini (food modelde YOK) ---\n")
sg <- summary(M_genonly)$coefficients$cond
print(round(sg[grep("infl|G_c", rownames(sg)),],5))

cat("\n===== 4. food modelde inflation_pct HIC yokken (M1 benzeri ama temiz) =====\n")
cat("(bu zaten M2 ama inflation_pct ana etki olmadan)\n")
M_foodclean <- glmmTMB(unrest_n ~ food_infl + food_infl2 + G_c + food_infl:G_c +
                         log_gdp + urban_pct + food_import_pct + undernourish_pct +
                         democracy + (1|iso3) + (1|yearf),
                       offset=logP, family=nbinom2, data=m)
cat("--- food x gini, inflation_pct modelde yok ---\n")
sf <- summary(M_foodclean)$coefficients$cond
print(round(sf[grep("food|G_c", rownames(sf)),],5))
cat("\n(Not: bu M_noinf ile ayni - food_infl:G_c = 0.00108 p=0.16)\n")

cat("\n===== 5. KRITIK: food_infl, inflation_pct'ten ARINDIRILINCA (orthogonalize) =====\n")
# food_infl'i inflation_pct uzerine regres et, residual = gida-spesifik kisim
m_cc <- m[!is.na(m$food_infl) & !is.na(m$inflation_pct),]
res_mod <- lm(food_infl ~ inflation_pct, data=m_cc)
m_cc$food_orth <- residuals(res_mod)
m_cc$food_orth2 <- m_cc$food_orth^2
m_cc$fo_Gc <- m_cc$food_orth * m_cc$G_c
M_orth <- glmmTMB(unrest_n ~ food_orth + food_orth2 + G_c + fo_Gc + inflation_pct +
                    log_gdp + urban_pct + food_import_pct + undernourish_pct +
                    democracy + (1|iso3) + (1|yearf),
                  offset=logP, family=nbinom2, data=m_cc)
cat("--- orthogonalized food (genel enflasyondan arinmis) x gini ---\n")
so <- summary(M_orth)$coefficients$cond
print(round(so[grep("food_orth|fo_Gc|G_c|inflation", rownames(so)),],5))


options(scipen=999)
m$logP <- log(m$population); m$G_c <- m$gini - mean(m$gini,na.rm=TRUE)

cat("===== DETRENDED SHOCK - standardize (z-skor) =====\n")
m$fs_z   <- as.numeric(scale(m$food_shock))
m$fs_z2  <- m$fs_z^2
m$fsz_Gc <- m$fs_z * m$G_c
M_detr_z <- glmmTMB(unrest_n ~ fs_z + fs_z2 + G_c + fsz_Gc +
                      log_gdp + urban_pct + food_import_pct + undernourish_pct +
                      inflation_pct + democracy + (1|iso3) + (1|yearf),
                    offset=logP, family=nbinom2, data=m)
s <- summary(M_detr_z)$coefficients$cond
cat("Detrended shock (standardized) x inequality:\n")
print(round(s["fsz_Gc", c("Estimate","Std. Error","Pr(>|z|)")],5))

cat("\n===== START-YEAR 2005 (teyit, tam deger) =====\n")
m05 <- m[m$Year >= 2005, ]
M_y05 <- glmmTMB(unrest_n ~ food_infl + food_infl2 + G_c + food_infl:G_c +
                   log_gdp + urban_pct + food_import_pct + undernourish_pct +
                   inflation_pct + democracy + (1|iso3) + (1|yearf),
                 offset=logP, family=nbinom2, data=m05)
s2 <- summary(M_y05)$coefficients$cond
cat("Start 2005, food_infl:G_c, N =", nobs(M_y05), ":\n")
print(round(s2["food_infl:G_c", c("Estimate","Std. Error","Pr(>|z|)")],5))

cat("\n===== GENEL ENFLASYON x GINI (Results paragrafi icin tam deger) =====\n")
m$infl_Gc <- m$inflation_pct * m$G_c; m$inflation2 <- m$inflation_pct^2
M_gen <- glmmTMB(unrest_n ~ inflation_pct + inflation2 + G_c + infl_Gc +
                   log_gdp + urban_pct + food_import_pct + undernourish_pct +
                   democracy + (1|iso3) + (1|yearf),
                 offset=logP, family=nbinom2, data=m)
s3 <- summary(M_gen)$coefficients$cond
cat("General inflation x inequality:\n")
print(round(s3["infl_Gc", c("Estimate","Std. Error","Pr(>|z|)")],5))


















options(scipen=999)
library(plm)
library(lmtest)
library(sandwich)

m$u_rate <- log1p((m$unrest_n / m$population) * 1e6)
m$G_c    <- m$gini - mean(m$gini, na.rm=TRUE)
m$fp_Gc  <- m$food_infl * m$G_c

dvars <- c("u_rate","food_infl","fp_Gc","G_c","log_gdp","urban_pct",
           "food_import_pct","undernourish_pct","inflation_pct","democracy")
d <- m[complete.cases(m[,dvars]), c("iso3","Year",dvars)]
d <- d[order(d$iso3, d$Year),]
pd <- pdata.frame(d, index=c("iso3","Year"))

fe <- plm(u_rate ~ food_infl + fp_Gc + G_c + log_gdp + urban_pct +
            food_import_pct + undernourish_pct + inflation_pct + democracy,
          data=pd, model="within", effect="twoways")

cat("===== LINEAR FE: interaction (fp_Gc) cesitli SE'lerle =====\n")
cat("--- (a) klasik SE ---\n")
print(round(coeftest(fe)["fp_Gc",],5))

cat("\n--- (b) cluster-robust (country) SE ---\n")
print(round(coeftest(fe, vcov=vcovHC(fe, type="HC1", cluster="group"))["fp_Gc",],5))

cat("\n--- (c) Driscoll-Kraay (CSD-robust) SE ---\n")
dk <- coeftest(fe, vcov=vcovSCC(fe, type="HC1", maxlag=2))
print(round(dk["fp_Gc",],5))

cat("\n===== tum DK katsayilari (interaction + ana terimler) =====\n")
print(round(coeftest(fe, vcov=vcovSCC(fe, type="HC1", maxlag=2))[c("food_infl","fp_Gc"),],5))

cat("\nN =", nrow(d), "| countries =", length(unique(d$iso3)), "| years =", length(unique(d$Year)), "\n")


options(scipen=999)
library(glmmTMB)
m$logP <- log(m$population); m$G_c <- m$gini - mean(m$gini,na.rm=TRUE)

cat("===== 1. NB modelinde yuksek-sayim ulkeleri cikar =====\n")
# en cok olayli ulkeler
ag <- aggregate(unrest_n ~ iso3, data=m, sum)
ag <- ag[order(-ag$unrest_n),]
cat("En cok olayli 6 ulke:\n"); print(head(ag,6))
top_ctry <- head(ag$iso3, 3)
cat("\nCikarilanlar:", paste(top_ctry, collapse=", "), "\n")
m_no <- m[!(m$iso3 %in% top_ctry),]
M_notop <- glmmTMB(unrest_n ~ food_infl + food_infl2 + G_c + food_infl:G_c +
                     log_gdp + urban_pct + food_import_pct + undernourish_pct +
                     inflation_pct + democracy + (1|iso3) + (1|yearf),
                   offset=logP, family=nbinom2, data=m_no)
s <- summary(M_notop)$coefficients$cond
cat("Interaction (top-3 cikarilmis), N =", nobs(M_notop), ":\n")
print(round(s["food_infl:G_c", c("Estimate","Std. Error","Pr(>|z|)")],5))

cat("\n===== 2. NB ana model marjinal etki - en yuksek Gini'de gercekten var mi =====\n")
# zaten tipping 57.5 idi; en eshitsiz ulkeler hangileri
hi_gini <- aggregate(gini ~ iso3, data=m, mean)
hi_gini <- hi_gini[order(-hi_gini$gini),]
cat("En eshitsiz 6 ulke (ortalama Gini):\n"); print(head(hi_gini,6))

cat("\n===== 3. Poisson (NB degil) ayni model - interaction =====\n")
M_pois <- glmmTMB(unrest_n ~ food_infl + food_infl2 + G_c + food_infl:G_c +
                    log_gdp + urban_pct + food_import_pct + undernourish_pct +
                    inflation_pct + democracy + (1|iso3) + (1|yearf),
                  offset=logP, family=poisson, data=m)
sp <- summary(M_pois)$coefficients$cond
cat("Poisson interaction:\n")
print(round(sp["food_infl:G_c", c("Estimate","Std. Error","Pr(>|z|)")],5))

cat("\n===== 4. NB ama log(unrest+1) lineer OLS karsilastirma (ayni outcome mantigii) =====\n")
m$log_unrest <- log1p(m$unrest_n)
# basit lm, FE dummy ile
M_loglin <- lm(log_unrest ~ food_infl + I(food_infl^2) + G_c + food_infl:G_c +
                 log_gdp + urban_pct + food_import_pct + undernourish_pct +
                 inflation_pct + democracy + factor(iso3) + factor(Year), data=m)
cl <- summary(M_loglin)$coefficients
cat("log(1+unrest) OLS interaction:\n")
print(round(cl["food_infl:G_c",],5))



options(scipen=999)
m$logP <- log(m$population); m$G_c <- m$gini - mean(m$gini,na.rm=TRUE)

base_fit <- function(dat) {
  glmmTMB(unrest_n ~ food_infl + food_infl2 + G_c + food_infl:G_c +
            log_gdp + urban_pct + food_import_pct + undernourish_pct +
            inflation_pct + democracy + (1|iso3) + (1|yearf),
          offset=logP, family=nbinom2, data=dat)
}

cat("===== 1. SADECE Guney Afrika cikar (en kritik test) =====\n")
M_noZAF <- base_fit(m[m$iso3 != "ZAF",])
s <- summary(M_noZAF)$coefficients$cond
v <- s["food_infl:G_c", c("Estimate","Std. Error","Pr(>|z|)")]
cat(sprintf("ZAF cikarilmis: b=%.5f SE=%.5f p=%.4f N=%d  %s\n",
            v[1],v[2],v[3],nobs(M_noZAF), ifelse(v[3]<0.05,"reject H0","COKTU")))

cat("\n===== 2. En esitsiz 3 ulke cikar (NAM, ZAF, SWZ) =====\n")
M_noHi <- base_fit(m[!(m$iso3 %in% c("NAM","ZAF","SWZ")),])
s <- summary(M_noHi)$coefficients$cond
v <- s["food_infl:G_c", c("Estimate","Std. Error","Pr(>|z|)")]
cat(sprintf("Hi-Gini-3 cikarilmis: b=%.5f SE=%.5f p=%.4f N=%d  %s\n",
            v[1],v[2],v[3],nobs(M_noHi), ifelse(v[3]<0.05,"reject H0","COKTU")))

cat("\n===== 3. LEAVE-ONE-COUNTRY-OUT JACKKNIFE (tum 40 ulke) =====\n")
ctrys <- unique(m$iso3[!is.na(m$food_infl)])
jack <- data.frame(iso3=character(), b=numeric(), p=numeric(), stringsAsFactors=FALSE)
for (cc in ctrys) {
  M <- tryCatch(base_fit(m[m$iso3 != cc,]), error=function(e) NULL)
  if (!is.null(M)) {
    s <- summary(M)$coefficients$cond
    if ("food_infl:G_c" %in% rownames(s)) {
      v <- s["food_infl:G_c", c("Estimate","Pr(>|z|)")]
      jack <- rbind(jack, data.frame(iso3=cc, b=v[1], p=v[2]))
    }
  }
}
jack <- jack[order(jack$p, decreasing=TRUE),]
cat("En kritik 8 ulke (cikarilinca p en cok artan = etkiyi en cok tasiyan):\n")
print(head(jack, 8), row.names=FALSE)
cat("\nOzet: jackknife p araligi [", round(min(jack$p),4), ",", round(max(jack$p),4), "]\n")
cat("Kac ulke cikariminda p>0.05 (etki kayboluyor):", sum(jack$p>0.05), "/", nrow(jack), "\n")
cat("Kac ulke cikariminda p<0.05 (etki saglam):", sum(jack$p<0.05), "/", nrow(jack), "\n")


options(scipen=999)
library(sf); library(rnaturalearth); library(rnaturalearthdata)
library(ggplot2); library(dplyr); library(scales); library(ggrepel)
sf_use_s2(FALSE)

library(ggrepel)
library(extrafont)
loadfonts(device = "win", quiet = TRUE)

font_import(pattern = "times", prompt = FALSE)
loadfonts(device = "win", quiet = TRUE)

gbar <- mean(m$gini, na.rm=TRUE)
m$G_c <- m$gini - gbar
M_full <- base_fit(m)
co  <- summary(M_full)$coefficients$cond
b1  <- co["food_infl","Estimate"]; b2 <- co["food_infl2","Estimate"]; b4 <- co["food_infl:G_c","Estimate"]
fpR <- median(m$food_infl, na.rm=TRUE)
turn_gini <- gbar - (b1 + 2*b2*fpR)/b4

ctry <- m %>% group_by(iso3) %>%
  summarise(gini = mean(gini, na.rm=TRUE), .groups="drop") %>%
  mutate(G_c = gini - gbar, slope = b1 + 2*b2*fpR + b4*G_c, irr = exp(slope))

afr <- ne_countries(scale="medium", continent="Africa", returnclass="sf") %>%
  left_join(ctry, by = c("iso_a3" = "iso3"))

anc  <- afr %>% filter(!is.na(slope))
miss <- is.na(anc$label_x) | is.na(anc$label_y)
if (any(miss)) {
  cc <- st_coordinates(st_point_on_surface(anc[miss, ]))
  anc$label_x[miss] <- cc[, 1]; anc$label_y[miss] <- cc[, 2]
}

Lmax <- 4.5
anc$len  <- rescale(abs(anc$slope), to = c(1.2, Lmax))
anc$xend <- anc$label_x
anc$yend <- anc$label_y + ifelse(anc$slope >= 0, anc$len, -anc$len)

lab <- st_drop_geometry(anc)
lab$lab_text <- sprintf("%s\n%+.3f", lab$name, lab$slope)   # ad + marjinal etki

p <- ggplot() +
  geom_sf(data = afr, aes(fill = gini), color = "white", linewidth = 0.2) +
  scale_fill_distiller("Gini\n(inequality)", palette = "YlOrBr",
                       direction = 1, na.value = "grey92") +
  geom_segment(data = anc,
               aes(x = label_x, y = label_y, xend = xend, yend = yend, color = slope),
               arrow = arrow(length = unit(0.11, "cm"), type = "closed"),
               linewidth = 0.7, lineend = "round") +
  scale_color_gradient2("Marginal effect of\na food price shock",
                        low = "#2166AC", mid = "grey60", high = "#B2182B", midpoint = 0) +
  geom_text_repel(data = lab,
                  aes(x = xend, y = yend, label = lab_text),
                  size = 2.3, lineheight = 0.85,
                  color = "grey15",
                  family = "Times New Roman", fontface = "bold",
                  bg.color = "white", bg.r = 0.18,
                  segment.size = 0.3, segment.color = "grey55",
                  min.segment.length = 0, box.padding = 0.45, point.padding = 0.15,
                  max.overlaps = Inf, seed = 42) +
  coord_sf(xlim = c(-21, 58), ylim = c(-37, 29), expand = FALSE) +
  theme_void(base_size = 12) +
  theme(legend.position = "right",
        legend.box.spacing = unit(1.1, "cm"),   
        legend.margin = margin(0, 0, 0, 8))

fig_dir <- "figures"

ggsave(file.path(fig_dir, "fig3_ssa_inequality_arrows.png"),
       p, width = 12, height = 10, dpi = 300, bg = "white", type = "cairo")

ggsave(file.path(fig_dir, "fig3_ssa_inequality_arrows.pdf"),
       p, width = 12, height = 10, device = cairo_pdf, bg = "white")


options(scipen=999)
library(glmmTMB)
m$logP <- log(m$population); m$G_c <- m$gini - mean(m$gini,na.rm=TRUE)

base_fit <- function(dat) {
  glmmTMB(unrest_n ~ food_infl + food_infl2 + G_c + food_infl:G_c +
            log_gdp + urban_pct + food_import_pct + undernourish_pct +
            inflation_pct + democracy + (1|iso3) + (1|yearf),
          offset=logP, family=nbinom2, data=dat)
}

cat("===== BASELINE (teyit) =====\n")
v <- summary(base_fit(m))$coefficients$cond["food_infl:G_c", c("Estimate","Std. Error","Pr(>|z|)")]
cat(sprintf("Full sample: b=%.5f SE=%.5f p=%.4f\n", v[1],v[2],v[3]))

cat("\n===== 1. SADECE Guney Afrika cikar (EN KRITIK) =====\n")
v <- summary(base_fit(m[m$iso3 != "ZAF",]))$coefficients$cond["food_infl:G_c", c("Estimate","Std. Error","Pr(>|z|)")]
cat(sprintf("ZAF cikarilmis: b=%.5f SE=%.5f p=%.4f  %s\n", v[1],v[2],v[3], ifelse(v[3]<0.05,"SAGLAM","COKTU")))

cat("\n===== 2. En esitsiz 4 guney ulke cikar (ZAF,NAM,SWZ,BWA) =====\n")
v <- summary(base_fit(m[!(m$iso3 %in% c("ZAF","NAM","SWZ","BWA")),]))$coefficients$cond["food_infl:G_c", c("Estimate","Std. Error","Pr(>|z|)")]
cat(sprintf("Guney-4 cikarilmis: b=%.5f SE=%.5f p=%.4f  %s\n", v[1],v[2],v[3], ifelse(v[3]<0.05,"SAGLAM","COKTU")))

cat("\n===== 3. Pozitif-egim tum ulkeler cikar (ZAF,NAM,SWZ,BWA - tipping ustu) =====\n")
cat("(yukaridaki ile ayni grup - pozitif marjinal etkili olanlar)\n")

cat("\n===== 4. LEAVE-ONE-COUNTRY-OUT JACKKNIFE (tum ulkeler) =====\n")
ctrys <- sort(unique(m$iso3[!is.na(m$food_infl)]))
jack <- data.frame(iso3=character(), b=numeric(), se=numeric(), p=numeric(), stringsAsFactors=FALSE)
for (cc in ctrys) {
  M <- tryCatch(base_fit(m[m$iso3 != cc,]), error=function(e) NULL)
  if (!is.null(M)) {
    s <- summary(M)$coefficients$cond
    if ("food_infl:G_c" %in% rownames(s)) {
      vv <- s["food_infl:G_c", c("Estimate","Std. Error","Pr(>|z|)")]
      jack <- rbind(jack, data.frame(iso3=cc, b=vv[1], se=vv[2], p=vv[3]))
    }
  }
}
jack <- jack[order(jack$p, decreasing=TRUE),]
cat("\nEtkiyi en cok tasiyan 10 ulke (cikarilinca p en cok artan):\n")
print(head(jack, 10), row.names=FALSE, digits=4)
cat("\n--- OZET ---\n")
cat("Jackknife b araligi: [", round(min(jack$b),5), ",", round(max(jack$b),5), "]\n")
cat("Jackknife p araligi: [", round(min(jack$p),4), ",", round(max(jack$p),4), "]\n")
cat("p<0.05 kalan (etki saglam):", sum(jack$p<0.05), "/", nrow(jack), "ulke cikariminda\n")
cat("p>0.05 olan (etki kayboluyor):", sum(jack$p>0.05), "/", nrow(jack), "ulke cikariminda\n")
cat("Tum b'ler pozitif mi:", all(jack$b>0), "\n")











options(scipen=999)
library(fixest)
m$u_rate <- log1p((m$unrest_n / m$population) * 1e6)
m$G_c    <- m$gini - mean(m$gini, na.rm=TRUE)

cat("===== mevcut enstrumanlar ve maruziyet degiskenleri =====\n")
cat("global_food_infl var mi:", "global_food_infl" %in% names(m), "| NA:", sum(is.na(m$global_food_infl)),"\n")
cat("food_import_pct var mi:", "food_import_pct" %in% names(m), "| NA:", sum(is.na(m$food_import_pct)),"\n")

imp_share <- tapply(m$food_import_pct, m$iso3, mean, na.rm=TRUE)
m$imp_base <- imp_share[as.character(m$iso3)]
m$imp_base_c <- m$imp_base - mean(m$imp_base, na.rm=TRUE)

m$Z_bartik  <- m$global_food_infl * (m$imp_base / 100)
m$Z_bartikG <- m$Z_bartik * m$G_c
m$fp_Gc     <- m$food_infl * m$G_c

cat("\n===== BARTIK IV: food_infl + interaction, iki enstruman =====\n")
iv_bartik <- tryCatch(
  feols(u_rate ~ log_gdp + urban_pct + food_import_pct + undernourish_pct +
          inflation_pct + democracy | iso3 + yearf |
          food_infl + fp_Gc ~ Z_bartik + Z_bartikG,
        data = m, cluster = ~iso3),
  error=function(e){cat("ERR:",conditionMessage(e),"\n");NULL})
if(!is.null(iv_bartik)){
  print(summary(iv_bartik))
  cat("\n--- first-stage F ---\n")
  print(fitstat(iv_bartik, ~ ivf1))
}

cat("\n===== KARSILASTIRMA: eski enstruman (direkt global) first-stage =====\n")
m$Zg  <- m$global_food_infl
m$ZgG <- m$Zg * m$G_c
iv_old <- tryCatch(
  feols(u_rate ~ log_gdp + urban_pct + food_import_pct + undernourish_pct +
          inflation_pct + democracy | iso3 + yearf |
          food_infl + fp_Gc ~ Zg + ZgG,
        data = m, cluster = ~iso3),
  error=function(e){cat("ERR:",conditionMessage(e),"\n");NULL})
if(!is.null(iv_old)){
  cat("--- eski: first-stage F ---\n")
  print(fitstat(iv_old, ~ ivf1))
  cat("eski interaction (fp_Gc):\n")
  print(round(coef(summary(iv_old))["fit_fp_Gc"],5))
}




options(scipen=999)
library(fixest)
m$logP <- log(m$population); m$G_c <- m$gini - mean(m$gini,na.rm=TRUE)
m$fp_Gc <- m$food_infl * m$G_c

cat("===== fixest gercek FE-NB (fenegbin) - duzeltilmis cikti =====\n")
M2_fe_fix <- fenegbin(unrest_n ~ food_infl + food_infl2 + G_c + fp_Gc +
                        log_gdp + urban_pct + food_import_pct + undernourish_pct +
                        inflation_pct + democracy + offset(logP) | iso3 + yearf,
                      data = m, cluster = ~iso3)
print(summary(M2_fe_fix))
cat("\n--- interaction (fp_Gc) ---\n")
ct <- M2_fe_fix$coeftable
print(round(ct["fp_Gc",],5))

cat("\n===== fixest PPML (FE-Poisson) teyit =====\n")
M2_ppml <- fepois(unrest_n ~ food_infl + food_infl2 + G_c + fp_Gc +
                    log_gdp + urban_pct + food_import_pct + undernourish_pct +
                    inflation_pct + democracy + offset(logP) | iso3 + yearf,
                  data = m, cluster = ~iso3)
ctp <- M2_ppml$coeftable
cat("PPML interaction (fp_Gc):\n")
print(round(ctp["fp_Gc",],5))

cat("\n===== OZET: uc FE yontemi =====\n")
cat("glmmTMB dummy-FE: b=0.00223 p=0.004 (belgedeki asil M2 bu)\n")
cat(sprintf("fixest FE-NB:     b=%.5f p=%.4f\n", ct["fp_Gc","Estimate"], ct["fp_Gc",4]))
cat(sprintf("fixest PPML:      b=%.5f p=%.4f\n", ctp["fp_Gc","Estimate"], ctp["fp_Gc",4]))






options(scipen=999)
library(fixest)
library(sandwich); library(lmtest)
m$logP <- log(m$population); m$G_c <- m$gini - mean(m$gini,na.rm=TRUE)
m$fp_Gc <- m$food_infl * m$G_c

cat("===== 1. fixest FE-NB: farkli SE turleri =====\n")
M <- fenegbin(unrest_n ~ food_infl + food_infl2 + G_c + fp_Gc +
                log_gdp + urban_pct + food_import_pct + undernourish_pct +
                inflation_pct + democracy + offset(logP) | iso3 + yearf, data = m)

cat("--- (a) standart (model-based) SE ---\n")
print(round(summary(M, vcov="iid")$coeftable["fp_Gc",],5))
cat("--- (b) cluster-robust by country ---\n")
print(round(summary(M, cluster=~iso3)$coeftable["fp_Gc",],5))
cat("--- (c) cluster by country AND year (two-way) ---\n")
print(round(summary(M, cluster=~iso3+yearf)$coeftable["fp_Gc",],5))

cat("\n===== 2. PPML ayni SE turleri =====\n")
MP <- fepois(unrest_n ~ food_infl + food_infl2 + G_c + fp_Gc +
               log_gdp + urban_pct + food_import_pct + undernourish_pct +
               inflation_pct + democracy + offset(logP) | iso3 + yearf, data = m)
cat("--- PPML iid ---\n");      print(round(summary(MP, vcov="iid")$coeftable["fp_Gc",],5))
cat("--- PPML cluster iso3 ---\n"); print(round(summary(MP, cluster=~iso3)$coeftable["fp_Gc",],5))

cat("\n===== 3. WILD CLUSTER BOOTSTRAP (40 kume icin dogru arac) =====\n")
# fixest fenegbin wild bootstrap dogrudan zor; glmmTMB dummy-FE uzerinden boottest dene
# basit: cluster-robust ama G ayarli (CR3 benzeri) - fixest ssc ayari
cat("--- fixest FE-NB, small-sample adjusted (ssc) ---\n")
print(round(summary(M, cluster=~iso3, ssc=ssc(adj=TRUE, cluster.adj=TRUE))$coeftable["fp_Gc",],5))

cat("\n===== 4. KAC KUME (cluster sayisi) =====\n")
cat("Ulke (cluster) sayisi:", length(unique(m$iso3[!is.na(m$fp_Gc)])), "\n")
cat("40 kume: cluster-robust SE asagi-sapmali olabilir, wild bootstrap onerilir\n")

cat("\n===== OZET =====\n")
cat("Soru: interaction cluster-robust SE altinda anlamli mi?\n")
cat("NB cluster p, PPML cluster p, ve small-sample ayarli p'ye bak\n")


options(scipen=999)
library(fixest); library(parallel)
m$logP <- log(m$population); m$G_c <- m$gini - mean(m$gini,na.rm=TRUE)
m$fp_Gc <- m$food_infl * m$G_c
mv <- c("unrest_n","food_infl","food_infl2","G_c","fp_Gc","log_gdp","urban_pct",
        "food_import_pct","undernourish_pct","inflation_pct","democracy","logP")
d <- m[complete.cases(m[,mv]),]
ctrys <- unique(d$iso3); nC <- length(ctrys)

f <- unrest_n ~ food_infl + food_infl2 + G_c + fp_Gc + log_gdp + urban_pct +
  food_import_pct + undernourish_pct + inflation_pct + democracy +
  offset(logP) | iso3 + yearf

b_ppml <- fepois(f, data=d)$coeftable["fp_Gc","Estimate"]
cat("PPML nokta:", round(b_ppml,5), "\n")

B <- 999
one_pp <- function(bb) {
  set.seed(2000+bb)
  samp <- sample(ctrys, nC, replace=TRUE)
  dl <- do.call(rbind, lapply(seq_along(samp), function(j){
    tmp <- d[d$iso3==samp[j],]; tmp$iso3 <- paste0(tmp$iso3,"_",j); tmp
  }))
  tryCatch(fepois(f, data=dl, glm.iter=50)$coeftable["fp_Gc","Estimate"],
           error=function(e) NA, warning=function(w) NA)
}

nco <- min(14, detectCores())
cl <- makeCluster(nco)
clusterEvalQ(cl, library(fixest))
clusterExport(cl, c("d","ctrys","nC","f"))
t0 <- Sys.time()
res <- parLapply(cl, 1:B, one_pp)
stopCluster(cl)
cat("Sure:", round(difftime(Sys.time(),t0,units="mins"),2), "dk\n")

boot_pp <- unlist(res); boot_pp <- boot_pp[!is.na(boot_pp)]
cat("Basarili:", length(boot_pp), "/", B, "\n")
if (length(boot_pp) > 100) {
  ci <- quantile(boot_pp, c(.025,.975))
  cat("PPML bootstrap SE:", round(sd(boot_pp),5), "\n")
  cat("95% CI: [", round(ci[1],5), ",", round(ci[2],5), "]\n")
  cat("p:", round(2*min(mean(boot_pp<=0), mean(boot_pp>=0)),4),
      ifelse(ci[1]>0 | ci[2]<0, "-> 0 DISINDA (anlamli)", "-> 0 ICINDE (anlamsiz)"), "\n")
} else {
  cat("PPML de yeterince converge etmedi - bu da kirilganlik isareti\n")
}

citation("rnaturalearth")
citation("sf")



options(scipen=999)
library(glmmTMB)
m$logP <- log(m$population); m$G_c <- m$gini - mean(m$gini,na.rm=TRUE)

base_terms <- "food_infl + food_infl2 + G_c + food_infl:G_c + log_gdp + urban_pct + food_import_pct + undernourish_pct + inflation_pct + democracy"

M_base <- glmmTMB(as.formula(paste("unrest_n ~", base_terms, "+ (1|iso3) + (1|yearf)")),
                  offset=logP, family=nbinom2, data=m)

cat("===== RQ1 TERCILE grouping LR test =====\n")

m$und_tercile <- cut(m$undernourish_pct, breaks=quantile(m$undernourish_pct, c(0,1/3,2/3,1), na.rm=TRUE),
                     include.lowest=TRUE, labels=c("low","mid","high"))
M_terc <- tryCatch(glmmTMB(as.formula(paste("unrest_n ~", base_terms,
                                            "+ food_infl:und_tercile + (1|iso3) + (1|yearf)")),
                           offset=logP, family=nbinom2, data=m),
                   error=function(e){cat("ERR:",conditionMessage(e),"\n");NULL})
if(!is.null(M_terc)){
  # ayni ornekleme indir (NA'lar yuzunden)
  lr_terc <- anova(M_base, M_terc)
  print(lr_terc)
  cat(sprintf("\nTERCILE: Chisq = %.3f, df = %d, p = %.4f\n",
              lr_terc$Chisq[2], lr_terc$Df[2], lr_terc$`Pr(>Chisq)`[2]))
}

cat("\n===== RQ1 NATURAL SPLINE LR test =====\n")
library(splines)
M_spl <- tryCatch(glmmTMB(as.formula(paste("unrest_n ~", base_terms,
                                           "+ food_infl:ns(undernourish_pct, df=3) + (1|iso3) + (1|yearf)")),
                          offset=logP, family=nbinom2, data=m),
                  error=function(e){cat("ERR:",conditionMessage(e),"\n");NULL})
if(!is.null(M_spl)){
  lr_spl <- anova(M_base, M_spl)
  print(lr_spl)
  cat(sprintf("\nSPLINE: Chisq = %.3f, df = %d, p = %.4f\n",
              lr_spl$Chisq[2], lr_spl$Df[2], lr_spl$`Pr(>Chisq)`[2]))
}

cat("\n===== OZET (Table 5 Panel B icin) =====\n")
if(!is.null(M_terc)) cat(sprintf("Tercile grouping: chi2(%d) = %.2f, p = %.3f\n",
                                 lr_terc$Df[2], lr_terc$Chisq[2], lr_terc$`Pr(>Chisq)`[2]))
if(!is.null(M_spl)) cat(sprintf("Natural spline: chi2(%d) = %.2f, p = %.3f\n",
                                lr_spl$Df[2], lr_spl$Chisq[2], lr_spl$`Pr(>Chisq)`[2]))


options(scipen=999)
m$logP <- log(m$population); m$G_c <- m$gini - mean(m$gini,na.rm=TRUE)


base_terms <- "food_infl + food_infl2 + G_c + food_infl:G_c + log_gdp + urban_pct + food_import_pct + undernourish_pct + inflation_pct + democracy"

M_base <- glmmTMB(as.formula(paste("unrest_n ~", base_terms, "+ (1|iso3) + (1|yearf)")),
                  offset=logP, family=nbinom2, data=m)

cat("===== RQ1 TERCILE grouping LR test =====\n")
m$und_tercile <- cut(m$undernourish_pct, breaks=quantile(m$undernourish_pct, c(0,1/3,2/3,1), na.rm=TRUE),
                     include.lowest=TRUE, labels=c("low","mid","high"))
M_terc <- tryCatch(glmmTMB(as.formula(paste("unrest_n ~", base_terms,
                                            "+ food_infl:und_tercile + (1|iso3) + (1|yearf)")),
                           offset=logP, family=nbinom2, data=m),
                   error=function(e){cat("ERR:",conditionMessage(e),"\n");NULL})
if(!is.null(M_terc)){
  lr_terc <- anova(M_base, M_terc)
  print(lr_terc)
  cat(sprintf("\nTERCILE: Chisq = %.3f, df = %d, p = %.4f\n",
              lr_terc$Chisq[2], lr_terc$Df[2], lr_terc$`Pr(>Chisq)`[2]))
}

cat("\n===== RQ1 NATURAL SPLINE LR test =====\n")
M_spl <- tryCatch(glmmTMB(as.formula(paste("unrest_n ~", base_terms,
                                           "+ food_infl:ns(undernourish_pct, df=3) + (1|iso3) + (1|yearf)")),
                          offset=logP, family=nbinom2, data=m),
                  error=function(e){cat("ERR:",conditionMessage(e),"\n");NULL})
if(!is.null(M_spl)){
  lr_spl <- anova(M_base, M_spl)
  print(lr_spl)
  cat(sprintf("\nSPLINE: Chisq = %.3f, df = %d, p = %.4f\n",
              lr_spl$Chisq[2], lr_spl$Df[2], lr_spl$`Pr(>Chisq)`[2]))
}

cat("\n===== OZET (Table 5 Panel B icin) =====\n")
if(!is.null(M_terc)) cat(sprintf("Tercile grouping: chi2(%d) = %.2f, p = %.3f\n",
                                 lr_terc$Df[2], lr_terc$Chisq[2], lr_terc$`Pr(>Chisq)`[2]))
if(!is.null(M_spl)) cat(sprintf("Natural spline: chi2(%d) = %.2f, p = %.3f\n",
                                lr_spl$Df[2], lr_spl$Chisq[2], lr_spl$`Pr(>Chisq)`[2]))

options(scipen=999)
library(glmmTMB); library(splines)
m$logP <- log(m$population)
gini_mean <- mean(m$gini, na.rm=TRUE)   # <-- eksik olan buydu
m$G_c <- m$gini - gini_mean

# RQ1 testleri icin ortak baz (dummy-FE, M2 yapisi)
base_rhs <- "food_infl * G_c + food_infl2 + log_gdp + urban_pct + food_import_pct + undernourish_pct + inflation_pct + democracy + iso3 + yearf + offset(logP)"

# ayni ornekleme sabitlemek icin once tum RQ1 degiskenlerini iceren complete-case
m$und_tercile <- cut(m$undernourish_pct,
                     breaks=quantile(m$undernourish_pct, c(0,1/3,2/3,1), na.rm=TRUE),
                     include.lowest=TRUE, labels=c("low","mid","high"))
rq_vars <- c("unrest_n","food_infl","food_infl2","G_c","log_gdp","urban_pct",
             "food_import_pct","undernourish_pct","inflation_pct","democracy",
             "und_tercile","logP","iso3","yearf")
dd <- m[complete.cases(m[, rq_vars]), ]
cat("RQ1 complete-case N =", nrow(dd), "\n\n")

# baz model (ayni ornekte)
M_base <- glmmTMB(as.formula(paste("unrest_n ~", base_rhs)),
                  family=nbinom2, data=dd)

cat("===== RQ1 TERCILE grouping LR test =====\n")
M_terc <- glmmTMB(as.formula(paste("unrest_n ~", base_rhs, "+ food_infl:und_tercile")),
                  family=nbinom2, data=dd)
lr_terc <- anova(M_base, M_terc)
print(lr_terc)
cat(sprintf("\n>> TERCILE: chi2(%d) = %.2f, p = %.4f\n",
            lr_terc$`Chi Df`[2], lr_terc$Chisq[2], lr_terc$`Pr(>Chisq)`[2]))

cat("\n===== RQ1 NATURAL SPLINE LR test =====\n")
M_spl <- glmmTMB(as.formula(paste("unrest_n ~", base_rhs, "+ food_infl:ns(undernourish_pct, df=3)")),
                 family=nbinom2, data=dd)
lr_spl <- anova(M_base, M_spl)
print(lr_spl)
cat(sprintf("\n>> SPLINE: chi2(%d) = %.2f, p = %.4f\n",
            lr_spl$`Chi Df`[2], lr_spl$Chisq[2], lr_spl$`Pr(>Chisq)`[2]))

cat("\n===== OZET (Table 5 Panel B) =====\n")
cat(sprintf("Tercile grouping (LR): chi2(%d) = %.2f, p = %.3f\n",
            lr_terc$`Chi Df`[2], lr_terc$Chisq[2], lr_terc$`Pr(>Chisq)`[2]))
cat(sprintf("Natural spline (LR):   chi2(%d) = %.2f, p = %.3f\n",
            lr_spl$`Chi Df`[2], lr_spl$Chisq[2], lr_spl$`Pr(>Chisq)`[2]))













options(scipen=999)
library(fixest)
m$logP <- log(m$population)
gini_mean <- mean(m$gini, na.rm=TRUE)
m$G_c <- m$gini - gini_mean
m$fp_Gc <- m$food_infl * m$G_c

M <- fenegbin(unrest_n ~ food_infl + food_infl2 + G_c + fp_Gc +
                log_gdp + urban_pct + food_import_pct + undernourish_pct +
                inflation_pct + democracy + offset(logP) | iso3 + yearf, data = m)

b  <- coef(M)
nm <- names(b)
# vcov'dan .theta'yi cikar, sadece katsayi blogunu al
V_iid <- vcov(M, vcov="iid")[nm, nm]
V_cl  <- vcov(M, cluster=~iso3)[nm, nm]

fp_med <- median(m$food_infl, na.rm=TRUE)
gini_seq <- seq(quantile(m$gini,.05,na.rm=TRUE), quantile(m$gini,.95,na.rm=TRUE), length.out=100)
gc_seq <- gini_seq - gini_mean

out <- data.frame(gini=gini_seq, slope=NA, se_iid=NA, se_cl=NA)
for (i in seq_along(gc_seq)) {
  g <- setNames(rep(0,length(b)), nm)
  g["food_infl"]  <- 1
  g["food_infl2"] <- 2*fp_med
  g["fp_Gc"]      <- gc_seq[i]
  est <- sum(g*b)
  out$slope[i]  <- est
  out$se_iid[i] <- sqrt(as.numeric(t(g) %*% V_iid %*% g))
  out$se_cl[i]  <- sqrt(as.numeric(t(g) %*% V_cl  %*% g))
}
out$lo_iid <- out$slope - 1.96*out$se_iid
out$hi_iid <- out$slope + 1.96*out$se_iid
out$lo_cl  <- out$slope - 1.96*out$se_cl
out$hi_cl  <- out$slope + 1.96*out$se_cl

zc <- approx(out$slope, out$gini, xout=0)$y
cat("Zero-crossing (slope=0) Gini =", round(zc,2), "\n")
cat("P25/P50/P75/P90 Gini:", round(quantile(m$gini,c(.25,.5,.75,.9),na.rm=TRUE),1), "\n")
cat("Model-based bandi sifiri DISLAYAN nokta:", sum(out$lo_iid>0 | out$hi_iid<0), "/100\n")
cat("Cluster bandi sifiri her yerde iceriyor mu:", all(out$lo_cl < 0 & out$hi_cl > 0), "\n")
cat("Cluster bandi sifiri DISLAYAN nokta:", sum(out$lo_cl>0 | out$hi_cl<0), "/100\n\n")

fig_dir <- "figures"
write.csv(out, file.path(fig_dir, "fig2_data.csv"), row.names=FALSE)
cat("fig2_data.csv yazildi\n")
print(head(out, 3))
print(tail(out, 3))







options(scipen=999)
library(fixest)
m$logP <- log(m$population)
gini_mean <- mean(m$gini, na.rm=TRUE)
m$G_c <- m$gini - gini_mean
m$fp_Gc <- m$food_infl * m$G_c

f <- unrest_n ~ food_infl + food_infl2 + G_c + fp_Gc +
  log_gdp + urban_pct + food_import_pct + undernourish_pct +
  inflation_pct + democracy + offset(logP) | iso3 + yearf

# Southern high-Gini 4 ulke disla
south4 <- c("ZAF","NAM","SWZ","BWA")
M_noS4 <- fenegbin(f, data = m[!(m$iso3 %in% south4), ])
ct <- M_noS4$coeftable["fp_Gc", ]
b <- ct["Estimate"]; se <- ct["Std. Error"]; p <- ct["Pr(>|z|)"]
ci_lo <- b - 1.96*se; ci_hi <- b + 1.96*se

cat("===== Excluding Southern high-Gini (ZAF, NAM, SWZ, BWA) =====\n")
cat(sprintf("beta = %.5f\n", b))
cat(sprintf("SE   = %.5f\n", se))
cat(sprintf("95%% CI = [%.4f, %.4f]\n", ci_lo, ci_hi))
cat(sprintf("p    = %.3f\n", p))
cat(sprintf("N    = %d\n", M_noS4$nobs))
cat("\nTable 5 satiri icin:\n")
cat(sprintf("Excluding high-inequality Southern African countries (ZAF, NAM, SWZ, BWA) | %.5f | %.5f | [%.4f, %.4f] | %.3f\n",
            b, se, ci_lo, ci_hi, p))

# karsilastirma: model-based de verelim (Table 5 model-based p sutunu kullaniyor)
M_iid <- summary(M_noS4, vcov="iid")$coeftable["fp_Gc",]
cat(sprintf("\n(model-based p = %.3f; Table 5 model-based p sutunu kullaniyorsa bunu koy)\n", M_iid["Pr(>|z|)"]))


