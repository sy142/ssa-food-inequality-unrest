# Data sources

The processed analytic panel (`ssa_food_unrest_panel.xlsx`) is derived from the following public sources. Raw data are not redistributed in this repository; download them from the providers below and run `code/01_build_panel.py` to reconstruct the panel.

| Variable(s) | Description | Source | Access |
|---|---|---|---|
| `unrest_n` | Annual count of protest and riot events, country-year | ACLED (Armed Conflict Location & Event Data) | https://acleddata.com — subject to ACLED terms of use; requires free registration and an API key |
| `gini` | Income inequality (Gini, disposable income) | SWIID (Standardized World Income Inequality Database) | https://fsolt.org/swiid/ |
| `gini_wb` | Income inequality (Gini), World Bank, used as fallback | World Bank World Development Indicators (`SI.POV.GINI`) | https://data.worldbank.org |
| `food_cpi`, `food_infl`, `food_shock` | Food consumer price index, year-over-year food inflation, and detrended food-price shock | FAOSTAT Consumer Prices (domain CP) | https://www.fao.org/faostat/ |
| `global_food_infl` | Leave-one-out global median of food inflation | Derived from FAOSTAT food inflation across all available countries | computed in `01_build_panel.py` |
| `gdp_pc` | GDP per capita, constant 2015 US$ | World Bank WDI (`NY.GDP.PCAP.KD`) | https://data.worldbank.org |
| `urban_pct` | Urban population, % of total | World Bank WDI (`SP.URB.TOTL.IN.ZS`) | https://data.worldbank.org |
| `food_import_pct` | Food imports, % of merchandise imports | World Bank WDI (`TM.VAL.FOOD.ZS.UN`) | https://data.worldbank.org |
| `undernourish_pct` | Prevalence of undernourishment | World Bank WDI (`SN.ITK.DEFC.ZS`) | https://data.worldbank.org |
| `inflation_pct` | Headline CPI inflation, annual % | World Bank WDI (`FP.CPI.TOTL.ZG`) | https://data.worldbank.org |
| `oop_health_pct` | Out-of-pocket health expenditure, % | World Bank WDI (`SH.XPD.OOPC.CH.ZS`) | https://data.worldbank.org |
| `population` | Total population (offset for count models) | World Bank WDI (`SP.POP.TOTL`) | https://data.worldbank.org |
| `gov_effect` | Government effectiveness estimate | World Bank Worldwide Governance Indicators (`GE.EST`) | https://info.worldbank.org/governance/wgi/ |
| `democracy` | Electoral democracy index (v2x_polyarchy) | V-Dem (Varieties of Democracy) v16, Country-Year Core | https://www.v-dem.net |
| `temp_anomaly` | Temperature change anomaly | FAOSTAT Temperature change (domain ET) | https://www.fao.org/faostat/ |

Country names were harmonised to ISO-3166 alpha-3 codes with the `country_converter` package. The panel spans 2000–2025; South Sudan (`SDN`) is excluded from the analytic sample. The panel is unbalanced (coverage varies by variable and country).

## Notes on redistribution

- **ACLED** data may not be redistributed; obtain them directly from ACLED under their terms of use.
- **SWIID** and **V-Dem** are downloaded as files from their providers and are not included here.
- World Bank and FAOSTAT series are retrieved programmatically by the build script (`01_build_panel.py`).
