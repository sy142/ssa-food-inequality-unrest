import os
import time
import logging
import numpy as np
import pandas as pd
import country_converter as coco

logging.getLogger("country_converter").setLevel(logging.CRITICAL)

YEAR_MIN, YEAR_MAX = 2000, 2025
YEARS = list(range(YEAR_MIN, YEAR_MAX + 1))
OUTDIR = "."
os.makedirs(OUTDIR, exist_ok=True)

cc = coco.CountryConverter()


def to_iso3(names):
    out = cc.convert(list(names), to="ISO3", not_found=None)
    if not isinstance(out, list):
        out = [out]
    return [x if (isinstance(x, str) and len(x) == 3 and x.isalpha()) else np.nan for x in out]


WB_INDICATORS = {
    "NY.GDP.PCAP.KD":     "gdp_pc",
    "SP.URB.TOTL.IN.ZS":  "urban_pct",
    "TM.VAL.FOOD.ZS.UN":  "food_import_pct",
    "SN.ITK.DEFC.ZS":     "undernourish_pct",
    "FP.CPI.TOTL.ZG":     "inflation_pct",
    "SP.POP.TOTL":        "population",
    "GE.EST":             "gov_effect",
    "SH.XPD.OOPC.CH.ZS":  "oop_health_pct",
    "SI.POV.GINI":        "gini_wb",
}


def _wb_long(series_codes, ssa, db, use_time=True):
    import wbgapi as wb
    kwargs = dict(economy=ssa, labels=False, skipBlanks=False, db=db)
    if use_time:
        kwargs["time"] = YEARS
    last = None
    for attempt in range(1, 6):
        try:
            raw = wb.data.DataFrame(series_codes, **kwargs).reset_index()
            break
        except Exception as e:
            last = e
            wait = 10 * attempt
            print(f"[WB]      attempt {attempt}/5 failed ({type(e).__name__}); retrying in {wait}s...")
            time.sleep(wait)
    else:
        raise last
    if "series" in raw.columns:
        long = raw.melt(id_vars=["economy", "series"], var_name="year", value_name="value")
    else:
        long = raw.melt(id_vars=["economy"], var_name="year", value_name="value")
        long["series"] = series_codes[0]
    long["year"] = long["year"].astype(str).str.replace("YR", "", regex=False)
    long = long[long["year"].str.fullmatch(r"\d{4}")]
    long["year"] = long["year"].astype(int)
    long = long[(long["year"] >= YEAR_MIN) & (long["year"] <= YEAR_MAX)]
    wide = long.pivot_table(index=["economy", "year"], columns="series", values="value").reset_index()
    return wide.rename(columns={"economy": "iso3"})


def get_gov_effect(ssa):
    import requests
    ssa_set = set(ssa)
    frames, page = [], 1
    try:
        while True:
            js = None
            for attempt in range(1, 6):
                r = requests.get(
                    "https://api.worldbank.org/v2/country/all/indicator/GE.EST",
                    params={"format": "json", "source": 3, "date": f"{YEAR_MIN}:{YEAR_MAX}",
                            "per_page": 20000, "page": page},
                    timeout=120)
                if r.status_code == 200:
                    js = r.json()
                    break
                print(f"[WGI]     attempt {attempt}/5 HTTP {r.status_code}; retrying...")
                time.sleep(10 * attempt)
            if js is None or not isinstance(js, list) or len(js) < 2 or not js[1]:
                break
            for row in js[1]:
                frames.append({"iso3": row.get("countryiso3code"),
                               "year": int(row["date"]),
                               "gov_effect": row.get("value")})
            if page >= js[0].get("pages", 1):
                break
            page += 1
    except Exception as e:
        print(f"[WGI]     GE.EST fetch failed ({e}) -> gov_effect skipped.")
        return None
    if not frames:
        print("[WGI]     GE.EST returned no data -> gov_effect skipped.")
        return None
    g = pd.DataFrame(frames)
    g = g[g["iso3"].isin(ssa_set)]
    g = g[(g.year >= YEAR_MIN) & (g.year <= YEAR_MAX)]
    g = g.groupby(["iso3", "year"], as_index=False)["gov_effect"].mean()
    print(f"[WGI]     {g.iso3.nunique()} countries (GE.EST)")
    return g


def get_worldbank():
    import wbgapi as wb
    ssa = list(wb.region.members("SSF"))
    wdi_codes = [k for k in WB_INDICATORS if k != "GE.EST"]
    df = _wb_long(wdi_codes, ssa, 2)
    ge = get_gov_effect(ssa)
    if ge is not None:
        df = df.merge(ge, on=["iso3", "year"], how="left")
    df = df.rename(columns=WB_INDICATORS)
    print(f"[WB]      {df.iso3.nunique()} countries x {df.year.nunique()} years")
    return df


FAO_CP_BULK_URL = "https://bulks-faostat.fao.org/production/ConsumerPriceIndices_E_All_Data_(Normalized).zip"
FAO_CATALOG_URL = "https://bulks-faostat.fao.org/production/datasets_E.xml"
MONTHS_12 = ["January", "February", "March", "April", "May", "June",
             "July", "August", "September", "October", "November", "December"]


def get_faostat_foodcpi(probe=True):
    import requests, zipfile, io
    r = requests.get(FAO_CP_BULK_URL, timeout=180)
    r.raise_for_status()
    zf = zipfile.ZipFile(io.BytesIO(r.content))
    csvs = [n for n in zf.namelist() if n.lower().endswith(".csv")
            and "flag" not in n.lower() and "symbol" not in n.lower()]
    name = csvs[0] if csvs else zf.namelist()[0]
    with zf.open(name) as f:
        data = f.read()
    try:
        raw = pd.read_csv(io.BytesIO(data), encoding="utf-8", low_memory=False)
    except UnicodeDecodeError:
        raw = pd.read_csv(io.BytesIO(data), encoding="latin-1", low_memory=False)

    if probe:
        print("[FAOSTAT] Items:", sorted(raw["Item"].dropna().unique()))
        print("[FAOSTAT] Months:", sorted(raw["Months"].dropna().unique()))

    food = raw[raw["Item"] == "Consumer Prices, Food Indices (2015 = 100)"].copy()
    food = food[food["Months"].isin(MONTHS_12)].copy()
    food = food.rename(columns={"Area": "country", "Year": "year", "Value": "food_cpi"})
    food["food_cpi"] = pd.to_numeric(food["food_cpi"], errors="coerce")
    food = food.dropna(subset=["food_cpi"])
    annual = food.groupby(["country", "year"], as_index=False)["food_cpi"].mean()
    annual["iso3"] = to_iso3(annual["country"])
    annual = annual.dropna(subset=["iso3"])
    annual["year"] = annual["year"].astype(int)
    annual = annual[(annual.year >= YEAR_MIN) & (annual.year <= YEAR_MAX)]
    annual = annual.groupby(["iso3", "year"], as_index=False)["food_cpi"].mean()
    print(f"[FAOSTAT] {annual.iso3.nunique()} countries x {annual.year.nunique()} years")
    return annual[["iso3", "year", "food_cpi"]]


def derive_food_price_shock(df):
    df = df.sort_values(["iso3", "year"]).reset_index(drop=True).copy()
    df["food_infl"] = df.groupby("iso3")["food_cpi"].pct_change() * 100
    df["global_food_infl"] = df.groupby("year")["food_infl"].transform("median")
    df["_log_fc"] = np.log(df["food_cpi"])

    def _resid(g):
        gg = g.dropna(subset=["_log_fc"])
        if len(gg) < 5:
            return pd.Series(np.nan, index=g.index)
        b1, b0 = np.polyfit(gg["year"].values.astype(float), gg["_log_fc"].values, 1)
        return g["_log_fc"] - (b0 + b1 * g["year"])

    df["food_shock"] = df.groupby("iso3", group_keys=False)[["year", "_log_fc"]].apply(_resid)
    return df.drop(columns=["_log_fc"])


def get_faostat_temperature(probe=True):
    import requests, zipfile, io
    import xml.etree.ElementTree as ET
    cat = requests.get(FAO_CATALOG_URL, timeout=120)
    cat.raise_for_status()
    url = None
    for ds in ET.fromstring(cat.content).findall("Dataset"):
        if "temperature" in (ds.findtext("DatasetName") or "").lower():
            url = ds.findtext("FileLocation")
            break
    if not url:
        print("[FAO-Temp] temperature dataset not found in catalog -> skipping.")
        return None
    r = requests.get(url, timeout=180)
    r.raise_for_status()
    zf = zipfile.ZipFile(io.BytesIO(r.content))
    csvs = [n for n in zf.namelist() if n.lower().endswith(".csv")
            and "flag" not in n.lower() and "symbol" not in n.lower()]
    name = csvs[0] if csvs else zf.namelist()[0]
    with zf.open(name) as f:
        data = f.read()
    try:
        raw = pd.read_csv(io.BytesIO(data), encoding="utf-8", low_memory=False)
    except UnicodeDecodeError:
        raw = pd.read_csv(io.BytesIO(data), encoding="latin-1", low_memory=False)

    if probe:
        print("[FAO-Temp] Elements:", sorted(raw["Element"].dropna().unique()))
        print("[FAO-Temp] Months:", sorted(raw["Months"].dropna().astype(str).unique())[:20])

    t = raw[raw["Element"].astype(str).str.contains("Temperature change", case=False, na=False)].copy()
    annual = t[t["Months"].astype(str).str.contains("Meteorological year", case=False, na=False)]
    if annual.empty:
        annual = t[t["Months"].isin(MONTHS_12)]
    annual = annual.rename(columns={"Area": "country", "Year": "year", "Value": "temp_anomaly"})
    annual["temp_anomaly"] = pd.to_numeric(annual["temp_anomaly"], errors="coerce")
    annual = annual.dropna(subset=["temp_anomaly"])
    annual["iso3"] = to_iso3(annual["country"])
    annual = annual.dropna(subset=["iso3"])
    annual["year"] = annual["year"].astype(int)
    annual = annual[(annual.year >= YEAR_MIN) & (annual.year <= YEAR_MAX)]
    annual = annual.groupby(["iso3", "year"], as_index=False)["temp_anomaly"].mean()
    print(f"[FAO-Temp] {annual.iso3.nunique()} countries x {annual.year.nunique()} years")
    return annual[["iso3", "year", "temp_anomaly"]]


SWIID_PATH = os.path.join(OUTDIR, "swiid_summary.csv")


def get_swiid(path=SWIID_PATH):
    if not os.path.exists(path):
        print(f"[SWIID]   {path} not found -> skipping (WB Gini used as fallback).")
        return None
    ext = os.path.splitext(path)[1].lower()
    if ext == ".csv":
        s = pd.read_csv(path)
    elif ext == ".dta":
        s = pd.read_stata(path)
    elif ext in (".rda", ".rdata"):
        import pyreadr
        objs = pyreadr.read_r(path)
        key = "swiid_summary" if "swiid_summary" in objs else list(objs.keys())[-1]
        s = objs[key]
    else:
        print(f"[SWIID]   unsupported file type {ext} -> skipping.")
        return None
    ccol = "country" if "country" in s.columns else ("country_name" if "country_name" in s.columns else None)
    gcol = "gini_disp" if "gini_disp" in s.columns else next((c for c in s.columns if "gini" in c.lower()), None)
    if ccol is None or gcol is None:
        print("[SWIID]   country/gini column not found -> skipping.")
        return None
    s = s.rename(columns={ccol: "country", gcol: "gini"})
    s["iso3"] = to_iso3(s["country"])
    s = s.dropna(subset=["iso3"])
    s = s[(s.year >= YEAR_MIN) & (s.year <= YEAR_MAX)]
    s = s.groupby(["iso3", "year"], as_index=False)["gini"].mean()
    print(f"[SWIID]   {s.iso3.nunique()} countries")
    return s[["iso3", "year", "gini"]]


VDEM_CSV = os.path.join(OUTDIR, "V-Dem-CY-Core-v16.csv")


def get_vdem(path=VDEM_CSV):
    if not os.path.exists(path):
        print(f"[V-Dem]   {path} not found -> skipping.")
        return None
    cols = ["country_name", "year", "v2x_polyarchy"]
    v = pd.read_csv(path, usecols=lambda c: c in cols)
    v["iso3"] = to_iso3(v["country_name"])
    v = v.dropna(subset=["iso3"])
    v = v[(v.year >= YEAR_MIN) & (v.year <= YEAR_MAX)]
    v = v.rename(columns={"v2x_polyarchy": "democracy"})
    v = v.groupby(["iso3", "year"], as_index=False)["democracy"].mean()
    print(f"[V-Dem]   {v.iso3.nunique()} countries")
    return v[["iso3", "year", "democracy"]]


ACLED_EMAIL = "sedatarslan@uludag.edu.tr"
ACLED_PASSWORD = "Salim.Sedat12"
ACLED_TOKEN_URL = "https://acleddata.com/oauth/token"
ACLED_READ_URL = "https://acleddata.com/api/acled/read"
ACLED_SSA_REGIONS = "1|2|3|4|5"


def get_acled():
    import requests
    email = ACLED_EMAIL or os.environ.get("ACLED_EMAIL") or os.environ.get("ACLED_USERNAME")
    password = ACLED_PASSWORD or os.environ.get("ACLED_PASSWORD")
    if not (email and password):
        print("[ACLED]   ACLED_EMAIL / ACLED_PASSWORD not set -> skipping.")
        return None

    tok = requests.post(ACLED_TOKEN_URL, timeout=60,
                        data={"username": email, "password": password,
                              "grant_type": "password", "client_id": "acled",
                              "scope": "authenticated"})
    tok.raise_for_status()
    headers = {"Authorization": f"Bearer {tok.json()['access_token']}"}

    frames, page = [], 1
    while page <= 1000:
        r = requests.get(ACLED_READ_URL, headers=headers, timeout=120, params={
            "_format": "json",
            "region": ACLED_SSA_REGIONS,
            "event_type": "Protests|Riots",
            "year": 1999, "year_where": ">",
            "fields": "country|year",
            "limit": 5000, "page": page,
        })
        r.raise_for_status()
        data = r.json().get("data", [])
        if not data:
            break
        frames.append(pd.DataFrame(data))
        page += 1
        time.sleep(0.5)
    if not frames:
        return None
    a = pd.concat(frames, ignore_index=True)
    a["iso3"] = to_iso3(a["country"])
    a = a.dropna(subset=["iso3"])
    a["year"] = a["year"].astype(int)
    unrest = a.groupby(["iso3", "year"]).size().reset_index(name="unrest_n")
    print(f"[ACLED]   {unrest.iso3.nunique()} countries, {int(unrest.unrest_n.sum())} events")
    return unrest


def build():
    wb    = get_worldbank()
    fao   = derive_food_price_shock(get_faostat_foodcpi(probe=True))
    temp  = get_faostat_temperature(probe=True)
    swid  = get_swiid()
    vdem  = get_vdem()
    acled = get_acled()

    iso_universe = sorted(wb.iso3.unique())
    grid = pd.MultiIndex.from_product([iso_universe, YEARS], names=["iso3", "year"]).to_frame(index=False)
    panel = grid.merge(wb, on=["iso3", "year"], how="left")
    for extra in (fao, temp, swid, vdem, acled):
        if extra is not None:
            panel = panel.merge(extra, on=["iso3", "year"], how="left")

    if "unrest_n" in panel.columns:
        panel["unrest_n"] = panel["unrest_n"].fillna(0)

    if "gini" in panel.columns:
        panel["gini"] = panel["gini"].fillna(panel.get("gini_wb"))
    else:
        panel["gini"] = panel.get("gini_wb")

    panel = panel.sort_values(["iso3", "year"])
    panel["gini"] = panel.groupby("iso3")["gini"].transform(
        lambda g: g.interpolate(method="linear", limit_direction="both").ffill().bfill())

    panel["Country"] = cc.convert(panel.iso3.tolist(), src="ISO3", to="name_short")
    panel = panel.rename(columns={"year": "Year"})
    lead = ["Country", "Year"]
    rest = [c for c in panel.columns if c not in lead + ["iso3"]]
    panel = panel[lead + rest + ["iso3"]]
    panel = panel.sort_values(["Country", "Year"]).reset_index(drop=True)

    out = os.path.join(OUTDIR, "ssa_food_unrest_panel.xlsx")
    try:
        panel.to_excel(out, index=False)
    except PermissionError:
        out = os.path.join(OUTDIR, f"ssa_food_unrest_panel_{time.strftime('%Y%m%d_%H%M%S')}.xlsx")
        panel.to_excel(out, index=False)
        print("Target file was locked (open in Excel?) -> wrote a new timestamped file.")
    out2 = os.path.join(OUTDIR, "veriler.xlsx")
    try:
        panel.to_excel(out2, index=False)
    except PermissionError:
        out2 = os.path.join(OUTDIR, f"veriler_{time.strftime('%Y%m%d_%H%M%S')}.xlsx")
        panel.to_excel(out2, index=False)
        print("veriler.xlsx was locked (open in Excel?) -> wrote a new timestamped file.")
    print(f"\nSaved -> {out}   shape={panel.shape}")
    print(f"Saved -> {out2}")
    print("Columns:", list(panel.columns))
    return panel


if __name__ == "__main__":
    build()


import os
import numpy as np
import pandas as pd
import country_converter as coco
import logging

logging.getLogger("country_converter").setLevel(logging.CRITICAL)
cc = coco.CountryConverter()

OUTDIR = "."
YEAR_MIN, YEAR_MAX = 2000, 2025

wgi_fp = os.path.join(OUTDIR, "WB_WGI_WIDEF.csv")
w = pd.read_csv(wgi_fp, low_memory=False)

ge = w[(w["INDICATOR"] == "GOV_WGI_GE") & (w["COMP_BREAKDOWN_1"] == "WGI_EST")].copy()

year_cols = [c for c in w.columns if str(c).isdigit()]
ge = ge.melt(id_vars=["REF_AREA"], value_vars=year_cols,
             var_name="Year", value_name="gov_effect")
ge["Year"] = ge["Year"].astype(int)
ge["gov_effect"] = pd.to_numeric(ge["gov_effect"], errors="coerce")

out = cc.convert(ge["REF_AREA"].tolist(), src="ISO3", to="ISO3", not_found=None)
ge["iso3"] = [x if (isinstance(x, str) and len(x) == 3 and x.isalpha()) else np.nan for x in out]
ge = ge.dropna(subset=["iso3"])
ge = ge[(ge.Year >= YEAR_MIN) & (ge.Year <= YEAR_MAX)]
ge = ge.groupby(["iso3", "Year"], as_index=False)["gov_effect"].mean()
print(f"[WGI-file] {ge.iso3.nunique()} countries x {ge.Year.nunique()} years, "
      f"range [{ge.gov_effect.min():.2f}, {ge.gov_effect.max():.2f}]")

panel = pd.read_excel(os.path.join(OUTDIR, "veriler.xlsx"))
if "gov_effect" in panel.columns:
    panel = panel.drop(columns=["gov_effect"])
panel = panel.merge(ge, on=["iso3", "Year"], how="left")

cols = [c for c in panel.columns if c != "iso3"] + ["iso3"]
panel = panel[cols]

for fn in ("veriler.xlsx", "ssa_food_unrest_panel.xlsx"):
    fp = os.path.join(OUTDIR, fn)
    try:
        panel.to_excel(fp, index=False)
    except PermissionError:
        import time
        fp = os.path.join(OUTDIR, fn.replace(".xlsx", f"_{time.strftime('%Y%m%d_%H%M%S')}.xlsx"))
        panel.to_excel(fp, index=False)
    print("Saved ->", fp)

print("shape:", panel.shape)
print("gov_effect coverage:", f"{panel.gov_effect.notna().mean()*100:.1f}%")
print("Columns:", list(panel.columns))















import os, time
import numpy as np
import pandas as pd
import country_converter as coco
import logging

logging.getLogger("country_converter").setLevel(logging.CRITICAL)
cc = coco.CountryConverter()

OUTDIR = "."
YEAR_MIN, YEAR_MAX = 2000, 2025
FAO_CP_BULK_URL = "https://bulks-faostat.fao.org/production/ConsumerPriceIndices_E_All_Data_(Normalized).zip"
MONTHS_12 = ["January","February","March","April","May","June",
             "July","August","September","October","November","December"]

def to_iso3(names):
    out = cc.convert(list(names), to="ISO3", not_found=None)
    if not isinstance(out, list):
        out = [out]
    return [x if (isinstance(x, str) and len(x) == 3 and x.isalpha()) else np.nan for x in out]

import requests, zipfile, io
r = requests.get(FAO_CP_BULK_URL, timeout=180); r.raise_for_status()
zf = zipfile.ZipFile(io.BytesIO(r.content))
name = [n for n in zf.namelist() if n.lower().endswith(".csv")
        and "flag" not in n.lower() and "symbol" not in n.lower()][0]
data = zf.open(name).read()
try:
    raw = pd.read_csv(io.BytesIO(data), encoding="utf-8", low_memory=False)
except UnicodeDecodeError:
    raw = pd.read_csv(io.BytesIO(data), encoding="latin-1", low_memory=False)

food = raw[raw["Item"] == "Consumer Prices, Food Indices (2015 = 100)"].copy()
food = food[food["Months"].isin(MONTHS_12)]
food = food.rename(columns={"Area":"country","Year":"year","Value":"food_cpi"})
food["food_cpi"] = pd.to_numeric(food["food_cpi"], errors="coerce")
food = food.dropna(subset=["food_cpi"])
ann = food.groupby(["country","year"], as_index=False)["food_cpi"].mean()
ann["iso3"] = to_iso3(ann["country"])
ann = ann.dropna(subset=["iso3"])
ann["year"] = ann["year"].astype(int)
ann = ann.groupby(["iso3","year"], as_index=False)["food_cpi"].mean()
ann = ann.sort_values(["iso3","year"])
ann["food_infl_all"] = ann.groupby("iso3")["food_cpi"].pct_change() * 100
print(f"[FAOSTAT] {ann.iso3.nunique()} countries for global LOO median")

def loo_median(s):
    v = s.to_numpy(dtype=float)
    out = np.full(len(v), np.nan)
    for i in range(len(v)):
        rest = np.delete(v, i); rest = rest[np.isfinite(rest)]
        if rest.size: out[i] = np.median(rest)
    return pd.Series(out, index=s.index)

ann["global_food_infl_loo"] = ann.groupby("year")["food_infl_all"].transform(loo_median)
glob = ann[["iso3","year","global_food_infl_loo"]].rename(columns={"year":"Year"})

fp = os.path.join(OUTDIR, "veriler.xlsx")
df = pd.read_excel(fp)
before = df["global_food_infl"].describe()
df = df.merge(glob, on=["iso3","Year"], how="left")
df["global_food_infl"] = df["global_food_infl_loo"]
df = df.drop(columns=["global_food_infl_loo"])
after = df["global_food_infl"].describe()

print("\nBEFORE:\n", before, "\n\nAFTER (204-country leave-one-out):\n", after)
print("coverage:", f"{df['global_food_infl'].notna().mean()*100:.1f}%")

for fn in ("veriler.xlsx", "ssa_food_unrest_panel.xlsx"):
    p = os.path.join(OUTDIR, fn)
    try:
        df.to_excel(p, index=False)
    except PermissionError:
        p = os.path.join(OUTDIR, fn.replace(".xlsx", f"_{time.strftime('%Y%m%d_%H%M%S')}.xlsx"))
        df.to_excel(p, index=False)
    print("Saved ->", p)
print("shape:", df.shape)







