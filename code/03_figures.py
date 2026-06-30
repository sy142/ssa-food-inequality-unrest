import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
import os

DATA_PATH = "data/ssa_food_unrest_panel.xlsx"
FIG_DIR   = "figures"
os.makedirs(FIG_DIR, exist_ok=True)

df = pd.read_excel(DATA_PATH)
df = df[df["iso3"] != "SDN"].copy()

region_map = {
    "BEN":"West","BFA":"West","CIV":"West","CPV":"West","GHA":"West","GIN":"West",
    "GMB":"West","GNB":"West","LBR":"West","MLI":"West","MRT":"West","NER":"West",
    "NGA":"West","SEN":"West","SLE":"West","TGO":"West",
    "BDI":"East","COM":"East","ERI":"East","ETH":"East","KEN":"East","MDG":"East",
    "MUS":"East","MWI":"East","MOZ":"East","RWA":"East","SYC":"East","SSD":"East",
    "TZA":"East","UGA":"East","ZMB":"East","ZWE":"East","SOM":"East",
    "AGO":"Central","CMR":"Central","CAF":"Central","TCD":"Central","COG":"Central",
    "COD":"Central","GNQ":"Central","GAB":"Central","STP":"Central",
    "BWA":"Southern","SWZ":"Southern","LSO":"Southern","NAM":"Southern","ZAF":"Southern",
}
df["region"] = df["iso3"].map(region_map)

mat = df.pivot_table(index="Country", columns="Year", values="unrest_n", aggfunc="sum")
order = mat.sum(axis=1).sort_values(ascending=False).index
mat = mat.loc[order]
years = mat.columns.values
countries = mat.index.values
n_c, n_y = mat.shape

reg_trend = df.groupby(["region","Year"])["unrest_n"].mean().reset_index()

yr = df.groupby("Year").agg(unrest=("unrest_n","mean"),
                            infl=("food_infl","median"),
                            gini=("gini","mean")).reset_index()
def z(x): return (x - np.nanmean(x)) / np.nanstd(x)
yr["z_unrest"] = z(yr["unrest"]); yr["z_infl"] = z(yr["infl"]); yr["z_gini"] = z(yr["gini"])

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": 0.6, "axes.edgecolor": "#444444",
})

fig = plt.figure(figsize=(17, 13))
gs = GridSpec(2, 2, width_ratios=[1, 0.62], height_ratios=[1, 1],
              wspace=0.16, hspace=0.28, left=0.12, right=0.95, top=0.94, bottom=0.08)
ax_main = fig.add_subplot(gs[:, 0])
ax_B    = fig.add_subplot(gs[0, 1])
ax_C    = fig.add_subplot(gs[1, 1])

cmap = plt.get_cmap("magma_r").copy()
cmap.set_bad("#f2f0f0")
vmax = np.nanpercentile(mat.values, 99)
norm = mcolors.LogNorm(vmin=1, vmax=vmax)
data_plot = np.where(np.isnan(mat.values), 0, mat.values)
data_masked = np.ma.masked_where(data_plot < 1, data_plot)

im = ax_main.imshow(data_masked, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")
ax_main.set_xticks(np.arange(n_y)); ax_main.set_xticklabels(years, rotation=90, fontsize=7)
ax_main.set_yticks(np.arange(n_c)); ax_main.set_yticklabels(countries, fontsize=7.5)
ax_main.set_xticks(np.arange(-.5, n_y, 1), minor=True)
ax_main.set_yticks(np.arange(-.5, n_c, 1), minor=True)
ax_main.grid(which="minor", color="white", linewidth=0.5)
ax_main.tick_params(which="minor", length=0); ax_main.tick_params(which="major", length=2)
for s in ["top","right"]: ax_main.spines[s].set_visible(False)
ax_main.set_xlabel("Year", fontsize=10)
ax_main.set_title("(A)  Social unrest by country and year",
                  fontsize=12, weight="bold", loc="left", pad=10)

cbar = fig.colorbar(im, ax=ax_main, orientation="vertical", fraction=0.02, pad=0.015, extend="max")
cbar.set_label("Events (count, log scale)", fontsize=8)
cbar.ax.tick_params(labelsize=7)

region_colors = {"West":"#1b9e77","East":"#d95f02","Central":"#7570b3","Southern":"#e7298a"}
for reg in ["West","East","Central","Southern"]:
    sub = reg_trend[reg_trend["region"]==reg].sort_values("Year")
    ax_B.plot(sub["Year"], sub["unrest_n"], marker="o", markersize=2.5,
              linewidth=1.6, color=region_colors[reg], label=reg)
ax_B.set_xlim(years.min()-0.5, years.max()+0.5)
ax_B.set_xticks(years[::2]); ax_B.set_xticklabels(years[::2], rotation=90, fontsize=7)
ax_B.set_ylabel("Mean events per country", fontsize=9)
ax_B.legend(fontsize=7.5, frameon=False, ncol=2, loc="upper left")
ax_B.grid(axis="y", color="#dddddd", linewidth=0.5)
for s in ["top","right"]: ax_B.spines[s].set_visible(False)
ax_B.set_title("(B)  Regional trends (UN subregions)",
               fontsize=12, weight="bold", loc="left", pad=8)

var_colors = {"z_unrest":"#7a1f3d","z_infl":"#e6a000","z_gini":"#2c6e9c"}
ax_C.axhline(0, color="#bbbbbb", linewidth=0.7)
ax_C.plot(yr["Year"], yr["z_unrest"], marker="o", ms=2.5, lw=1.8,
          color=var_colors["z_unrest"], label="Social unrest")
ax_C.plot(yr["Year"], yr["z_infl"], marker="s", ms=2.5, lw=1.8,
          color=var_colors["z_infl"], label="Food inflation")
ax_C.plot(yr["Year"], yr["z_gini"], marker="^", ms=2.5, lw=1.8,
          color=var_colors["z_gini"], label="Inequality (Gini)")
ax_C.set_xlim(years.min()-0.5, years.max()+0.5)
ax_C.set_xticks(years[::2]); ax_C.set_xticklabels(years[::2], rotation=90, fontsize=7)
ax_C.set_ylabel("Standardised (z-score)", fontsize=9)
ax_C.set_xlabel("Year", fontsize=10)
ax_C.legend(fontsize=7.5, frameon=False, ncol=1, loc="upper left")
ax_C.grid(axis="y", color="#dddddd", linewidth=0.5)
for s in ["top","right"]: ax_C.spines[s].set_visible(False)
ax_C.set_title("(C)  Key variables over time (panel means, standardised)",
               fontsize=12, weight="bold", loc="left", pad=8)

out_png = os.path.join(FIG_DIR, "fig1_unrest_overview.png")
out_pdf = os.path.join(FIG_DIR, "fig1_unrest_overview.pdf")
fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
print("saved:", out_png)
plt.close(fig)

# figure 2

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
import matplotlib as mpl

fig_dir = "figures"
df = pd.read_csv(fig_dir + "/fig2_data.csv")

try:
    mpl.rcParams['font.family'] = 'Times New Roman'
except Exception:
    pass
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['axes.linewidth'] = 0.8

gini = df['gini'].values
slope = df['slope'].values
lo_iid, hi_iid = df['lo_iid'].values, df['hi_iid'].values
lo_cl,  hi_cl  = df['lo_cl'].values,  df['hi_cl'].values

zc = 57.53 
p25, p50, p75, p90 = 40.0, 43.8, 48.2, 55.0

fig, ax = plt.subplots(figsize=(7.5, 5.2))

ax.fill_between(gini, lo_cl, hi_cl, color="#9aa0a6", alpha=0.30,
                label="95% CI (country-clustered)", linewidth=0)
ax.fill_between(gini, lo_iid, hi_iid, color="#1f6f8b", alpha=0.40,
                label="95% CI (model-based)", linewidth=0)
ax.plot(gini, slope, color="#0b3d52", linewidth=2.0, label="Marginal effect", zorder=5)

ax.axhline(0, color="black", linewidth=0.8, linestyle="-", zorder=4)
ax.axvline(zc, color="#7a7a7a", linewidth=0.9, linestyle="--", zorder=3)
ax.text(zc+0.3, ax.get_ylim()[1]*0.92, f"Sign change\n(Gini = {zc:.1f})",
        fontsize=8.5, color="#4a4a4a", va="top", ha="left")

ymin = min(lo_cl)*1.05
for p, lab in [(p25,"P25"),(p50,"P50"),(p75,"P75"),(p90,"P90")]:
    ax.plot([p,p],[ymin, ymin+0.004], color="black", linewidth=0.8, clip_on=False)
    ax.text(p, ymin-0.004, lab, fontsize=7.5, ha="center", va="top", color="#333333")

ax.set_xlabel("Income inequality (Gini index)", fontsize=11)
ax.set_ylabel("Marginal effect of a food price shock\non log unrest", fontsize=11)
ax.set_xlim(gini.min(), gini.max())
ax.tick_params(labelsize=9)

ax.text(0.025, 0.05,
        "Two-way FE negative binomial\nN = 770; 40 countries, 23 years",
        transform=ax.transAxes, fontsize=8, va="bottom", ha="left",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                  edgecolor="#cccccc", linewidth=0.6))

ax.legend(loc="upper left", fontsize=8.5, frameon=True, framealpha=0.9,
          edgecolor="#cccccc")

for s in ["top","right"]:
    ax.spines[s].set_visible(False)

plt.tight_layout()
plt.savefig(fig_dir + "/fig2_marginal_effect.png", dpi=300, bbox_inches="tight")
plt.savefig(fig_dir + "/fig2_marginal_effect.pdf", bbox_inches="tight")
print("Figure 2 kaydedildi (png + pdf)")
print(f"Model-based bant (dar) ve cluster-robust bant (genis) cizildi")
print(f"Zero-crossing: {zc}, percentile P25-P90 isaretlendi")

