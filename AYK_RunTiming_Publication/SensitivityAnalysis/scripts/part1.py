import sys; sys.path.insert(0,'.')
import os, numpy as np, pandas as pd, sa_core as s

THRESHOLD = 0.7
OUTDIR = "/sessions/nice-sweet-cori/mnt/AYK_RunTiming/AYK_RunTiming_Publication/SensitivityAnalysis/1_Individual_Assignment_Concentration"
os.makedirs(os.path.join(OUTDIR,"data"), exist_ok=True)
os.makedirs(os.path.join(OUTDIR,"figures"), exist_ok=True)

edges = s.load_edges()
# group membership matrix (G x E)
G = s.GROUP_ORDER
grp = edges["group"].to_numpy()
group_idx = {g:i for i,g in enumerate(G)}
Gmat = np.zeros((len(G), len(edges)))
for e,g in enumerate(grp):
    if isinstance(g,str) and g in group_idx:
        Gmat[group_idx[g], e] = 1.0
in_group_edge = np.array([isinstance(g,str) for g in grp])

per_fish = []
for yr in s.YEARS:
    nat = s.load_natal(yr)
    iso = nat["natal_iso"].to_numpy(dtype=float)
    resc = s.rescaled_matrix(edges, iso, THRESHOLD)     # (E x F)
    total_all = resc.sum(axis=0)                         # (F,)
    grouped   = resc[in_group_edge,:].sum(axis=0)        # mass within the 6 groups
    gsums = Gmat @ resc                                  # (G x F) mass per group
    dom_i = gsums.argmax(axis=0)
    dom_mass = gsums.max(axis=0)
    dom_name = [G[i] for i in dom_i]
    n_groups_hit = (gsums > 0).sum(axis=0)               # # of groups receiving mass
    n_edges_hit = (resc > 0).sum(axis=0)                 # # of stream edges retained
    with np.errstate(divide="ignore", invalid="ignore"):
        prop_dom_of_grouped = np.where(grouped>0, dom_mass/grouped, np.nan)
        prop_dom_of_all     = np.where(total_all>0, dom_mass/total_all, np.nan)
        prop_outside        = np.where(total_all>0, (total_all-grouped)/total_all, np.nan)
    df = pd.DataFrame({
        "year": yr,
        "fish_id": nat.get("Fish_id", pd.Series([np.nan]*len(nat))).values,
        "oto_num": nat.get("OtoNum", pd.Series([np.nan]*len(nat))).values,
        "DOY": nat["DOY"].values,
        "natal_iso": iso,
        "dominant_group": dom_name,
        "prop_in_dominant_group_of_grouped": prop_dom_of_grouped,
        "prop_in_dominant_group_of_all": prop_dom_of_all,
        "prop_outside_regional_groups": prop_outside,
        "n_groups_with_mass": n_groups_hit,
        "n_stream_edges_retained": n_edges_hit,
    })
    per_fish.append(df)

pf = pd.concat(per_fish, ignore_index=True)
pf.to_csv(os.path.join(OUTDIR,"data","per_individual_group_concentration.csv"), index=False)

# ---------------- summary ----------------
def pct(x): return round(100*x,1)
main = pf["prop_in_dominant_group_of_grouped"].dropna()
allm = pf["prop_in_dominant_group_of_all"].dropna()
summary = {
    "n_individuals": len(pf),
    "metric": "proportion of each fish's assignment in its single dominant regional group",
    "PRIMARY_within_grouped_mean": round(main.mean(),4),
    "PRIMARY_within_grouped_median": round(main.median(),4),
    "PRIMARY_within_grouped_q25": round(main.quantile(.25),4),
    "PRIMARY_within_grouped_q75": round(main.quantile(.75),4),
    "PRIMARY_within_grouped_min": round(main.min(),4),
    "pct_fish_dominant>=0.5_grouped": pct((main>=0.5).mean()),
    "pct_fish_dominant>=0.7_grouped": pct((main>=0.7).mean()),
    "pct_fish_dominant>=0.9_grouped": pct((main>=0.9).mean()),
    "pct_fish_dominant==1.0_grouped": pct((main>=0.999).mean()),
    "SECONDARY_of_all_mass_mean": round(allm.mean(),4),
    "SECONDARY_of_all_mass_median": round(allm.median(),4),
    "mean_prop_outside_groups": round(pf["prop_outside_regional_groups"].mean(),4),
    "mean_n_groups_with_mass": round(pf["n_groups_with_mass"].mean(),3),
    "mean_n_stream_edges_retained": round(pf["n_stream_edges_retained"].mean(),1),
}
pd.DataFrame([summary]).T.rename(columns={0:"value"}).to_csv(os.path.join(OUTDIR,"data","summary_statistics.csv"))

# summary by dominant group and by year
pf.groupby("dominant_group").agg(
    n_fish=("prop_in_dominant_group_of_grouped","size"),
    mean_prop_in_dom=("prop_in_dominant_group_of_grouped","mean"),
    median_prop_in_dom=("prop_in_dominant_group_of_grouped","median"),
).round(4).to_csv(os.path.join(OUTDIR,"data","by_dominant_group.csv"))
pf.groupby("year").agg(
    n_fish=("prop_in_dominant_group_of_grouped","size"),
    mean_prop_in_dom=("prop_in_dominant_group_of_grouped","mean"),
    median_prop_in_dom=("prop_in_dominant_group_of_grouped","median"),
    mean_prop_outside=("prop_outside_regional_groups","mean"),
).round(4).to_csv(os.path.join(OUTDIR,"data","by_year.csv"))

import json
print(json.dumps(summary, indent=2))
print("\nby dominant group:")
print(pf.groupby("dominant_group")["prop_in_dominant_group_of_grouped"].agg(["size","mean","median"]).round(3))
