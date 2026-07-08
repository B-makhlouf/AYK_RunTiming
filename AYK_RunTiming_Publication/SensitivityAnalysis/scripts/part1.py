"""Part 1 (LENGTH-BASED): for each individual fish, of the total stream length (km)
its assignment covers after the 0.7 threshold, what fraction lies within a single
(dominant) regional group. Retained edge = rescaled value >= threshold (>0)."""
import sys; sys.path.insert(0,'.')
import os, numpy as np, pandas as pd, sa_core as s

THRESHOLD = 0.7
OUTDIR = "/sessions/nice-sweet-cori/mnt/AYK_RunTiming/AYK_RunTiming_Publication/SensitivityAnalysis/1_Individual_Assignment_Concentration"
os.makedirs(os.path.join(OUTDIR,"data"), exist_ok=True)
os.makedirs(os.path.join(OUTDIR,"figures"), exist_ok=True)

edges = s.load_edges()
G = s.GROUP_ORDER
group_idx = {g:i for i,g in enumerate(G)}
km = edges["length_km"].to_numpy()                       # (E,)
grp = edges["group"].to_numpy()
# group membership matrix weighted by length: Gmat[g,e] = km_e if edge in group g
Gmat = np.zeros((len(G), len(edges)))
for e,g in enumerate(grp):
    if isinstance(g,str) and g in group_idx:
        Gmat[group_idx[g], e] = km[e]
in_group_edge = np.array([isinstance(g,str) for g in grp])
km_in_group = km * in_group_edge

per_fish = []
for yr in s.YEARS:
    nat = s.load_natal(yr)
    iso = nat["natal_iso"].to_numpy(dtype=float)
    resc = s.rescaled_matrix(edges, iso, THRESHOLD)      # (E x F)
    retained = (resc > 0).astype(float)                  # binary edge membership
    total_km   = retained.T @ km                         # (F,) total assigned km (all edges)
    grouped_km = retained.T @ km_in_group                # km within the 6 groups
    gkm = Gmat @ retained                                # (G x F) km per group
    dom_i = gkm.argmax(axis=0)
    dom_km = gkm.max(axis=0)
    dom_name = [G[i] for i in dom_i]
    n_groups_hit = (gkm > 0).sum(axis=0)
    n_edges_hit = retained.sum(axis=0).astype(int)
    with np.errstate(divide="ignore", invalid="ignore"):
        prop_dom_of_all     = np.where(total_km>0, dom_km/total_km, np.nan)
        prop_dom_of_grouped = np.where(grouped_km>0, dom_km/grouped_km, np.nan)
        prop_outside        = np.where(total_km>0, (total_km-grouped_km)/total_km, np.nan)
    per_fish.append(pd.DataFrame({
        "year": yr,
        "fish_id": nat.get("Fish_id", pd.Series([np.nan]*len(nat))).values,
        "oto_num": nat.get("OtoNum", pd.Series([np.nan]*len(nat))).values,
        "DOY": nat["DOY"].values,
        "natal_iso": iso,
        "total_assigned_km": total_km,
        "dominant_group": dom_name,
        "dominant_group_km": dom_km,
        "prop_km_in_dominant_group_of_all": prop_dom_of_all,
        "prop_km_in_dominant_group_of_grouped": prop_dom_of_grouped,
        "prop_km_outside_regional_groups": prop_outside,
        "n_groups_covered": n_groups_hit,
        "n_stream_edges_retained": n_edges_hit,
    }))

pf = pd.concat(per_fish, ignore_index=True)
pf.to_csv(os.path.join(OUTDIR,"data","per_individual_group_concentration.csv"), index=False)

def pct(x): return round(100*x,1)
main = pf["prop_km_in_dominant_group_of_all"].dropna()
summary = {
    "n_individuals": len(pf),
    "metric": "proportion of each fish's assigned stream length (km) within its single dominant regional group, threshold 0.7",
    "mean": round(main.mean(),4),
    "median": round(main.median(),4),
    "q25": round(main.quantile(.25),4),
    "q75": round(main.quantile(.75),4),
    "min": round(main.min(),4),
    "max": round(main.max(),4),
    "pct_fish>=0.5": pct((main>=0.5).mean()),
    "pct_fish>=0.7": pct((main>=0.7).mean()),
    "pct_fish>=0.9": pct((main>=0.9).mean()),
    "pct_fish==1.0": pct((main>=0.999).mean()),
    "mean_of_grouped_only": round(pf["prop_km_in_dominant_group_of_grouped"].mean(),4),
    "mean_km_outside_share": round(pf["prop_km_outside_regional_groups"].mean(),4),
    "mean_total_assigned_km": round(pf["total_assigned_km"].mean(),1),
    "mean_n_groups_covered": round(pf["n_groups_covered"].mean(),3),
    "mean_n_stream_edges_retained": round(pf["n_stream_edges_retained"].mean(),1),
}
pd.DataFrame([summary]).T.rename(columns={0:"value"}).to_csv(os.path.join(OUTDIR,"data","summary_statistics.csv"))
pf.groupby("dominant_group").agg(
    n_fish=("prop_km_in_dominant_group_of_all","size"),
    mean_prop=("prop_km_in_dominant_group_of_all","mean"),
    median_prop=("prop_km_in_dominant_group_of_all","median"),
).round(4).to_csv(os.path.join(OUTDIR,"data","by_dominant_group.csv"))
pf.groupby("year").agg(
    n_fish=("prop_km_in_dominant_group_of_all","size"),
    mean_prop=("prop_km_in_dominant_group_of_all","mean"),
    median_prop=("prop_km_in_dominant_group_of_all","median"),
    mean_total_km=("total_assigned_km","mean"),
).round(4).to_csv(os.path.join(OUTDIR,"data","by_year.csv"))
import json; print(json.dumps(summary, indent=2))
