import sys; sys.path.insert(0,'.')
import os, numpy as np, pandas as pd, sa_core as s
from cumdist import cumulative_distribution, closure_protection

THRESHOLDS=[0.5,0.6,0.7,0.8,0.9]
BASE=0.7
OUTDIR="/sessions/nice-sweet-cori/mnt/AYK_RunTiming/AYK_RunTiming_Publication/SensitivityAnalysis/2_Threshold_Sensitivity"
os.makedirs(os.path.join(OUTDIR,"data"),exist_ok=True)
os.makedirs(os.path.join(OUTDIR,"figures"),exist_ok=True)

edges=s.load_edges()

# common DOY grid for cumulative-distribution comparison
GRID=np.arange(149,206,3)

all_year=[]      # closure protection per year/group/threshold
all_summary=[]   # mean/min/max per group/threshold
cum_curves=[]    # group cumulative % per threshold on common grid (avg across years)

for thr in THRESHOLDS:
    cum=cumulative_distribution(edges,thr)
    # ---- closure protection ----
    clus=closure_protection(cum)
    clus["threshold"]=thr
    all_year.append(clus.copy())
    summ=clus.groupby("cluster",as_index=False).agg(
        mean_pct=("closure_protection_pct","mean"),
        min_pct=("closure_protection_pct","min"),
        max_pct=("closure_protection_pct","max"))
    summ["threshold"]=thr
    all_summary.append(summ)
    # ---- group cumulative distribution curves ----
    d=cum[~cum["mgmt_river"].isin(["Johnson","Lower Kusko"])].copy()
    d["cluster"]=d["mgmt_river"].map(s.CLUSTER_MAP)
    d=d[d["cluster"].notna()]
    grp=d.groupby(["year","cluster","doy"],as_index=False)["cumulative_production"].sum()
    for (yr,cl),g in grp.groupby(["year","cluster"]):
        g=g.sort_values("doy")
        final=g["cumulative_production"].max()
        pct=g["cumulative_production"].values/final*100
        yi=np.interp(GRID,g["doy"].values,pct,left=np.nan,right=100.0)
        for doy,val in zip(GRID,yi):
            cum_curves.append((thr,cl,yr,doy,val))

by_year=pd.concat(all_year,ignore_index=True)
by_year["cluster_num"]=by_year["cluster"].str.replace("Group ","").astype(int)
by_year=by_year.sort_values(["threshold","cluster_num","year"])
by_year[["threshold","year","cluster","closure_protection_pct","cluster_total_production","cluster_closure_production","n_units"]].to_csv(
    os.path.join(OUTDIR,"data","closure_protection_by_year_and_threshold.csv"),index=False)

summ=pd.concat(all_summary,ignore_index=True)
summ["cluster_num"]=summ["cluster"].str.replace("Group ","").astype(int)
# wide table: mean protection per group x threshold
wide=summ.pivot(index="cluster",columns="threshold",values="mean_pct").round(2)
wide.index=pd.CategoricalIndex(wide.index,categories=s.GROUP_ORDER,ordered=True)
wide=wide.sort_index()
wide.to_csv(os.path.join(OUTDIR,"data","closure_protection_mean_by_group_x_threshold.csv"))
summ.sort_values(["cluster_num","threshold"])[["cluster","threshold","mean_pct","min_pct","max_pct"]].round(3).to_csv(
    os.path.join(OUTDIR,"data","closure_protection_summary_by_group_and_threshold.csv"),index=False)

cc=pd.DataFrame(cum_curves,columns=["threshold","cluster","year","doy","cum_pct"])
cc_avg=cc.groupby(["threshold","cluster","doy"],as_index=False)["cum_pct"].mean()
cc_avg.to_csv(os.path.join(OUTDIR,"data","group_cumulative_distribution_by_threshold.csv"),index=False)

print("=== Mean closure protection (%) by group x threshold (0.7 = baseline) ===")
print(wide.to_string())
# how much does baseline change?
base=summ[summ.threshold==BASE].set_index("cluster")["mean_pct"]
print("\n=== Max deviation from baseline (0.7) across thresholds, by group ===")
for cl in s.GROUP_ORDER:
    vals=summ[summ.cluster==cl].set_index("threshold")["mean_pct"]
    dev=(vals-base[cl]).abs().max()
    print(f"  {cl}: baseline={base[cl]:.1f}%  max abs dev={dev:.1f} pts  range=[{vals.min():.1f},{vals.max():.1f}]")
# ranking stability
print("\n=== Group ranking (most->least protected) at each threshold ===")
for thr in THRESHOLDS:
    r=summ[summ.threshold==thr].sort_values("mean_pct",ascending=False)["cluster"].tolist()
    print(f"  thr {thr}: {' > '.join(r)}")
