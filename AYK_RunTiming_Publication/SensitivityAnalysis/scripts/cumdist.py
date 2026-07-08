import sys; sys.path.insert(0,'.')
import numpy as np, pandas as pd, sa_core as s

def cumulative_distribution(edges, threshold, interval_days=3):
    """Replicate R/03: per (year, doy, mgmt_river) cumulative_production."""
    out = []
    for yr in s.YEARS:
        nat = s.load_natal(yr)
        dmin, dmax = np.floor(nat["DOY"].min()), np.ceil(nat["DOY"].max())
        seq = list(np.arange(dmin, dmax+1e-9, interval_days))
        if seq[-1] < dmax:
            seq.append(dmax)
        for i, doy in enumerate(seq, start=1):
            sub = nat[nat["DOY"] <= doy]
            if len(sub)==0: continue
            bas = s.basin_assign_sum(edges, sub, threshold)
            mp = s.mgmt_production(edges, bas)
            mp["year"]=yr; mp["doy"]=doy
            out.append(mp.rename(columns={"total_production":"cumulative_production"}))
    return pd.concat(out, ignore_index=True)

def closure_protection(cum, start=153, end=162, mindoy=149, maxdoy=205):
    """Replicate R/06 Analysis 4 cluster closure protection by year."""
    d = cum.copy()
    d = d[~d["mgmt_river"].isin(["Johnson","Lower Kusko"])]
    d["cluster"] = d["mgmt_river"].map(s.CLUSTER_MAP)
    d = d[d["cluster"].notna()]
    d = d[(d["doy"]>=mindoy)&(d["doy"]<=maxdoy)]
    recs=[]
    for (yr,mgmt,cl), g in d.groupby(["year","mgmt_river","cluster"]):
        g=g.sort_values("doy")
        total = g["cumulative_production"].max()
        cw = g[(g["doy"]>=start)&(g["doy"]<=end)]
        cwp = g[g["doy"]==cw["doy"].max()]["cumulative_production"].iloc[0] if len(cw)>0 else 0.0
        recs.append((yr,mgmt,cl,total,cwp))
    unit=pd.DataFrame(recs,columns=["year","mgmt_river","cluster","total_annual_production","closure_window_production"])
    clus=unit.groupby(["year","cluster"],as_index=False).agg(
        cluster_total_production=("total_annual_production","sum"),
        cluster_closure_production=("closure_window_production","sum"),
        n_units=("mgmt_river","count"))
    clus["closure_protection_pct"]=clus["cluster_closure_production"]/clus["cluster_total_production"]*100
    return clus

if __name__=="__main__":
    edges=s.load_edges()
    cum=cumulative_distribution(edges,0.7)
    clus=closure_protection(cum)
    ref=pd.read_csv(f"{s.PUB}/outputs/Cluster Analysis Results/CSV_Exports/analysis4_cluster_closure_protection_by_year.csv")
    ref["cluster"]=ref["cluster"].str.replace("Cluster","Group")
    m=ref.merge(clus,on=["year","cluster"],suffixes=("_ref","_mine"))
    m["diff_pct"]=(m["closure_protection_pct_ref"]-m["closure_protection_pct_mine"]).abs()
    m["diff_tot"]=(m["cluster_total_production_ref"]-m["cluster_total_production_mine"]).abs()
    print("VALIDATION B: closure protection by year (threshold 0.7)")
    print("  matched rows:",len(m),"of ref",len(ref))
    print("  max abs diff protection_pct:",m["diff_pct"].max())
    print("  max abs diff cluster_total_production:",m["diff_tot"].max())
    print(m.sort_values("diff_pct",ascending=False).head(5)[["year","cluster","closure_protection_pct_ref","closure_protection_pct_mine","diff_pct"]].to_string(index=False))
