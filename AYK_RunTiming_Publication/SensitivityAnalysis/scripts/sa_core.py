"""Python replication of the Kuskokwim Chinook natal-origin assignment engine
(mirrors R/00_setup.R perform_assignment + aggregation). Uses only edge
attributes (no geometry), read from the shapefile .dbf via pyshp."""
import os
import numpy as np
import pandas as pd
import shapefile  # pyshp

PUB = "/sessions/nice-sweet-cori/mnt/AYK_RunTiming/AYK_RunTiming_Publication"
EDGES_SHP = os.path.join(PUB, "data", "spatial", "Spatial Data", "KuskoUSGS_HUC_joined.shp")
NATAL_DIR = os.path.join(PUB, "data", "natal_origins")

MIN_ERROR = 0.0006
MIN_STREAM_ORDER = 3
YEARS = [2017, 2018, 2019, 2020, 2021, 2022]

# Authoritative regional-group map from R/06_cluster_analysis.R.
# S. Fork Kusko -> Group 2 (combined w/ Upper Kusko Main); Johnson/Lower Kusko excluded.
CLUSTER_MAP = {
    "N. Fork Kusko": "Group 1", "E. Fork Kuskokwim River": "Group 1", "E. Fork Kuskokwim": "Group 1",
    "S. Fork Kusko": "Group 2", "Upper Kusko Main": "Group 2", "Big River": "Group 2",
    "Swift": "Group 2", "Tatlawiksuk": "Group 2",
    "Stony": "Group 3",
    "Takotna and Nixon Fork": "Group 4", "Holitna": "Group 4", "George": "Group 4", "Hoholitna": "Group 4",
    "Middle Kusko Main": "Group 5", "Oskakawlik": "Group 5", "Holokuk": "Group 5",
    "Kisaralik": "Group 5", "Tuluksak": "Group 5",
    "Aniak": "Group 6", "Kwethluk": "Group 6",
}
GROUP_ORDER = [f"Group {i}" for i in range(1, 7)]
Q_JUN11, Q_JUN21, Q_JUL01 = 163, 173, 183

def load_edges():
    r = shapefile.Reader(EDGES_SHP)
    flds = [f[0] for f in r.fields[1:]]
    idx = {n: i for i, n in enumerate(flds)}
    cols = ["iso_pred", "isose_pred", "Str_Order", "UniPh2oNoE", "SPAWNING_C", "mgmt_river", "Shape_Leng"]
    rows = [[rec[idx[c]] for c in cols] for rec in r.iterRecords()]
    df = pd.DataFrame(rows, columns=cols)
    for c in ["iso_pred", "isose_pred", "Str_Order", "UniPh2oNoE", "SPAWNING_C", "Shape_Leng"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["mgmt_river"] = df["mgmt_river"].astype(str).str.strip()
    df = df[df["Str_Order"] >= MIN_STREAM_ORDER].reset_index(drop=True)
    df.loc[df["mgmt_river"] == "Johnson", "mgmt_river"] = "Lower Kusko"
    df["length_km"] = df["Shape_Leng"] / 1000.0
    df["StreamOrderPrior"] = 1.0
    df["PresencePrior"] = np.where(df["Str_Order"].isin([6, 7, 8]) & (df["SPAWNING_C"] == 0), 0.0, 1.0)
    df["pid_prior"] = df["UniPh2oNoE"].astype(float)
    isose_mod = np.where(df["isose_pred"] < MIN_ERROR, MIN_ERROR, df["isose_pred"])
    within_site = 0.0003133684 / 1.96
    analyt = 0.00011 / 2
    df["error"] = np.sqrt(isose_mod**2 + within_site**2 + analyt**2)
    df["prior_prod"] = df["pid_prior"] * df["StreamOrderPrior"] * df["PresencePrior"]
    df["group"] = df["mgmt_river"].map(CLUSTER_MAP)
    return df

def load_natal(year):
    fp = os.path.join(NATAL_DIR, f"{year}_Kusko_Natal_Origins_Genetics_CPUE.csv")
    d = pd.read_csv(fp)
    d = d[d["natal_iso"].notna() & d["dailyCPUEprop"].notna()].reset_index(drop=True)
    return d

def rescaled_matrix(edges, iso_vals, threshold):
    iso_pred = edges["iso_pred"].to_numpy()[:, None]
    err = edges["error"].to_numpy()[:, None]
    prior = edges["prior_prod"].to_numpy()[:, None]
    iso_o = np.asarray(iso_vals, dtype=float)[None, :]
    coef = 1.0 / np.sqrt(2 * np.pi * err**2)
    assign = coef * np.exp(-((iso_o - iso_pred) ** 2) / (2 * err**2)) * prior
    assign_norm = assign / assign.sum(axis=0, keepdims=True)
    rescaled = assign_norm / assign_norm.max(axis=0, keepdims=True)
    rescaled[rescaled < threshold] = 0.0
    return rescaled

def basin_assign_sum(edges, natal_df, threshold):
    iso_vals = natal_df["natal_iso"].to_numpy(dtype=float)
    co = pd.to_numeric(natal_df["COratio"], errors="coerce").to_numpy()
    rescaled = rescaled_matrix(edges, iso_vals, threshold)
    co_safe = np.where(np.isnan(co), 0.0, co)
    return (rescaled * co_safe[None, :]).sum(axis=1)

def mgmt_production(edges, bas):
    tmp = edges.copy()
    tmp["basin_assign"] = bas
    tmp = tmp[(tmp["mgmt_river"] != "") & (tmp["mgmt_river"].notna()) & (tmp["mgmt_river"] != "nan")]
    g = tmp.groupby("mgmt_river", as_index=False)["basin_assign"].sum()
    return g.rename(columns={"basin_assign": "total_production"})

def quartile_of(doy):
    if doy <= Q_JUN11: return "Q1"
    if doy <= Q_JUN21: return "Q2"
    if doy <= Q_JUL01: return "Q3"
    return "Q4"
