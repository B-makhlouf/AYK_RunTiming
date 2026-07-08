import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR="/sessions/nice-sweet-cori/mnt/AYK_RunTiming/AYK_RunTiming_Publication/SensitivityAnalysis/2_Threshold_Sensitivity"
GROUPS=[f"Group {i}" for i in range(1,7)]
gcolors=dict(zip(GROUPS, plt.cm.viridis(np.linspace(0,0.92,6))))
THRESHOLDS=[0.5,0.6,0.7,0.8,0.9]
tcolors=dict(zip(THRESHOLDS, plt.cm.plasma(np.linspace(0,0.85,5))))

summ=pd.read_csv(os.path.join(OUTDIR,"data","closure_protection_summary_by_group_and_threshold.csv"))

# ---- FIG 1: closure protection vs threshold (headline) ----
fig,ax=plt.subplots(figsize=(9,6))
for g in GROUPS:
    d=summ[summ.cluster==g].sort_values("threshold")
    ax.plot(d.threshold,d.mean_pct,marker="o",lw=2.2,color=gcolors[g],label=g)
    ax.fill_between(d.threshold,d.min_pct,d.max_pct,color=gcolors[g],alpha=0.10)
ax.axvline(0.7,color="grey",ls="--",lw=1.5)
ax.text(0.7,ax.get_ylim()[1]*0.97,"  paper baseline (0.7)",color="grey",va="top",fontsize=9)
ax.set_xlabel("Assignment threshold")
ax.set_ylabel("Mean front-end closure protection (%)\n(June 1–11, avg. across 2017–2022)")
ax.set_title("Closure protection by regional group is robust to the assignment threshold",fontsize=12,fontweight="bold")
ax.set_xticks(THRESHOLDS)
ax.legend(frameon=False,title="Regional group",ncol=2)
ax.spines[["top","right"]].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(OUTDIR,"figures","closure_protection_vs_threshold.png"),dpi=200)
plt.close(fig)

# ---- FIG 2: same info as grouped bars (very easy to read) ----
fig,ax=plt.subplots(figsize=(10,6))
x=np.arange(len(GROUPS)); w=0.16
for i,thr in enumerate(THRESHOLDS):
    vals=[summ[(summ.cluster==g)&(summ.threshold==thr)]["mean_pct"].iloc[0] for g in GROUPS]
    bars=ax.bar(x+(i-2)*w,vals,w,color=tcolors[thr],
                label=f"{thr}"+("  (baseline)" if thr==0.7 else ""),
                edgecolor="black" if thr==0.7 else "none",linewidth=1.2 if thr==0.7 else 0)
ax.set_xticks(x); ax.set_xticklabels(GROUPS)
ax.set_ylabel("Mean front-end closure protection (%)")
ax.set_xlabel("Regional group")
ax.set_title("Front-end closure protection by group across assignment thresholds",fontsize=12,fontweight="bold")
ax.legend(frameon=False,title="Threshold",ncol=5,loc="upper right")
ax.spines[["top","right"]].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(OUTDIR,"figures","closure_protection_bars_by_threshold.png"),dpi=200)
plt.close(fig)

# ---- FIG 3: group cumulative distribution curves, faceted, threshold lines ----
cc=pd.read_csv(os.path.join(OUTDIR,"data","group_cumulative_distribution_by_threshold.csv"))
fig,axes=plt.subplots(2,3,figsize=(14,8),sharex=True,sharey=True)
for ax,g in zip(axes.ravel(),GROUPS):
    for thr in THRESHOLDS:
        d=cc[(cc.cluster==g)&(cc.threshold==thr)].sort_values("doy")
        ax.plot(d.doy,d.cum_pct,lw=2 if thr==0.7 else 1.3,color=tcolors[thr],
                ls="-" if thr==0.7 else "-",alpha=1 if thr==0.7 else 0.85,
                label=f"{thr}"+(" (baseline)" if thr==0.7 else ""))
    ax.axvspan(153,162,color="grey",alpha=0.15)
    ax.axhline(50,color="grey",ls=":",lw=0.8)
    ax.set_title(g,fontsize=11)
    ax.spines[["top","right"]].set_visible(False)
for ax in axes[-1]: ax.set_xlabel("Day of year")
for ax in axes[:,0]: ax.set_ylabel("Cumulative % of group's run")
axes[0,2].legend(frameon=False,title="Threshold",fontsize=8)
fig.suptitle("Group cumulative run-timing distributions across assignment thresholds\n(shaded band = June 1–11 front-end closure window)",
             fontsize=12,fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig(os.path.join(OUTDIR,"figures","group_cumulative_distribution_by_threshold.png"),dpi=200)
plt.close(fig)
print("part2 figures written")
