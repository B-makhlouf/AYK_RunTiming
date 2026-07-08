import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR="/sessions/nice-sweet-cori/mnt/AYK_RunTiming/AYK_RunTiming_Publication/SensitivityAnalysis/1_Individual_Assignment_Concentration"
pf=pd.read_csv(os.path.join(OUTDIR,"data","per_individual_group_concentration.csv"))
x=pf["prop_in_dominant_group_of_all"].dropna().values
mean_v, med_v = x.mean(), np.median(x)

fig,axes=plt.subplots(1,2,figsize=(13,5.2))

# Panel A: histogram
ax=axes[0]
ax.hist(x*100, bins=np.arange(30,101,5), color="#4C78A8", edgecolor="white")
ax.axvline(mean_v*100,color="#E45756",lw=2,ls="-",label=f"Mean = {mean_v*100:.1f}%")
ax.axvline(med_v*100,color="#222222",lw=2,ls="--",label=f"Median = {med_v*100:.1f}%")
ax.set_xlabel("% of an individual's assignment in its single dominant regional group")
ax.set_ylabel("Number of individual fish")
ax.set_title("A. Distribution across all 1,465 individuals",fontsize=12,loc="left")
ax.legend(frameon=False)
ax.spines[["top","right"]].set_visible(False)

# Panel B: ECDF with thresholds annotated
ax=axes[1]
xs=np.sort(x)*100
ys=np.arange(1,len(xs)+1)/len(xs)*100
ax.plot(xs,ys,color="#4C78A8",lw=2)
for t in [50,70,90]:
    frac=(x*100>=t).mean()*100
    ax.axvline(t,color="grey",ls=":",lw=1)
    ax.annotate(f"{frac:.0f}% of fish\n≥ {t}%",xy=(t,100-frac),xytext=(t+1,max(3,100-frac-6)),
                fontsize=9,color="#333333")
ax.set_xlabel("% of an individual's assignment in its single dominant regional group")
ax.set_ylabel("Cumulative % of individuals")
ax.set_title("B. Cumulative distribution",fontsize=12,loc="left")
ax.spines[["top","right"]].set_visible(False)

fig.suptitle("Individual assignment concentration in a single regional group (threshold = 0.7)",
             fontsize=13,fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.96])
fig.savefig(os.path.join(OUTDIR,"figures","individual_concentration_distribution.png"),dpi=200)
plt.close(fig)

# Panel: by dominant group boxplot
fig,ax=plt.subplots(figsize=(8,5))
groups=sorted(pf["dominant_group"].dropna().unique())
data=[pf.loc[pf["dominant_group"]==g,"prop_in_dominant_group_of_all"].dropna()*100 for g in groups]
bp=ax.boxplot(data,labels=groups,patch_artist=True,medianprops=dict(color="black",lw=2))
for patch in bp["boxes"]: patch.set_facecolor("#8FBBD9")
ax.axhline(mean_v*100,color="#E45756",ls="--",lw=1.5,label=f"Overall mean = {mean_v*100:.1f}%")
ax.set_ylabel("% of assignment in dominant group")
ax.set_xlabel("Dominant regional group")
ax.set_title("Concentration by dominant regional group (threshold = 0.7)",fontsize=12)
ax.legend(frameon=False); ax.spines[["top","right"]].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(OUTDIR,"figures","concentration_by_group.png"),dpi=200)
plt.close(fig)
print("figures written")
print(f"mean={mean_v:.4f} median={med_v:.4f}")
for t in [50,70,90]: print(f">= {t}%: {(x*100>=t).mean()*100:.1f}% of fish")
