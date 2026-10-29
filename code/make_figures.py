import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter
RES=sys.argv[1]; O=sys.argv[2]+'/'; os.makedirs(O,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
C={'Minor':'#9ecae1','Moderate':'#4292c6','Major':'#08306b'}
# Figure 1 (from class-level counts)
K=pd.read_csv(RES+'/class_level_contingency.csv')
lab={'Synergistic':'Synergistic','MetabolismInhibition':'Metabolism inhibition','MetabolismInduction':'Metabolism induction','ProteinBinding':'Protein binding','Additive':'Additive'}
names=[lab.get(x,'Other*') for x in K.primary_mechanism][::-1]; M=K[['Minor','Moderate','Major']].values[::-1].astype(float); n=M.sum(1); P=M/n[:,None]*100
fig,ax=plt.subplots(figsize=(6.3,2.9)); left=np.zeros(len(names))
for j,l in enumerate(['Minor','Moderate','Major']):
    ax.barh(names,P[:,j],left=left,color=C[l],label=l,edgecolor='white',linewidth=.5); left+=P[:,j]
for i in range(len(names)): ax.text(101,i,f'n = {int(n[i])}',va='center',fontsize=8)
ax.set_xlim(0,115); ax.set_xticks([0,25,50,75,100]); ax.set_xlabel('Share of pairs (%)')
ax.legend(ncol=3,frameon=False,loc='upper center',bbox_to_anchor=(.45,1.16)); fig.tight_layout(); [fig.savefig(O+'Figure1.'+e,dpi=300) for e in ('png','tiff')]
# Figure 2
E=pd.read_csv(RES+'/protein_level_results.csv').sort_values('OR')
fig,ax=plt.subplots(figsize=(6.3,4.4)); y=np.arange(len(E)); ax.axvline(1,color='grey',lw=.8,ls='--')
spec=[('OR','lo','hi','o','#08306b',0.27,'Pair-level'),('lef_OR','lef_lo','lef_hi','^','#d94801',0.09,'Without leflunomide'),('conf_OR','conf_lo','conf_hi','s','#6baed6',-0.09,'Conflicting pairs excluded'),('mm_OR','mm_lo','mm_hi','D','#31a354',-0.27,'Drug random intercept')]
for o,l,h,mk,col,dy,lb in spec:
    ax.errorbar(E[o],y+dy,xerr=[E[o]-E[l],E[h]-E[o]],fmt=mk,color=col,ms=4 if mk=='o' else 3.3,capsize=1.5,lw=.9,label=lb)
ax.set_xscale('log'); ax.set_yticks(y); ax.set_yticklabels([f'{g} (n = {k})' for g,k in zip(E.gene,E.n)])
ax.set_xlim(0.04,30); ax.xaxis.set_major_locator(FixedLocator([0.05,0.1,0.2,0.5,1,2,5,10,20])); ax.xaxis.set_major_formatter(FixedFormatter(['0.05','0.1','0.2','0.5','1','2','5','10','20'])); ax.xaxis.set_minor_locator(FixedLocator([]))
ax.set_xlabel('Odds ratio for a Major grade (log scale)')
ax.legend(frameon=False,fontsize=7.5,loc='upper center',ncol=2,bbox_to_anchor=(0.42,1.15)); fig.tight_layout(); [fig.savefig(O+'Figure2.'+e,dpi=300) for e in ('png','tiff')]
# Figure 3
V=pd.read_csv(RES+'/external_validation_by_grade.csv')
def wil(k,n,z=1.96):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return 100*(c-h),100*(c+h)
fig,axs=plt.subplots(1,2,figsize=(6.3,2.7))
for ax,(t,kc,nc) in zip(axs,[('(a) DrugBank clinically proven set','drugbank_confirmed','drugbank_mapped'),('(b) KEGG DDI set','kegg_confirmed','kegg_mapped')]):
    v=100*V[kc]/V[nc]; w=[wil(k,n) for k,n in zip(V[kc],V[nc])]
    ax.bar(V.grade,v,color=['#08306b','#4292c6','#9ecae1','#bdbdbd'],yerr=[v-[a for a,b in w],[b for a,b in w]-v],capsize=3,error_kw={'lw':.8})
    for i,(k,n) in enumerate(zip(V[kc],V[nc])): ax.text(i,w[i][1]+1,f'{k}/{n}',ha='center',fontsize=7.5)
    ax.set_title(t,fontsize=9,loc='left'); ax.set_ylabel('Pairs confirmed (%)'); ax.set_ylim(0,max(b for a,b in w)*1.18)
fig.tight_layout(); [fig.savefig(O+'Figure3.'+e,dpi=300) for e in ('png','tiff')]
print('ok')
