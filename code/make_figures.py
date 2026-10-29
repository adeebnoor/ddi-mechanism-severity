# Usage: python make_figures.py results figures
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys, os
RES=sys.argv[1] if len(sys.argv)>1 else 'results'
O=(sys.argv[2] if len(sys.argv)>2 else 'figures')+'/'
os.makedirs(O,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
C={'Minor':'#9ecae1','Moderate':'#4292c6','Major':'#08306b'}
# Fig1
T={'Synergistic':(9,278,259),'Metabolism inhibition':(1,52,198),'Metabolism induction':(3,52,54),'Protein binding':(0,6,20),'Additive':(2,15,1),'Other*':(0,4,9)}
names=list(T)[::-1]; M=np.array([T[k] for k in names],float); n=M.sum(1); P=M/n[:,None]*100
fig,ax=plt.subplots(figsize=(6.3,2.9)); left=np.zeros(len(names))
for j,lab in enumerate(['Minor','Moderate','Major']):
    ax.barh(names,P[:,j],left=left,color=C[lab],label=lab,edgecolor='white',linewidth=.5); left+=P[:,j]
for i,k in enumerate(names): ax.text(101,i,f'n = {int(n[i])}',va='center',fontsize=8)
ax.set_xlim(0,115); ax.set_xticks([0,25,50,75,100]); ax.set_xlabel('Share of pairs (%)')
ax.legend(ncol=3,frameon=False,loc='upper center',bbox_to_anchor=(.45,1.16)); fig.tight_layout(); fig.savefig(O+'Figure1.png',dpi=300); fig.savefig(O+'Figure1.tiff',dpi=300)
# Fig2 forest
E=pd.read_csv(RES+'/protein_level_results.csv').sort_values('OR')
fig,ax=plt.subplots(figsize=(6.3,3.9)); y=np.arange(len(E))
ax.axvline(1,color='grey',lw=.8,ls='--')
ax.errorbar(E.OR,y,xerr=[E.OR-E.lo,E.hi-E.OR],fmt='o',color='#08306b',ms=5,capsize=2,label='Primary analysis')
ax.errorbar(E.single_OR,y+0.22,xerr=[E.single_OR-E.single_lo,E.single_hi-E.single_OR],fmt='s',color='#6baed6',ms=3.5,capsize=1.5,lw=.8,label='Single-protein pairs only')
m=E.lodo_OR.notna()
ax.errorbar(E.lodo_OR[m],y[m.values]-0.22,xerr=[(E.lodo_OR-E.lodo_lo)[m],(E.lodo_hi-E.lodo_OR)[m]],fmt='^',color='#d94801',ms=3.5,capsize=1.5,lw=.8,label='Most frequent drug removed')
ax.set_xscale('log'); ax.set_yticks(y); ax.set_yticklabels([f'{g} (n = {k})' for g,k in zip(E.gene,E.n)])
ax.set_xlabel('Odds ratio for a Major grade (log scale)'); ax.set_xlim(0.015,30)
from matplotlib.ticker import FixedLocator, FixedFormatter
ax.xaxis.set_major_locator(FixedLocator([0.02,0.05,0.1,0.2,0.5,1,2,5,10,20])); ax.xaxis.set_major_formatter(FixedFormatter(['0.02','0.05','0.1','0.2','0.5','1','2','5','10','20'])); ax.xaxis.set_minor_locator(FixedLocator([]))
ax.legend(frameon=False,fontsize=7.5,loc='upper center',ncol=3,bbox_to_anchor=(0.4,1.12)); fig.tight_layout(); fig.savefig(O+'Figure2.png',dpi=300); fig.savefig(O+'Figure2.tiff',dpi=300)
# Fig3 external validation
dbk={'Major':(11,44),'Moderate':(37,246),'Minor':(2,26),'Unknown':(23,856)}
kg={'Major':(20,57),'Moderate':(95,268),'Minor':(4,28),'Unknown':(85,890)}
fig,axs=plt.subplots(1,2,figsize=(6.3,2.6),sharey=False)
for ax,(t,d) in zip(axs,[('(a) DrugBank clinically proven set',dbk),('(b) KEGG DDI set',kg)]):
    labs=list(d); v=[100*a/b for a,b in d.values()]
    ax.bar(labs,v,color=['#08306b','#4292c6','#9ecae1','#bdbdbd'])
    for i,(a,b) in enumerate(d.values()): ax.text(i,100*a/b+.8,f'{a}/{b}',ha='center',fontsize=7.5)
    ax.set_title(t,fontsize=9,loc='left'); ax.set_ylabel('Pairs confirmed (%)'); ax.set_ylim(0,max(v)*1.25)
fig.tight_layout(); fig.savefig(O+'Figure3.png',dpi=300); fig.savefig(O+'Figure3.tiff',dpi=300)
print('ok')
