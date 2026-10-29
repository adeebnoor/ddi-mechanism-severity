"""External checks (Table 5, Figure 3) from results/external_validation_by_grade.csv.
Usage: python code/external_checks.py results"""
import sys, numpy as np, pandas as pd
from scipy.stats import norm, fisher_exact
V=pd.read_csv(sys.argv[1]+'/external_validation_by_grade.csv').set_index('grade')
def wilson(k,n,z=1.96):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return 100*(c-h),100*(c+h)
def ca(k,n,s):
    k,n,s=map(lambda x:np.asarray(x,float),(k,n,s)); N=n.sum(); p=k.sum()/N; sb=(n*s).sum()/N
    z=(k*(s-sb)).sum()/np.sqrt(p*(1-p)*(n*(s-sb)**2).sum()); return z,2*norm.sf(abs(z))
for src in ['drugbank','kegg']:
    k=V[src+'_confirmed']; n=V[src+'_mapped']
    for g in V.index: print(src,g,f'{k[g]}/{n[g]}',f'{100*k[g]/n[g]:.1f}%',tuple(round(x,1) for x in wilson(k[g],n[g])))
    G=['Major','Moderate','Minor']; kg=k[G].sum(); ng=n[G].sum()
    print(src,'graded vs ungraded',fisher_exact([[kg,ng-kg],[k['Unknown'],n['Unknown']-k['Unknown']]]))
    print(src,'trend among graded (scores 3,2,1)',ca(k[G],n[G],[3,2,1]))
    print(src,'trend incl. ungraded (scores 4..1)',ca(k[G+['Unknown']],n[G+['Unknown']],[4,3,2,1]))
