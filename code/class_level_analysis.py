"""Class-level tests (Tables 1, 2, S4) from results/class_level_contingency.csv.
Usage: python code/class_level_analysis.py results"""
import sys, numpy as np, pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, random_table
from scipy.stats.contingency import odds_ratio
K=pd.read_csv(sys.argv[1]+'/class_level_contingency.csv').set_index('primary_mechanism')[['Minor','Moderate','Major']]
M=K.values; n=M.sum()
chi,p,dof,_=chi2_contingency(M,correction=False)
rt=random_table(M.sum(1),M.sum(0),seed=20261029); sims=rt.rvs(20000)
ex=sum((((t-np.outer(t.sum(1),t.sum(0))/n)**2)/(np.outer(t.sum(1),t.sum(0))/n)).sum()>=chi for t in sims)
print(f'chi2({dof}) = {chi:.2f}, asymptotic p = {p:.2e}, Monte Carlo exceedances = {ex}/20000, Cramer V = {np.sqrt(chi/(n*2)):.3f}')
tot=M.sum(0); rows={'MetabolismInhibition':('Major',),'Synergistic':('Moderate',),'ProteinBinding':('Minor','Moderate'),'MetabolismInduction':('Major',)}
def bh(p):
    p=np.asarray(p,float); m=len(p); o=np.argsort(p); q=p[o]*m/np.arange(1,m+1); q=np.minimum.accumulate(q[::-1])[::-1]; r=np.empty(m); r[o]=np.minimum(q,1); return r
res=[]
for name,outc in rows.items():
    r=K.loc[name]; a=sum(r[c] for c in outc); b=r.sum()-a; c=sum(tot[list(K.columns).index(o)] for o in outc)-a; d=n-r.sum()-c
    OR=a*d/(b*c); se=np.sqrt(1/a+1/b+1/c+1/d); pv=fisher_exact([[a,b],[c,d]])[1]
    rr=(a/(a+b))/(c/(c+d)); s2=np.sqrt(1/a-1/(a+b)+1/c-1/(c+d)); lo=np.exp(np.log(rr)-1.96*s2) if rr>1 else np.exp(np.log(rr)+1.96*s2)
    ev=lambda x:(x if x>=1 else 1/x)+np.sqrt((x if x>=1 else 1/x)*((x if x>=1 else 1/x)-1))
    res.append([name,'/'.join(outc),a,b,c,d,OR,np.exp(np.log(OR)-1.96*se),np.exp(np.log(OR)+1.96*se),pv,rr,ev(rr),ev(lo)])
# H5: transporter inhibition 8/8 Major (pooled in 'Other' in the contingency file)
a,b=8,0; c=tot[2]-8; d=n-8-c; r5=odds_ratio([[a,b],[c,d]],kind='conditional')
res.append(['TransporterInhibition','Major',a,b,c,d,np.inf,r5.confidence_interval()[0],np.inf,fisher_exact([[a,b],[c,d]])[1],np.nan,np.nan,np.nan])
T=pd.DataFrame(res,columns=['mechanism','outcome','a','b','c','d','OR','lo','hi','p','RR','E_point','E_limit']); T['q']=bh(T.p)
pd.set_option('display.width',200); print(T.round(4).to_string(index=False))
