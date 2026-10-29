"""Protein-level analysis v2: pair-level Fisher, full leave-one-drug-out (both arms),
key-drug exclusions, conflicting-attribution exclusion, and a multiple-membership
drug random-intercept logistic model (statsmodels BinomialBayesMixedGLM)."""
import sys, numpy as np, pandas as pd, scipy.sparse as sp
from scipy.stats import fisher_exact
from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM
D=sys.argv[1]; OUT=sys.argv[2]
e=pd.read_csv(D+'/enzyme_pair_severity.csv')
fix={'Q4U2R8':'SLC22A6 (OAT1)','Q8TCC7':'SLC22A8 (OAT3)','Q9Y6L6':'SLCO1B1 (OATP1B1)'}
e['gene']=[fix.get(u,g) for u,g in zip(e.uniprot,e.enzyme_gene)]
e['pair']=[tuple(sorted([a,b])) for a,b in zip(e.CID_A,e.CID_B)]
g=e[e.severity!='Unknown']
inh=g[g.direction.isin(['inhibition','transporter_inhibition'])]
ind=g[g.direction.isin(['induction','transporter_induction'])]
pe=inh[['pair','gene','severity']].drop_duplicates(['pair','gene'])
sev=pe.drop_duplicates('pair').set_index('pair').severity
allp=sorted(set(pe.pair)); maj={p:int(sev[p]=='Major') for p in allp}
cnt=pe.groupby('gene').pair.nunique(); tested=list(cnt[cnt>=10].sort_values(ascending=False).index)
names=dict(zip(e.CID_A,e.drug_A)); names.update(dict(zip(e.CID_B,e.drug_B)))
# conflicting: pair has both inhibition and induction attributed to same protein
conf=set(pe.merge(ind[['pair','gene']].drop_duplicates(),on=['pair','gene']).pair)
def tab(P,O):
    a=sum(maj[p] for p in P); b=len(P)-a; c=sum(maj[p] for p in O); d=len(O)-c; return a,b,c,d
def orci(a,b,c,d):
    h=0.5 if min(a,b,c,d)==0 else 0
    A,B,C,Dd=a+h,b+h,c+h,d+h; o=A*Dd/(B*C); se=np.sqrt(1/A+1/B+1/C+1/Dd)
    return o,np.exp(np.log(o)-1.96*se),np.exp(np.log(o)+1.96*se),fisher_exact([[a,b],[c,d]])[1]
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=p[o]*n/np.arange(1,n+1); q=np.minimum.accumulate(q[::-1])[::-1]; r=np.empty(n); r[o]=np.minimum(q,1); return r
drugs=sorted({d for p in allp for d in p}); di={d:i for i,d in enumerate(drugs)}
Z=sp.lil_matrix((len(allp),len(drugs)))
for i,p in enumerate(allp):
    for d in p: Z[i,di[d]]=1
Z=Z.toarray()
y=np.array([maj[p] for p in allp])
# top Major-contributing drugs overall
contrib=pd.Series({d:sum(maj[p] for p in allp if d in p) for d in drugs}).sort_values(ascending=False)
LEF=[d for d in drugs if 'Leflunomide' in str(names.get(d))][0]; MTX=[d for d in drugs if 'Methotrexate' in str(names.get(d))][0]
def re_fit(sub,ind_cols):
    sub=sorted(sub); dr=sorted({d for p in sub for d in p}); ix={d:i for i,d in enumerate(dr)}
    Zs=np.zeros((len(sub),len(dr)))
    for i,p in enumerate(sub):
        for d in p: Zs[i,ix[d]]=1
    yy=np.array([maj[p] for p in sub])
    X=np.column_stack([np.ones(len(sub))]+[[1 if p in c else 0 for p in sub] for c in ind_cols])
    f=BinomialBayesMixedGLM(yy,X,Zs,np.zeros(len(dr),dtype=int),vcp_p=1.0,fe_p=2.0).fit_vb()
    return [(np.exp(f.fe_mean[j]),np.exp(f.fe_mean[j]-1.96*f.fe_sd[j]),np.exp(f.fe_mean[j]+1.96*f.fe_sd[j])) for j in range(1,X.shape[1])]
psets={g:set(pe[pe.gene==g].pair) for g in tested}
ntested=pe[pe.gene.isin(tested)].groupby('pair').gene.nunique(); single=set(ntested[ntested==1].index)
joint=dict(zip(tested,re_fit(set(allp),[psets[g] for g in tested])))
rows=[]
for gname in tested:
    P=set(pe[pe.gene==gname].pair); O=set(allp)-P
    a,b,c,d=tab(P,O); o,lo,hi,p=orci(a,b,c,d)
    # full LODO over all drugs
    ors=[];sig=0;tot=0
    for dd in drugs:
        k={q for q in allp if dd not in q}; P2=P&k; O2=k-P
        if len(P2)<5: continue
        r=orci(*tab(P2,O2)); ors.append(r[0]); tot+=1; sig+=r[3]<0.05
    def excl(ds):
        k={q for q in allp if not any(x in q for x in ds)}; P2=P&k; O2=k-P
        if len(P2)<5: return (np.nan,)*4+(len(P2),)
        r=orci(*tab(P2,O2)); return r+(len(P2),)
    r_lef=excl([LEF]); r_lm=excl([LEF,MTX]); r_top5=excl(list(contrib.index[:5]))
    kc=[q for q in allp if q not in conf]; P3=P&set(kc); O3=set(kc)-P
    r_conf=orci(*tab(P3,O3))+(len(P3),) if len(P3)>=5 else (np.nan,)*4+(len(P3),)
    # multiple-membership random intercept model
    x=np.array([1 if q in P else 0 for q in allp])
    X=np.column_stack([np.ones(len(y)),x])
    m=BinomialBayesMixedGLM(y,X,Z,np.zeros(Z.shape[1],dtype=int),vcp_p=1.0,fe_p=2.0)
    f=m.fit_vb()
    b1=f.fe_mean[1]; s1=f.fe_sd[1]
    mm=(np.exp(b1),np.exp(b1-1.96*s1),np.exp(b1+1.96*s1))
    re_lef=re_fit({q for q in allp if LEF not in q},[P])[0]
    re_conf=re_fit(set(allp)-conf,[P])[0]
    Ps=P&single; Os=(set(allp)&single)-Ps; r_s=orci(*tab(Ps,Os))
    sub=pe[pe.gene==gname]; dc=pd.concat([sub.pair.map(lambda t:t[0]),sub.pair.map(lambda t:t[1])]).value_counts()
    top=dc.index[0]; k={q for q in allp if top not in q}; P4=P&k
    r_top=orci(*tab(P4,k-P)) if len(P4)>1 else (np.nan,)*4
    extra=dict(re_lef_OR=re_lef[0],re_lef_lo=re_lef[1],re_lef_hi=re_lef[2],re_conf_OR=re_conf[0],re_conf_lo=re_conf[1],re_conf_hi=re_conf[2],
        joint_OR=joint[gname][0],joint_lo=joint[gname][1],joint_hi=joint[gname][2],
        single_n=len(Ps),single_OR=r_s[0],single_lo=r_s[1],single_hi=r_s[2],single_p=r_s[3],
        topdrug=str(names.get(top)),topdrug_n=int(dc.iloc[0]),top_OR=r_top[0],top_lo=r_top[1],top_hi=r_top[2],top_p=r_top[3])
    rows.append(dict(**extra,gene=gname,n=len(P),major=a,pct=100*a/len(P),OR=o,lo=lo,hi=hi,p=p,
        lodo_min=min(ors),lodo_max=max(ors),lodo_sig=sig,lodo_tot=tot,
        lef_OR=r_lef[0],lef_lo=r_lef[1],lef_hi=r_lef[2],lef_p=r_lef[3],lef_n=r_lef[4],
        lm_OR=r_lm[0],lm_lo=r_lm[1],lm_hi=r_lm[2],lm_p=r_lm[3],lm_n=r_lm[4],
        top5_OR=r_top5[0],top5_p=r_top5[3],
        conf_OR=r_conf[0],conf_lo=r_conf[1],conf_hi=r_conf[2],conf_p=r_conf[3],conf_n=r_conf[4],
        mm_OR=mm[0],mm_lo=mm[1],mm_hi=mm[2],drug_sd=float(np.exp(f.vcp_mean[0]))))
T=pd.DataFrame(rows); T['q']=bh(T.p)
T.to_csv(OUT,index=False)
pd.set_option('display.width',260)
print(T[['gene','n','major','OR','lo','hi','p','q','lodo_min','lodo_max','lodo_sig','lodo_tot','lef_OR','lef_p','lm_OR','lm_lo','lm_hi','lm_p','top5_OR','top5_p','conf_n','conf_OR','conf_p','mm_OR','mm_lo','mm_hi','drug_sd']].round(3).to_string())
print('pairs',len(allp),'major',y.sum(),'drugs',len(drugs),'conflicting pairs',len(conf))
print('top5 Major contributors',[(names.get(d),int(contrib[d])) for d in contrib.index[:5]])
print('leflunomide Major',contrib[LEF],'methotrexate Major',contrib[MTX])
# null model drug SD
m0=BinomialBayesMixedGLM(y,np.ones((len(y),1)),Z,np.zeros(Z.shape[1],dtype=int),vcp_p=1.0,fe_p=2.0).fit_vb()
mtx=[q for q in allp if MTX in q]; Po=psets['SLC22A6 (OAT1)']
a1=sum(maj[q] for q in mtx if q in Po); n1=sum(1 for q in mtx if q in Po); c1=sum(maj[q] for q in mtx if q not in Po)
print('within methotrexate: OAT1',a1,'/',n1,'other',c1,'/',len(mtx)-n1,'Fisher p',fisher_exact([[a1,n1-a1],[c1,len(mtx)-n1-c1]])[1])
print('null drug SD',np.exp(m0.vcp_mean[0]))
