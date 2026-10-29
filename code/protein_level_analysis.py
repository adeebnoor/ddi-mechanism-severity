import pandas as pd, numpy as np, sys, json
from scipy.stats import fisher_exact
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); r=np.empty(n); q=p[o]*n/np.arange(1,n+1)
    q=np.minimum.accumulate(q[::-1])[::-1]; r[o]=np.minimum(q,1); return r
D=sys.argv[1]
e=pd.read_csv(D+'/enzyme_pair_severity.csv')
# relabel by UniProt accession (Swiss-Prot verified)
fix={'Q4U2R8':'SLC22A6 (OAT1)','Q8TCC7':'SLC22A8 (OAT3)','Q9Y6L6':'SLCO1B1 (OATP1B1)'}
e['gene']=[fix.get(u,g) for u,g in zip(e.uniprot,e.enzyme_gene)]
e['pair']=[tuple(sorted([a,b])) for a,b in zip(e.CID_A,e.CID_B)]
inh=e[e.direction.isin(['inhibition','transporter_inhibition'])&(e.severity!='Unknown')]
pe=inh[['pair','gene','severity','drug_A','drug_B','CID_A','CID_B']].drop_duplicates(['pair','gene'])
sev=pe.drop_duplicates('pair').set_index('pair').severity
allp=set(pe.pair)
cnt=pe.groupby('gene').pair.nunique().sort_values(ascending=False)
tested=list(cnt[cnt>=10].index)
print('eligible',cnt[cnt>=10].to_dict()); print('unique non-Unknown inhibition pairs',len(allp))
ntp=len(set(pe[pe.gene.isin(tested)].pair)); print('pairs with >=1 tested enzyme',ntp)
def stats(P,O):
    a=sum(sev[p]=='Major' for p in P); b=len(P)-a; c=sum(sev[p]=='Major' for p in O); d=len(O)-c
    t=[[a,b],[c,d]]; p=fisher_exact(t)[1]
    if min(a,b,c,d)==0: aa,bb,cc,dd=[x+0.5 for x in (a,b,c,d)]
    else: aa,bb,cc,dd=a,b,c,d
    OR=aa*dd/(bb*cc); se=np.sqrt(1/aa+1/bb+1/cc+1/dd)
    return dict(n=len(P),major=a,pct=100*a/len(P),ncomp=len(O),majcomp=c,OR=OR,lo=np.exp(np.log(OR)-1.96*se),hi=np.exp(np.log(OR)+1.96*se),p=p)
def evalue(OR,lo,hi):
    # E-value for OR (rare-outcome not assumed: use sqrt transform for common outcome)
    def ev(r):
        if r<1: r=1/r
        return r+np.sqrt(r*(r-1))
    r=np.sqrt(OR) if True else OR
    ci = np.sqrt(lo) if OR>1 else np.sqrt(hi)
    e_pt=ev(r); e_ci = 1.0 if (lo<=1<=hi) else ev(ci)
    return e_pt,e_ci
rows=[]
for g in tested:
    P=set(pe[pe.gene==g].pair); O=allp-P
    s=stats(P,O); s['gene']=g; rows.append(s)
T=pd.DataFrame(rows); T['q']=bh(T.p)
ev=[evalue(r.OR,r.lo,r.hi) for r in T.itertuples()]
T['E_pt']=[x[0] for x in ev]; T['E_ci']=[x[1] for x in ev]
# sensitivity 1: single-enzyme pairs (pairs attributed to exactly one tested enzyme)
ntested=pe[pe.gene.isin(tested)].groupby('pair').gene.nunique()
single=set(ntested[ntested==1].index)
s1=[]
for g in tested:
    P=set(pe[pe.gene==g].pair)&single; O=(allp&single)-P
    s1.append(stats(P,O) if len(P)>0 else None)
T['single_n']=[x['n'] for x in s1]; T['single_OR']=[x['OR'] for x in s1]; T['single_lo']=[x['lo'] for x in s1]; T['single_hi']=[x['hi'] for x in s1]; T['single_p']=[x['p'] for x in s1]
# sensitivity 2: drop the most frequent drug within enzyme set (both arms)
s2=[];top=[]
for g in tested:
    sub=pe[pe.gene==g]
    drugs=pd.concat([sub.CID_A,sub.CID_B]).value_counts()
    cid=drugs.index[0]; share=drugs.iloc[0]/sub.pair.nunique()
    name=pd.concat([sub[sub.CID_A==cid].drug_A,sub[sub.CID_B==cid].drug_B]).iloc[0]
    keep={p for p in allp if cid not in p}
    P=set(sub.pair)&keep; O=keep-set(sub.pair)
    s2.append(stats(P,O) if len(P)>1 else None); top.append((name,int(drugs.iloc[0]),round(100*share,1)))
T['topdrug']=[t[0] for t in top]; T['topdrug_n']=[t[1] for t in top]; T['topdrug_pct']=[t[2] for t in top]
T['lodo_n']=[x['n'] if x else None for x in s2]; T['lodo_OR']=[x['OR'] if x else None for x in s2]; T['lodo_lo']=[x['lo'] if x else None for x in s2]; T['lodo_hi']=[x['hi'] if x else None for x in s2]; T['lodo_p']=[x['p'] if x else None for x in s2]
T=T.sort_values('OR',ascending=False)
pd.set_option('display.width',250)
print(T.round(4).to_string())
T.to_csv(sys.argv[2],index=False)
print('single-enzyme pairs',len(allp&single))
# baseline
print('overall major among all non-Unknown inhibition pairs', sum(sev[p]=='Major' for p in allp), len(allp))
