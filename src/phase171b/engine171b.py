from common171b import *
import numpy as np, math
from collections import defaultdict
sys.path.insert(0,str(O/'runtime_dependencies'))
from numba import njit, prange

def metadata():
    rows=[r for r in csvread(A/'discovery_validation_assignment_preview.csv') if r['metadata_eligible']=='True']
    genes=csvread(A/'continuous_signature_definition.csv')
    pids=sorted({r['pert_id'] for r in rows})
    queries={'G2_R2_CAP8':np.array([float(r['main_R2_L2_weight']) for r in genes]),'G2_R1':np.array([float(r['sensitivity_R1_L2_weight']) for r in genes])}
    gi={r['HGNC_ID']:i for i,r in enumerate(genes)}
    for filename in ['discrete_signature_membership_preview.csv','pathway_signature_membership_preview.csv']:
        for r in csvread(A/filename):
            q=r['query_id']
            if q not in queries:queries[q]=np.zeros(len(genes),dtype=np.float64)
            queries[q][gi[r['HGNC_ID']]]=float(r['weight'])
    strata=[np.array([i for i,g in enumerate(genes) if g['null_stratum']==s],dtype=np.int64) for s in sorted({g['null_stratum'] for g in genes})]
    return rows,genes,pids,queries,strata

def csr(groups):
    sizes=np.array([len(g) for g in groups],dtype=np.int64)
    return np.array([i for g in groups for i in g],dtype=np.int64),np.r_[0,np.cumsum(sizes)]

def hierarchy(rows,pids):
    ck=[(r['pert_id'],r['fold'],r['base_cell_id'],r['cell_id'],float(r['dose_uM']),float(r['time_h'])) for r in rows]
    keys=sorted(set(ck)); index={k:i for i,k in enumerate(keys)}
    groups=[[] for _ in keys]
    for i,k in enumerate(ck):groups[index[k]].append(i)
    bk=sorted({k[:3] for k in keys});bi={k:i for i,k in enumerate(bk)}
    bg=[[] for _ in bk]
    for i,k in enumerate(keys):bg[bi[k[:3]]].append(i)
    pk=[(p,f) for f in ['DISCOVERY','VALIDATION'] for p in pids];pi={k:i for i,k in enumerate(pk)}
    pg=[[] for _ in pk]
    for i,k in enumerate(bk):pg[pi[k[:2]]].append(i)
    return dict(context_keys=keys,base_keys=bk,compound_keys=pk,profile_to_context=csr(groups),context_to_base=csr(bg),base_to_compound=csr(pg),context_groups=groups,base_groups=bg,compound_groups=pg)

@njit(cache=False)
def grouped_median(values,index,ptr):
    out=np.full((len(ptr)-1,values.shape[1]),np.nan)
    for g in range(len(ptr)-1):
        n=ptr[g+1]-ptr[g]
        if n==1:
            out[g,:]=values[index[ptr[g]],:]
            continue
        for j in range(values.shape[1]):
            v=np.empty(n,dtype=np.float64);k=0
            for a in range(ptr[g],ptr[g+1]):
                x=values[index[a],j]
                if np.isfinite(x):v[k]=x;k+=1
            if k:
                v=np.sort(v[:k]);m=k//2
                out[g,j]=v[m] if k%2 else (v[m-1]+v[m])/2
    return out

def aggregate(values,h):
    c=grouped_median(values,*h['profile_to_context'])
    b=grouped_median(c,*h['context_to_base'])
    p=grouped_median(b,*h['base_to_compound'])
    return c,b,p

def scalar_score(z,w):
    # Independent scalar product and norms, math.fsum has no BLAS dependency.
    ids=[i for i,x in enumerate(w) if x!=0]
    if any(not math.isfinite(float(z[i])) for i in ids):return float('nan')
    numerator=math.fsum(float(z[i])*float(w[i]) for i in ids)
    zn=math.sqrt(math.fsum(float(z[i])*float(z[i]) for i in ids))
    wn=math.sqrt(math.fsum(float(w[i])*float(w[i]) for i in ids))
    return numerator/(zn*wn) if zn and wn else float('nan')

def independent_median(values):
    x=sorted(float(v) for v in values if math.isfinite(float(v)))
    if not x:return float('nan')
    n=len(x);return x[n//2] if n%2 else (x[n//2-1]+x[n//2])/2

def independent_aggregate(row_score_pairs):
    groups=defaultdict(list)
    for r,s in row_score_pairs:groups[(r['pert_id'],r['fold'],r['base_cell_id'],r['cell_id'],float(r['dose_uM']),float(r['time_h']))].append(s)
    contexts={k:independent_median(v) for k,v in groups.items()}
    bases=defaultdict(list)
    for k,v in contexts.items():bases[k[:3]].append(v)
    bases={k:independent_median(v) for k,v in bases.items()}
    compounds=defaultdict(list)
    for k,v in bases.items():compounds[k[:2]].append(v)
    return contexts,bases,{k:independent_median(v) for k,v in compounds.items()}

def vector_scores(X,queries):
    names=list(queries);out=np.full((len(X),len(names)),np.nan,dtype=np.float64)
    for qi,q in enumerate(names):
        w=queries[q];ids=np.flatnonzero(w);v=w[ids];wn=np.sqrt(np.sum(v*v,dtype=np.float64))
        for lo in range(0,len(X),256):
            z=np.asarray(X[lo:lo+256])[:,ids]
            zn=np.sqrt(np.sum(z*z,axis=1,dtype=np.float64))
            with np.errstate(invalid='ignore',divide='ignore'):out[lo:lo+len(z),qi]=(z@v)/(zn*wn)
        print('PROFILE_QUERY_COMPLETE',q,flush=True)
    return out
