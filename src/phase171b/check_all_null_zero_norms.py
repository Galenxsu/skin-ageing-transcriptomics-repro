from engine171b import *
from null171b import plans
import time

@njit(cache=False)
def zero_rows(z,ids):
    count=0
    for i in range(z.shape[0]):
        allzero=True
        for j in ids:
            if not z[i,j]:allzero=False;break
        if allzero:count+=1
    return count

rows,genes,pids,queries,strata=metadata()
zs=np.empty((107201,len(strata)),dtype=np.int64)
for lo in range(0,107201,256):
    x=np.load(O/'MATRIX/frozen_107201x8227_float64.npy',mmap_mode='r');z=x[lo:lo+256]==0
    for j,idx in enumerate(strata):zs[lo:lo+len(z),j]=z[:,idx].sum(axis=1)
    x._mmap.close()
records=[]
for name,qs,B,seed in plans():
    if name=='CONTINUOUS':continue
    q=qs[0];w=queries[q];needed=np.array([np.count_nonzero(w[idx]) for idx in strata]);risk=np.flatnonzero(np.all(zs>=needed,axis=1))
    if len(risk)==0:
        records.append(dict(query=q,B=B,seed=seed,at_risk_profiles=0,iterations_checked=B,zero_norm_profile_iterations=0,proof='STRATUM_CARDINALITY_EXCLUDES_ZERO_NORM',seconds=0));continue
    x=np.load(O/'MATRIX/frozen_107201x8227_float64.npy',mmap_mode='r');mask=x[risk]==0;x._mmap.close()
    rng=np.random.Generator(np.random.PCG64(seed));fail=0;tic=time.perf_counter()
    for iteration in range(B):
        active=[]
        for idx in strata:
            perm=rng.permutation(idx)
            if q.startswith('G2_E'):
                count=np.count_nonzero(w[idx]);active.extend(perm[:count])
            else:active.extend(idx[w[perm]!=0])
        fail+=int(zero_rows(mask,np.asarray(active,dtype=np.int64)))
        if iteration%25000==0:print('NULL_ZERO_NORM_AUDIT',q,iteration,'/',B,'failures',fail,flush=True)
    records.append(dict(query=q,B=B,seed=seed,at_risk_profiles=len(risk),iterations_checked=B,zero_norm_profile_iterations=fail,proof='ALL_FROZEN_RANDOM_MEMBERS_CHECKED_NO_SCORES_COMPUTED',seconds=time.perf_counter()-tic))
    csvwrite(O/'QA/ALL_FROZEN_NULL_ZERO_NORM_AUDIT.csv',records)
csvwrite(O/'QA/ALL_FROZEN_NULL_ZERO_NORM_AUDIT.csv',records)
assert all(r['zero_norm_profile_iterations']==0 for r in records),'Frozen null includes zero-norm profiles: must stop and preserve affected query; no masked dropping permitted'
event('ALL_FROZEN_NULL_PROFILE_NORMS_PROVEN_NONZERO',queries=len(records),same_evaluable_mask=True)
