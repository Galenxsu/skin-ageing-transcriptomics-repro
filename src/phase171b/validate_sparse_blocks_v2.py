from null171b import *
from fast_sparse_blocks_v2 import sparse_scores_v2

rows,genes,pids,queries,strata=metadata();X=np.load(O/'MATRIX/frozen_107201x8227_float64.npy',mmap_mode='r');checks=[]
for name,qs,B,seed in plans():
    if name=='CONTINUOUS':continue
    q=qs[0];rng=np.random.Generator(np.random.PCG64(seed));W,_=permutations(rng,queries,strata,qs,32)
    tic=time.perf_counter();a=sparse_scores_v2(X,W[q]);new_seconds=time.perf_counter()-tic
    tic=time.perf_counter();b=sparse_scores(X,W[q]);old_seconds=time.perf_counter()-tic
    exact=bool(np.array_equal(a,b,equal_nan=True));assert exact,'Scheduling-only optimization did not preserve exact arithmetic'
    checks.append(dict(query=q,signature_block=256,null_block=32,exact_all_profile_null_values=exact,values_compared=int(a.size),v2_seconds=new_seconds,v1_seconds=old_seconds))
    print(checks[-1],flush=True)
js(O/'QA/SPARSE_BLOCK_SCHEDULING_V2_VALIDATION.json',dict(utc=now(),pass_all=True,kernel_sha256=sha(O/'SCRIPTS/fast_sparse_blocks_v2.py'),checks=checks,scientific_rule_changes='NONE',change='Schedule existing 256-signature blocks in one compiled parallel loop; per-profile arithmetic and null block32 unchanged.'))
