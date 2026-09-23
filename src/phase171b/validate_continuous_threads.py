from null171b import *
from fast_null_kernels import continuous_scores,parallel_aggregate

rows,genes,pids,queries,strata=metadata();h=hierarchy(rows,pids)
X=np.load(O/'MATRIX/frozen_107201x8227_float64.npy',mmap_mode='r')
rng=np.random.Generator(np.random.PCG64(171260901));Ws,_=permutations(rng,queries,strata,['G2_R2_CAP8','G2_R1'],32)
obs=np.load(O/'RESULTS/observed_compound_scores.npy');names=jr(O/'RESULTS/query_order.json');checks=[]
for q,W in Ws.items():
    with threadpool_limits(limits=1):
        tic=time.perf_counter();a=continuous_scores(X,W);elapsed=time.perf_counter()-tic
        b=null_profile_scores(X,W,True)
    exact=bool(np.array_equal(a,b,equal_nan=True));diff=float(np.nanmax(np.abs(a-b)))
    _,_,ca=parallel_aggregate(a,h);_,_,cb=aggregate(b,h)
    oo=obs[:,names.index(q)][:,None]
    counts=bool(np.array_equal(np.sum(ca<=oo,axis=1),np.sum(cb<=oo,axis=1)) and np.array_equal(np.sum(ca>=oo,axis=1),np.sum(cb>=oo,axis=1)))
    assert exact and counts
    checks.append(dict(query=q,exact_profile_match=exact,exact_extreme_counts=counts,max_abs_diff=diff,parallel_seconds=elapsed))
    print(checks[-1],flush=True)
js(O/'QA/CONTINUOUS_THREAD_VALIDATION.json',dict(utc=now(),checks=checks,pass_all=True,kernel_sha256=sha(O/'SCRIPTS/fast_null_kernels.py'),scientific_rule_changes='NONE'))
