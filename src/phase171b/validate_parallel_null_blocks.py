from null171b import *
from concurrent.futures import ThreadPoolExecutor
from numba import set_num_threads,threading_layer

rows,genes,pids,queries,strata=metadata();h=hierarchy(rows,pids);X=np.load(O/'MATRIX/frozen_107201x8227_float64.npy',mmap_mode='r')
def calculate(q,W):
    set_num_threads(8)
    s=continuous_scores(X,W) if q.startswith('G2_R') else sparse_scores(X,W)
    _,_,c=parallel_aggregate(s,h)
    return s,c

# Compile before starting concurrent callers; no alternative statistical rule.
parallel_aggregate(np.zeros((len(rows),1),dtype=np.float64),h)
layer=threading_layer();assert layer!='workqueue','Concurrent scheduling not supported by this runtime; retain sequential driver'
checks=[]
with threadpool_limits(limits=1):
    for q,seed in [('G2_R2_CAP8',171260901),('G2_E25',171261000),('G2_E300',171261006),('REACTOME_KERATINIZATION',171263001)]:
        rng=np.random.Generator(np.random.PCG64(seed));jobs=[permutations(rng,queries,strata,[q],32)[0][q] for _ in range(3)]
        tic=time.perf_counter()
        with ThreadPoolExecutor(max_workers=3) as pool:parallel=list(pool.map(lambda W:calculate(q,W),jobs))
        parallel_seconds=time.perf_counter()-tic;tic=time.perf_counter()
        for j,W in enumerate(jobs):
            profile,compound=calculate(q,W);exact_p=np.array_equal(profile,parallel[j][0],equal_nan=True);exact_c=np.array_equal(compound,parallel[j][1],equal_nan=True)
            assert exact_p and exact_c,'Concurrent scheduling changed arithmetic'
        sequential_seconds=time.perf_counter()-tic
        checks.append(dict(query=q,blocks=3,iterations_per_block=32,profile_values_per_block=107201*32,exact_profile_and_compound_values=True,parallel_seconds=parallel_seconds,sequential_seconds=sequential_seconds))
        print(checks[-1],flush=True)
js(O/'QA/PARALLEL_NULL_BLOCK_VALIDATION.json',dict(utc=now(),pass_all=True,threading_layer=layer,workers=3,Numba_threads_per_worker=8,BLAS_threads=1,checks=checks,scientific_rule_changes='NONE',execution='Random members generated serially in frozen iteration/stratum order; score blocks concurrently; commit checkpoints in original iteration order.'))
