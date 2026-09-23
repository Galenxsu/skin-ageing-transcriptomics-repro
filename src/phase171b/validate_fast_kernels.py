from null171b import *
from fast_null_kernels import sparse_scores,parallel_aggregate

rows,genes,pids,queries,strata=metadata();h=hierarchy(rows,pids)
X=np.load(O/'MATRIX/frozen_107201x8227_float64.npy',mmap_mode='r')
obs=np.load(O/'RESULTS/observed_compound_scores.npy');names=jr(O/'RESULTS/query_order.json')
reports=[]
for name,qs,B,seed in plans():
    rng=np.random.Generator(np.random.PCG64(seed));Ws,hashes=permutations(rng,queries,strata,qs,32)
    for q in qs:
        with threadpool_limits(limits=1):
            t=time.perf_counter();ref=null_profile_scores(X,Ws[q],name=='CONTINUOUS');oldtime=time.perf_counter()-t
            t=time.perf_counter();val=ref.copy() if name=='CONTINUOUS' else sparse_scores(X,Ws[q]);newtime=time.perf_counter()-t
        t=time.perf_counter();rc,rb,rp=aggregate(ref,h);serialtime=time.perf_counter()-t
        t=time.perf_counter();vc,vb,vp=parallel_aggregate(val,h);parallel_time=time.perf_counter()-t
        for level,a,b in [('PROFILE',val,ref),('CONTEXT',vc,rc),('BASE',vb,rb),('COMPOUND',vp,rp)]:
            finite=np.isfinite(a)&np.isfinite(b);d=np.abs(a-b);ok=np.all(d[finite]<=1e-12+1e-10*np.abs(b[finite]))
            stateok=np.array_equal(np.isnan(a),np.isnan(b)) and np.array_equal(np.isposinf(a),np.isposinf(b)) and np.array_equal(np.isneginf(a),np.isneginf(b))
            signok=np.array_equal(np.sign(a[finite]),np.sign(b[finite]))
            maxij=np.unravel_index(np.nanargmax(d),d.shape)
            detail=compare(float(a[maxij]),float(b[maxij]),'direction')
            reports.append(dict(query=q,level=level,finite_comparisons=int(finite.sum()),all_tolerance_pass=bool(ok),all_state_pass=bool(stateok),all_direction_pass=bool(signok),reference_seconds=oldtime,accelerated_seconds=newtime,serial_aggregation_seconds=serialtime,parallel_aggregation_seconds=parallel_time,max_difference_index=str(maxij),**detail))
            assert ok and stateok and signok,'FLOAT_COMPARISON_FAILURE'
        o=obs[:,names.index(q)][:,None]
        assert np.array_equal(np.sum(vp<=o,axis=1),np.sum(rp<=o,axis=1)),'LOWER_EXTREME_COUNT_DISCORDANCE'
        assert np.array_equal(np.sum(vp>=o,axis=1),np.sum(rp>=o,axis=1)),'UPPER_EXTREME_COUNT_DISCORDANCE'
        print('ACCELERATION_VALIDATED',q,'reference',round(oldtime,3),'candidate',round(newtime,3),'aggregate',round(serialtime,3),round(parallel_time,3),flush=True)
csvwrite(O/'QA/ACCELERATED_IMPLEMENTATION_COMPARISON.csv',reports)
js(O/'QA/ACCELERATED_IMPLEMENTATION_VALIDATION.json',dict(utc=now(),pass_all=True,queries=len(queries),iterations_per_query=32,scientific_rule_changes='NONE',all_null_extreme_counts_exact=True,kernel_sha256=sha(O/'SCRIPTS/fast_null_kernels.py'),comparison_csv_sha256=sha(O/'QA/ACCELERATED_IMPLEMENTATION_COMPARISON.csv')))
