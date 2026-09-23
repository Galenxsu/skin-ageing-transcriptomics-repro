from null171b import *

rows,genes,pids,queries,strata=metadata();h=hierarchy(rows,pids)
X=np.load(O/'MATRIX/frozen_107201x8227_float64.npy',mmap_mode='r')
byfold=defaultdict(list);contexts=defaultdict(int)
for r in rows:byfold[(r['pert_id'],r['fold'])].append(r)
for k in h['context_keys']:contexts[k[:2]]+=1
support=np.array([len(byfold[k])>=3 and contexts[k]>=2 for k in h['compound_keys']])
checks=[];old=set();times=[]
for name,qs,B,seed in plans():
    p=O/'CHECKPOINTS'/f'{name}.json'
    if not p.exists():continue
    cp=jr(p);assert cp['next_iteration']==32
    old.add(cp['engine_sha256']);rng=np.random.Generator(np.random.PCG64(seed))
    W,hashes=permutations(rng,queries,strata,qs,32)
    mh=csvread(O/'NULL'/f'{name}_iteration_member_hashes.csv')
    assert [r['permutation_member_sha256'] for r in mh]==hashes
    assert rng.bit_generator.state==cp['rng_state']
    for q in qs:
        tic=time.perf_counter()
        with threadpool_limits(limits=8):
            s=null_profile_scores(X,W[q],q in ['G2_R2_CAP8','G2_R1'])
            _,_,c=aggregate(s,h)
        c[~support]=np.nan
        times.append(dict(query=q,iterations=32,seconds=time.perf_counter()-tic))
        for fold in (['DISCOVERY','VALIDATION'] if q=='G2_R2_CAP8' else ['DISCOVERY']):
            sl=slice(0,len(pids)) if fold=='DISCOVERY' else slice(len(pids),2*len(pids))
            expected=c[sl].T.copy();actual=np.load(O/'NULL'/f'{q}_{fold}.npy',mmap_mode='r')[:32]
            digest=hashlib.sha256(expected.astype('<f8').tobytes()).hexdigest()
            checks.append(dict(query=q,fold=fold,exact_array_match=bool(np.array_equal(actual,expected,equal_nan=True)),output_hash_match=digest==cp['chunk_hashes'][q+'_'+fold],iteration_hashes_match=True,rng_state_match=True))
    print('REPLAY_VERIFIED',name,flush=True)
result=dict(utc=now(),new_engine_sha256=ENGINE_HASH_AT_IMPORT,previous_checkpoint_engine_sha256=sorted(old),all_existing_chunks_reproduced_exactly=all(r['exact_array_match'] and r['output_hash_match'] for r in checks),checks=checks,timing=times,incident='Tool session reported -1, but first null process continued; metadata setup was optimized before process completion. The old process recorded file hash after source edit. Its start-event hash identifies loaded source. Both existing 32-iteration chunks were replayed from frozen seed; exact values, hashes, permutation hashes and RNG state verified. Original logs retained.',scientific_rule_changes='NONE',repair='Capture engine hash once at import; verify unchanged before checkpoint commit; metadata lookup optimization does not change members/order/math/random process.')
js(O/'QA/NULL_IMPLEMENTATION_TRANSITION_VERIFICATION.json',result)
assert result['all_existing_chunks_reproduced_exactly']
event('NULL_CHECKPOINT_PROVENANCE_REPLAY_VERIFIED',report_sha256=sha(O/'QA/NULL_IMPLEMENTATION_TRANSITION_VERIFICATION.json'))
print(json.dumps(result),flush=True)
