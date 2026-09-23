from null171b import *
from concurrent.futures import ThreadPoolExecutor
from collections import deque
from numba import set_num_threads
import traceback

DRIVER_HASH=sha(Path(__file__))

def run_parallel():
    handoff=jr(O/'QA/SCHEDULER_HANDOFF_VERIFICATION.json')
    assert handoff['all_pass'] and handoff['new_driver_sha256']==DRIVER_HASH
    assert jr(O/'QA/PARALLEL_NULL_BLOCK_VALIDATION.json')['pass_all']
    assert jr(O/'QA/FULL_INDEPENDENT_OBSERVED_VALIDATION.json')['pass_all']
    assert all(int(r['zero_norm_profile_iterations'])==0 for r in csvread(O/'QA/ALL_FROZEN_NULL_ZERO_NORM_AUDIT.csv'))
    for r in csvread(A/'FREEZE_MANIFEST.csv'):assert sha(A/r['relative_path'])==r['sha256']
    rows,genes,pids,queries,strata=metadata();h=hierarchy(rows,pids)
    X=np.load(O/'MATRIX/frozen_107201x8227_float64.npy',mmap_mode='r');names=jr(O/'RESULTS/query_order.json');obs=np.load(O/'RESULTS/observed_compound_scores.npy')
    bypert=defaultdict(list);byfold=defaultdict(list);cx=defaultdict(int)
    for i,r in enumerate(rows):bypert[r['pert_id']].append(i);byfold[r['pert_id'],r['fold']].append(i)
    for key in h['context_keys']:cx[key[:2]]+=1
    support=np.array([len(byfold[k])>=3 and cx[k]>=2 for k in h['compound_keys']]);npert=len(pids)
    fixed={r['sig_id'] for r in csvread(O/'QA/FROZEN_INDEPENDENT_AUDIT_SAMPLE.csv')};ix=[i for i,r in enumerate(rows) if r['sig_id'] in fixed]
    complete={p for p,ii in bypert.items() if all(rows[i]['sig_id'] in fixed for i in ii)}
    # Warm compiled aggregation before concurrent calls. No random generator use.
    parallel_aggregate(np.zeros((len(rows),1),dtype=np.float64),h)
    def calculate(job):
        set_num_threads(8);lo,W,hashes,state,qs=job;out={};audits=[];tic=time.perf_counter()
        for q in qs:
            s=continuous_scores(X,W[q]) if q in ['G2_R2_CAP8','G2_R1'] else sparse_scores(X,W[q])
            assert np.all(np.isfinite(s)),'Unexpected nonfinite null profile: no dropping or imputation permitted'
            _,_,c=parallel_aggregate(s,h)
            if lo==0:
                tests,refs=independent_null_audit(X,W[q][:,:2],s[:,:2],rows,h,c[:,:2],q,lo,ix,complete);audits.extend(tests)
                for k,key in enumerate(h['compound_keys']):
                    if key in refs:
                        _,_,ind=independent_aggregate([(rows[i],scalar_score(X[i],queries[q])) for i in byfold[key]])
                        v=obs[k,names.index(q)];vr=ind[key];a=np.asarray(refs[key])
                        assert (int(np.sum(c[k,:2]<=v)),int(np.sum(c[k,:2]>=v)))==(int(np.sum(a<=vr)),int(np.sum(a>=vr))),'EXACT_NULL_EXTREME_COUNT_DISCORDANCE'
            c[~support]=np.nan;out[q]=c
        assert all(t['numeric_pass'] and t['decision_pass'] for t in audits),'Independent null pilot disagreement'
        return lo,hashes,state,out,audits,time.perf_counter()-tic
    event('PARALLEL_NULL_CONTINUATION_START',pid=os.getpid(),driver_sha256=DRIVER_HASH,scoring_kernel_sha256=ENGINE_HASH_AT_IMPORT,workers=3,scientific_rule_changes='NONE')
    for name,qs,B,seed in plans():
        cpfile=O/'CHECKPOINTS'/f'{name}.json';rng=np.random.Generator(np.random.PCG64(seed));start=0
        paths={(q,f):O/'NULL'/f'{q}_{f}.npy' for q in qs for f in (['DISCOVERY','VALIDATION'] if q=='G2_R2_CAP8' else ['DISCOVERY'])}
        if cpfile.exists():
            cp=jr(cpfile);assert cp['seed']==seed and cp['B']==B
            assert cp['engine_sha256'] in (DRIVER_HASH,ENGINE_HASH_AT_IMPORT,'3bf81fc57933bf83f42fd05abf297adddcda0d6c3c75fe945fbc5cfa9d5baec3')
            start=cp['next_iteration'];rng.bit_generator.state=cp['rng_state']
            maps={k:np.load(p,mmap_mode='r+') for k,p in paths.items()}
            for (q,f),mm in maps.items():
                data=np.asarray(mm[cp['chunk_start']:cp['chunk_end']]);assert hashlib.sha256(data.astype('<f8').tobytes()).hexdigest()==cp['chunk_hashes'][q+'_'+f]
            if start==B:
                for mm in maps.values():mm._mmap.close()
                continue
        else:
            assert not any(p.exists() for p in paths.values()),'Uncheckpointed files require recovery audit'
            maps={k:np.lib.format.open_memmap(p,mode='w+',dtype=np.float64,shape=(B,npert)) for k,p in paths.items()}
        def new_job(lo):
            W,hashes=permutations(rng,queries,strata,qs,min(32,B-lo));return lo,W,hashes,json.loads(json.dumps(rng.bit_generator.state)),qs
        with threadpool_limits(limits=1),ThreadPoolExecutor(max_workers=3) as pool:
            pending=deque();next_submit=start
            while next_submit<B and len(pending)<3:
                pending.append(pool.submit(calculate,new_job(next_submit)));next_submit+=min(32,B-next_submit)
            while pending:
                lo,hashes,state,out,audits,seconds=pending.popleft().result();n=len(hashes);chunkhash={}
                for (q,fold),mm in maps.items():
                    sl=slice(0,npert) if fold=='DISCOVERY' else slice(npert,2*npert);values=out[q][sl].T.copy();mm[lo:lo+n]=values;mm.flush();chunkhash[q+'_'+fold]=hashlib.sha256(values.astype('<f8').tobytes()).hexdigest()
                if audits:csvwrite(O/'QA'/f'INDEPENDENT_NULL_PILOT_{name}.csv',audits)
                assert sha(Path(__file__))==DRIVER_HASH and sha(O/'SCRIPTS/null171b.py')==ENGINE_HASH_AT_IMPORT
                with open(O/'NULL'/f'{name}_iteration_member_hashes.csv','a',encoding='utf-8',newline='') as f:
                    wr=csv.writer(f)
                    if lo==0:wr.writerow(['iteration_zero_based','permutation_member_sha256'])
                    wr.writerows((lo+j,d) for j,d in enumerate(hashes))
                cp=dict(name=name,seed=seed,B=B,next_iteration=lo+n,rng_state=state,engine_sha256=DRIVER_HASH,scoring_kernel_sha256=ENGINE_HASH_AT_IMPORT,input_manifest_sha256=INPUT_MANIFEST_HASH,parameter_sha256=PARAMETER_HASH,chunk_start=lo,chunk_end=lo+n,chunk_hashes=chunkhash,utc=now(),elapsed_seconds=seconds,workers=3,commit_order='FROZEN_ITERATION_ORDER')
                js(cpfile,cp)
                with open(O/'NULL/null_execution_log.jsonl','a',encoding='utf-8') as f:f.write(json.dumps(cp)+'\n')
                js(O/'CURRENT_STATUS.json',dict(status='RUNNING_FROZEN_NULL_NO_FINAL_RESULTS',**cp))
                if (lo//32)%25==0:print('NULL_PARALLEL_CHECKPOINT',name,lo+n,'/',B,flush=True)
                if next_submit<B:
                    pending.append(pool.submit(calculate,new_job(next_submit)));next_submit+=min(32,B-next_submit)
        for mm in maps.values():mm._mmap.close()
    js(O/'CURRENT_STATUS.json',dict(status='NULL_COMPLETED_INDEPENDENT_FULL_COUNT_AND_FINAL_AUDITS_PENDING',utc=now()))
    event('ALL_FROZEN_NULL_COMPLETED',driver_sha256=DRIVER_HASH)

if __name__=='__main__':
    js(O/'NULL_RUN_PROCESS.json',dict(pid=os.getpid(),utc=now(),engine_sha256=DRIVER_HASH))
    try:run_parallel()
    except BaseException as e:
        report=dict(status='STOPPED_NULL_EXECUTION_ERROR_NO_FINAL_CLASSIFICATION',utc=now(),error_type=type(e).__name__,error=str(e),traceback=traceback.format_exc())
        js(O/'NULL_EXECUTION_ERROR.json',report);js(O/'CURRENT_STATUS.json',report);event('PARALLEL_NULL_STOPPED',**report);raise
