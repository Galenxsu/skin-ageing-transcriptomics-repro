from engine171b import *
from numeric_rules import compare
from threadpoolctl import threadpool_limits
from fast_null_kernels import continuous_scores,sparse_scores,parallel_aggregate
import time, argparse
ENGINE_HASH_AT_IMPORT=sha(Path(__file__))
INPUT_MANIFEST_HASH=sha(O/'EXECUTION_INPUT_HASH_BEFORE.csv')
PARAMETER_HASH=hashlib.sha256((sha(A/'FREEZE_MANIFEST.csv')+'|'+sha(O/'ADDENDUM_SHA256_MANIFEST.csv')).encode()).hexdigest()

def plans():
    return [('CONTINUOUS', ['G2_R2_CAP8','G2_R1'],99999,171260901)]+[(q,[q],299999,171261000+i) for i,q in enumerate(['G2_E25','G2_E50','G2_E75','G2_E100','G2_E150','G2_E200','G2_E300'])]+[(q,[q],99999,171262000+i) for i,q in enumerate(['HALLMARK_MTORC1_SIGNALING','KEGG_LYSOSOME'])]+[(q,[q],199999,171263000+i) for i,q in enumerate(['GOBP_EPIDERMIS_DEVELOPMENT','REACTOME_KERATINIZATION','REACTOME_AUTOPHAGY'])]

def permutations(rng,queries,strata,qs,n):
    mats={q:np.zeros((8227,n),dtype=np.float64) for q in qs};hashes=[]
    for j in range(n):
        digest=hashlib.sha256()
        for idx in strata:
            perm=rng.permutation(idx)
            digest.update(perm.astype('<i8').tobytes())
            for q in qs:
                original=queries[q]
                if q.startswith('G2_E'):
                    up=int(np.count_nonzero(original[idx]>0));down=int(np.count_nonzero(original[idx]<0))
                    k=int(q.split('E')[1]);mats[q][perm[:up],j]=.5/k;mats[q][perm[up:up+down],j]=-.5/k
                else:mats[q][idx,j]=original[perm]
        hashes.append(digest.hexdigest())
    return mats,hashes

def null_profile_scores(X,W,continuous):
    n=W.shape[1];out=np.full((len(X),n),np.nan,dtype=np.float64)
    wn=np.sqrt(np.sum(W*W,axis=0,dtype=np.float64))
    if continuous:
        for lo in range(0,len(X),256):
            z=np.asarray(X[lo:lo+256]);zn=np.sqrt(np.sum(z*z,axis=1,dtype=np.float64))
            out[lo:lo+len(z)]=(z@W)/(zn[:,None]*wn[None,:])
    else:
        # Sparse members only; no full-expression reweighting or imputation.
        for j in range(n):
            ids=np.flatnonzero(W[:,j]);w=W[ids,j]
            for lo in range(0,len(X),256):
                z=np.asarray(X[lo:lo+256])[:,ids]
                zn=np.sqrt(np.sum(z*z,axis=1,dtype=np.float64))
                out[lo:lo+len(z),j]=(z@w)/(zn*wn[j])
    return out

def independent_null_audit(X,W,scores,rows,h,qscores,q,iteration_start,ix,complete):
    tests=[];refscores=np.empty((len(ix),W.shape[1]),dtype=np.float64)
    compounds={}
    for j in range(W.shape[1]):
        for k,i in enumerate(ix):
            v=scalar_score(X[i],W[:,j]);refscores[k,j]=v
            tests.append(dict(query=q,iteration=iteration_start+j,level='PROFILE',key=rows[i]['sig_id'],**compare(float(scores[i,j]),v,'direction')))
        _,_,ip=independent_aggregate([(rows[i],refscores[k,j]) for k,i in enumerate(ix) if rows[i]['pert_id'] in complete])
        for k,key in enumerate(h['compound_keys']):
            if key in ip:
                tests.append(dict(query=q,iteration=iteration_start+j,level='COMPOUND',key='|'.join(key),**compare(float(qscores[k,j]),ip[key],'direction')))
                compounds.setdefault(key,[]).append(ip[key])
    return tests,compounds

def run(pilot_only=False):
    print('NULL_PREPARATION_START',flush=True)
    assert jr(O/'QA/INDEPENDENT_OBSERVED_SUMMARY.json')['failed']==0
    assert jr(O/'QA/ACCELERATED_IMPLEMENTATION_VALIDATION.json')['pass_all']
    assert jr(O/'QA/CONTINUOUS_THREAD_VALIDATION.json')['pass_all']
    rows,genes,pids,queries,strata=metadata();h=hierarchy(rows,pids)
    print('NULL_METADATA_LINKED',len(rows),len(pids),flush=True)
    X=np.load(O/'MATRIX/frozen_107201x8227_float64.npy',mmap_mode='r')
    fixed={r['sig_id'] for r in csvread(O/'QA/FROZEN_INDEPENDENT_AUDIT_SAMPLE.csv')}
    ix=[i for i,r in enumerate(rows) if r['sig_id'] in fixed]
    bypert=defaultdict(list)
    byfold=defaultdict(list)
    for r in rows:bypert[r['pert_id']].append(r);byfold[(r['pert_id'],r['fold'])].append(r)
    complete={p for p,rr in bypert.items() if all(r['sig_id'] in fixed for r in rr)}
    obs=np.load(O/'RESULTS/observed_compound_scores.npy');names=jr(O/'RESULTS/query_order.json')
    context_counts=defaultdict(int)
    for k in h['context_keys']:context_counts[k[:2]]+=1
    support=np.array([len(byfold[(p,f)])>=3 and context_counts[(p,f)]>=2 for p,f in h['compound_keys']])
    npert=len(pids)
    event('NULL_ENGINE_START',pilot_only=pilot_only,engine_sha256=ENGINE_HASH_AT_IMPORT,numpy_version=np.__version__)
    for name,qs,B,seed in plans():
        checkpoint=O/'CHECKPOINTS'/f'{name}.json'
        rng=np.random.Generator(np.random.PCG64(seed));start=0
        paths={(q,f):O/'NULL'/f'{q}_{f}.npy' for q in qs for f in (['DISCOVERY','VALIDATION'] if q=='G2_R2_CAP8' else ['DISCOVERY'])}
        if checkpoint.exists():
            cp=jr(checkpoint);assert cp['seed']==seed and cp['B']==B
            if cp['engine_sha256']!=ENGINE_HASH_AT_IMPORT:
                transition=jr(O/'QA/NULL_IMPLEMENTATION_TRANSITION_VERIFICATION.json')
                assert transition['new_engine_sha256']==ENGINE_HASH_AT_IMPORT and transition['all_existing_chunks_reproduced_exactly']
                assert cp['engine_sha256'] in transition['previous_checkpoint_engine_sha256']
            rng.bit_generator.state=cp['rng_state'];start=cp['next_iteration']
            maps={k:np.load(p,mmap_mode='r+') for k,p in paths.items()}
            for (q,f),mm in maps.items():
                saved=np.asarray(mm[cp['chunk_start']:cp['chunk_end']])
                assert hashlib.sha256(saved.astype('<f8').tobytes()).hexdigest()==cp['chunk_hashes'][q+'_'+f],'Checkpoint output hash mismatch'
            if start==B:continue
            if pilot_only and start>=32:continue
        else:
            assert not any(p.exists() for p in paths.values()),'Uncheckpointed null file: stop for recovery audit'
            maps={k:np.lib.format.open_memmap(p,mode='w+',dtype=np.float64,shape=(B,npert)) for k,p in paths.items()}
        stop=min(B,32) if pilot_only else B
        for lo in range(start,stop,32):
            tic=time.perf_counter();n=min(32,stop-lo)
            W,hashes=permutations(rng,queries,strata,qs,n)
            chunkoutputs={};audits=[]
            with threadpool_limits(limits=1):
                for q in qs:
                    scores=continuous_scores(X,W[q]) if q in ['G2_R2_CAP8','G2_R1'] else sparse_scores(X,W[q])
                    _,_,comp=parallel_aggregate(scores,h)
                    if lo==0:
                        audit,reference=independent_null_audit(X,W[q][:,:2],scores[:,:2],rows,h,comp[:,:2],q,lo,ix,complete)
                        audits.extend(audit)
                        # Extreme counts in the two-iteration audit are exact, not tolerance based.
                        for k,key in enumerate(h['compound_keys']):
                            if key in reference:
                                a=np.asarray(reference[key]);v=float(obs[k,names.index(q)])
                                # Observed reference recomputed independently from original expression.
                                pairs=[(r,scalar_score(X[i],queries[q])) for i,r in enumerate(rows) if (r['pert_id'],r['fold'])==key]
                                _,_,ip=independent_aggregate(pairs);vref=ip[key]
                                lower1=int(np.count_nonzero(comp[k,:2]<=v));lower2=int(np.count_nonzero(a<=vref))
                                upper1=int(np.count_nonzero(comp[k,:2]>=v));upper2=int(np.count_nonzero(a>=vref))
                                if (lower1,upper1)!=(lower2,upper2):raise RuntimeError('EXACT_NULL_EXTREME_COUNT_DISCORDANCE')
                    comp[~support,:]=np.nan
                    for (qq,fold),mm in maps.items():
                        if qq!=q:continue
                        sl=slice(0,npert) if fold=='DISCOVERY' else slice(npert,2*npert)
                        values=comp[sl].T.copy();mm[lo:lo+n]=values
                        chunkoutputs[q+'_'+fold]=hashlib.sha256(values.astype('<f8').tobytes()).hexdigest()
            if audits:
                csvwrite(O/'QA'/f'INDEPENDENT_NULL_PILOT_{name}.csv',audits)
                if any(not(t['numeric_pass'] and t['decision_pass']) for t in audits):
                    event('INDEPENDENT_NULL_PILOT_FAILED',query=name)
                    raise RuntimeError('INDEPENDENT_NULL_PILOT_DISCORDANCE')
            for mm in maps.values():mm.flush()
            with open(O/'NULL'/f'{name}_iteration_member_hashes.csv','a',encoding='utf-8',newline='') as f:
                wr=csv.writer(f)
                if lo==0:wr.writerow(['iteration_zero_based','permutation_member_sha256'])
                wr.writerows((lo+j,d) for j,d in enumerate(hashes))
            assert sha(Path(__file__))==ENGINE_HASH_AT_IMPORT,'Running implementation file changed; do not commit checkpoint'
            cp=dict(name=name,seed=seed,B=B,next_iteration=lo+n,rng_state=rng.bit_generator.state,engine_sha256=ENGINE_HASH_AT_IMPORT,input_manifest_sha256=INPUT_MANIFEST_HASH,parameter_sha256=PARAMETER_HASH,chunk_start=lo,chunk_end=lo+n,chunk_hashes=chunkoutputs,utc=now(),elapsed_seconds=time.perf_counter()-tic)
            js(checkpoint,cp)
            js(O/'CURRENT_STATUS.json',dict(status='RUNNING_FROZEN_NULL_NO_FINAL_RESULTS',**cp))
            with open(O/'NULL/null_execution_log.jsonl','a',encoding='utf-8') as f:f.write(json.dumps(cp)+'\n')
            if lo==0 or (lo//32)%25==0: print('NULL_CHECKPOINT',name,lo+n,'/',B,'seconds',round(cp['elapsed_seconds'],3),flush=True)
        if pilot_only:break
    event('NULL_ENGINE_EXIT',pilot_only=pilot_only)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--pilot-only',action='store_true');args=parser.parse_args();run(args.pilot_only)
