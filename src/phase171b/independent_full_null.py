from common171b import *
import numpy as np, math
from collections import defaultdict
sys.path.insert(0,str(O/'runtime_dependencies'))
from numba import njit
from numeric_rules import compare

@njit(cache=False,fastmath=False)
def reference_group(values,indices,boundaries):
    out=np.full((len(boundaries)-1,values.shape[1]),np.nan)
    for g in range(len(boundaries)-1):
        n=boundaries[g+1]-boundaries[g]
        for j in range(values.shape[1]):
            ordered=np.empty(n,dtype=np.float64);count=0
            for a in range(boundaries[g],boundaries[g+1]):
                value=values[indices[a],j]
                if not math.isfinite(value):continue
                k=count
                while k>0 and ordered[k-1]>value:
                    ordered[k]=ordered[k-1];k-=1
                ordered[k]=value;count+=1
            if count:
                mid=count//2
                out[g,j]=ordered[mid] if count%2 else (ordered[mid-1]+ordered[mid])*0.5
    return out

def independent_groups(keys):
    mapping=defaultdict(list)
    for i,k in enumerate(keys):mapping[k].append(i)
    labels=sorted(mapping);indices=[];boundaries=[0]
    for k in labels:indices.extend(mapping[k]);boundaries.append(len(indices))
    return labels,(np.asarray(indices,dtype=np.int64),np.asarray(boundaries,dtype=np.int64))

def run():
    requested={r['pert_id'] for r in csvread(O/'QA/CRITICAL_COMPOUND_INDEPENDENT_NULL_AUDIT_SET.csv')}
    allrows=[r for r in csvread(A/'discovery_validation_assignment_preview.csv') if r['metadata_eligible']=='True']
    selected=[i for i,r in enumerate(allrows) if r['pert_id'] in requested];rows=[allrows[i] for i in selected]
    genes=csvread(A/'continuous_signature_definition.csv');w=np.asarray([float(r['main_R2_L2_weight']) for r in genes],dtype=np.float64)
    X=np.empty((len(selected),8227),dtype=np.float64)
    for lo in range(0,len(selected),256):
        source=np.load(O/'MATRIX/frozen_107201x8227_float64.npy',mmap_mode='r');X[lo:lo+256]=source[selected[lo:lo+256]];source._mmap.close()
    ck,pc=independent_groups([(r['pert_id'],r['fold'],r['base_cell_id'],r['cell_id'],float(r['dose_uM']),float(r['time_h'])) for r in rows])
    bk,cb=independent_groups([k[:3] for k in ck]);pk,bp=independent_groups([k[:2] for k in bk])
    def chain(s):return reference_group(reference_group(reference_group(s,*pc),*cb),*bp)
    counts=defaultdict(int);ccounts=defaultdict(int)
    for r in rows:counts[r['pert_id'],r['fold']]+=1
    for k in ck:ccounts[k[:2]]+=1
    support=np.array([counts[k]>=3 and ccounts[k]>=2 for k in pk])
    zn=np.sqrt(np.einsum('ij,ij->i',X,X,optimize=False))
    wn=math.sqrt(float(np.einsum('i,i->',w,w,optimize=False)))
    observed=chain((np.einsum('ij,j->i',X,w,optimize=False)/(zn*wn))[:,None])[:,0]
    stats={(r['pert_id'],r['fold']):r for r in csvread(O/'RESULTS/all_query_statistics.csv') if r['query_id']=='G2_R2_CAP8'}
    tests=[]
    for i,k in enumerate(pk):
        if support[i]:
            t=compare(float(stats[k]['score']),float(observed[i]),'direction');tests.append(dict(pert_id=k[0],fold=k[1],**t));assert t['numeric_pass'] and t['decision_pass']
    csvwrite(O/'QA/CRITICAL_OBSERVED_INDEPENDENT_RECOMPUTATION.csv',tests)
    B=99999;seed=171260901;rng=np.random.Generator(np.random.PCG64(seed))
    strata=[]
    for label in sorted({r['null_stratum'] for r in genes}):strata.append(np.asarray([i for i,r in enumerate(genes) if r['null_stratum']==label],dtype=np.int64))
    expected_hashes=csvread(O/'NULL/CONTINUOUS_iteration_member_hashes.csv');assert len(expected_hashes)==B
    pids=sorted({r['pert_id'] for r in allrows});pindex={p:i for i,p in enumerate(pids)}
    first={f:np.load(O/'NULL'/f'G2_R2_CAP8_{f}.npy',mmap_mode='r') for f in ['DISCOVERY','VALIDATION']}
    valid=np.zeros(len(pk),dtype=np.int64);lower=valid.copy();upper=valid.copy()
    dest=O/'QA/CRITICAL_NULL_CONTINUOUS_COMPARISONS_FLOAT64.npy';assert not dest.exists(),'Existing independent run requires explicit verified checkpoint recovery'
    pairs_file=open(dest,'wb')
    np.lib.format.write_array_header_2_0(pairs_file,dict(descr=np.dtype(np.float64).str,fortran_order=False,shape=(B,len(pk),4)))
    engine_hash=sha(Path(__file__));input_hash=sha(O/'EXECUTION_INPUT_HASH_BEFORE.csv')
    for lo in range(0,B,32):
        n=min(32,B-lo);weights=np.empty((8227,n),dtype=np.float64)
        for j in range(n):
            digest=hashlib.sha256()
            for idx in strata:
                selected_order=rng.permutation(idx);weights[idx,j]=w[selected_order];digest.update(selected_order.astype('<i8').tobytes())
            assert digest.hexdigest()==expected_hashes[lo+j]['permutation_member_sha256'],'NULL_MEMBER_SEQUENCE_DISCORDANCE'
        norms=np.sqrt(np.einsum('ij,ij->j',weights,weights,optimize=False))
        s=np.empty((len(rows),n),dtype=np.float64)
        for start in range(0,len(rows),256):
            s[start:start+256]=np.einsum('ij,jk->ik',X[start:start+256],weights,optimize=False)/(zn[start:start+256,None]*norms[None,:])
        ref=chain(s);ref[~support]=np.nan
        primary=np.stack([first[k[1]][lo:lo+n,pindex[k[0]]] for k in pk])
        same_state=np.array_equal(np.isnan(primary),np.isnan(ref)) and np.array_equal(np.isposinf(primary),np.isposinf(ref)) and np.array_equal(np.isneginf(primary),np.isneginf(ref))
        finite=np.isfinite(primary)&np.isfinite(ref);delta=np.abs(primary-ref)
        direction_ok=np.array_equal(np.sign(primary[finite]),np.sign(ref[finite]))
        tolerance_ok=np.all(delta[finite]<=1e-12+1e-10*np.abs(ref[finite]))
        with np.errstate(divide='ignore',invalid='ignore'):relative=np.where(ref==0,np.where(delta==0,0,np.inf),delta/np.abs(ref))
        comparison_block=np.stack([primary.T,ref.T,delta.T,relative.T],axis=2).astype(np.float64)
        comparison_block.tofile(pairs_file);pairs_file.flush()
        if not(same_state and direction_ok and tolerance_ok):
            js(O/'QA/INDEPENDENT_FULL_NULL_DISCORDANCE.json',dict(utc=now(),start=lo,end=lo+n,state_match=same_state,direction_match=direction_ok,tolerance_match=bool(tolerance_ok),status='DIRECTION_DISCORDANCE' if not direction_ok else 'NUMERICAL_OR_STATE_DISCORDANCE'))
            raise RuntimeError('Stop final classification: independent full null discordance')
        valid+=np.isfinite(ref).sum(axis=1);lower+=(ref<=observed[:,None]).sum(axis=1);upper+=(ref>=observed[:,None]).sum(axis=1)
        assert sha(Path(__file__))==engine_hash
        js(O/'CHECKPOINTS/INDEPENDENT_CRITICAL_MAIN_NULL.json',dict(utc=now(),next_iteration=lo+n,chunk_start=lo,chunk_end=lo+n,B=B,seed=seed,rng_state=rng.bit_generator.state,valid=valid.tolist(),lower=lower.tolist(),upper=upper.tolist(),engine_sha256=engine_hash,input_manifest_sha256=input_hash,output_chunk_sha256=hashlib.sha256(comparison_block.tobytes()).hexdigest(),status='COMPLETE'))
        if lo%3200==0:print('INDEPENDENT_CRITICAL_NULL',lo+n,'/',B,flush=True)
    pairs_file.close()
    assert rng.bit_generator.state==jr(O/'CHECKPOINTS/CONTINUOUS.json')['rng_state']
    count_results=[]
    for i,k in enumerate(pk):
        expected=stats[k];ok=(int(valid[i]),int(lower[i]),int(upper[i]))==(int(expected['null_valid_iterations']),int(expected['lower_extreme_count']),int(expected['upper_extreme_count']))
        count_results.append(dict(pert_id=k[0],fold=k[1],valid_path2=int(valid[i]),lower_path2=int(lower[i]),upper_path2=int(upper[i]),valid_path1=expected['null_valid_iterations'],lower_path1=expected['lower_extreme_count'],upper_path1=expected['upper_extreme_count'],exact_counts_match=ok))
    csvwrite(O/'QA/INDEPENDENT_COMPLETE_NULL_EXTREME_COUNTS.csv',count_results)
    assert all(r['exact_counts_match'] for r in count_results),'EXACT_NULL_EXTREME_COUNT_DISCORDANCE'
    csvwrite(O/'QA/CRITICAL_NULL_COMPARISON_AXIS.csv',[dict(index=i,pert_id=k[0],fold=k[1]) for i,k in enumerate(pk)])
    js(O/'QA/INDEPENDENT_FULL_NULL_VALIDATION.json',dict(utc=now(),status='PASS',compounds=len(requested),folds=len(pk),B=B,seed=seed,permutation_hashes_match=True,counts_exact=True,continuous_comparison_sha256=sha(dest),comparison_last_axis=['path1_comparison_value','path2_reference_value','absolute_difference','relative_difference']))
    event('INDEPENDENT_FULL_NULL_COUNTS_VERIFIED',critical_compounds=len(requested),iterations=B)

if __name__=='__main__':run()
