from engine171b import *
from numeric_rules import compare
import time
from threadpoolctl import threadpool_limits

def run():
    rows,genes,pids,queries,strata=metadata();h=hierarchy(rows,pids)
    X=np.load(O/'MATRIX/frozen_107201x8227_float64.npy',mmap_mode='r')
    assert sha(O/'MATRIX/frozen_107201x8227_float64.npy')==jr(O/'RESULTS/EXPRESSION_INTEGRITY_SUMMARY.json')['matrix_sha256']
    fixed=csvread(O/'QA/FROZEN_INDEPENDENT_AUDIT_SAMPLE.csv');sids={r['sig_id'] for r in fixed}
    ix=[i for i,r in enumerate(rows) if r['sig_id'] in sids]
    counts=defaultdict(int);acounts=defaultdict(int)
    for r in rows:counts[r['pert_id']]+=1
    for i in ix:acounts[rows[i]['pert_id']]+=1
    complete={p for p,n in acounts.items() if n==counts[p]}
    t=time.perf_counter()
    with threadpool_limits(limits=8):
        scores=vector_scores(X,queries)
    contexts,bases,compounds=aggregate(scores,h)
    tests=[];names=list(queries)
    for qi,q in enumerate(names):
        independent={i:scalar_score(X[i],queries[q]) for i in ix}
        for i in ix:
            result=compare(float(scores[i,qi]),independent[i],decision='direction')
            tests.append(dict(level='PROFILE',query=q,key=rows[i]['sig_id'],**result))
        ic,ib,ip=independent_aggregate([(rows[i],independent[i]) for i in ix if rows[i]['pert_id'] in complete])
        for level,keys,array,reference in [('CONTEXT',h['context_keys'],contexts,ic),('BASE_CELL',h['base_keys'],bases,ib),('COMPOUND',h['compound_keys'],compounds,ip)]:
            for i,k in enumerate(keys):
                if k in reference:tests.append(dict(level=level,query=q,key='|'.join(map(str,k)),**compare(float(array[i,qi]),reference[k],decision='direction')))
        print('INDEPENDENT_OBSERVED_QUERY_COMPLETE',q,flush=True)
    csvwrite(O/'QA/INDEPENDENT_OBSERVED_NUMERICAL_COMPARISONS.csv',tests)
    failed=[r for r in tests if not (r['numeric_pass'] and r['decision_pass'])]
    summary=dict(utc=now(),comparisons=len(tests),failed=len(failed),fixed_signatures=len(ix),complete_audit_compounds=len(complete),runtime_seconds=time.perf_counter()-t)
    js(O/'QA/INDEPENDENT_OBSERVED_SUMMARY.json',summary)
    if failed:
        event('INDEPENDENT_OBSERVED_FAILED',**summary)
        raise RuntimeError('Independent observed scoring disagreement; no null or classification permitted')
    np.save(O/'RESULTS/observed_profile_scores.npy',scores)
    np.save(O/'RESULTS/observed_context_scores.npy',contexts)
    np.save(O/'RESULTS/observed_base_scores.npy',bases)
    np.save(O/'RESULTS/observed_compound_scores.npy',compounds)
    js(O/'RESULTS/query_order.json',names)
    csvwrite(O/'RESULTS/compound_fold_order.csv',[dict(index=i,pert_id=k[0],fold=k[1]) for i,k in enumerate(h['compound_keys'])])
    csvwrite(O/'RESULTS/context_order.csv',[dict(index=i,pert_id=k[0],fold=k[1],base_cell_id=k[2],cell_id=k[3],dose_uM=k[4],time_h=k[5]) for i,k in enumerate(h['context_keys'])])
    csvwrite(O/'RESULTS/base_cell_order.csv',[dict(index=i,pert_id=k[0],fold=k[1],base_cell_id=k[2]) for i,k in enumerate(h['base_keys'])])
    event('OBSERVED_SCORES_INDEPENDENTLY_VERIFIED',**summary)
    print(json.dumps(summary),flush=True)

if __name__=='__main__':run()
