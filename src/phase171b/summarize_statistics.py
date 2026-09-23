from engine171b import *
from numeric_rules import compare
from fractions import Fraction

def bh_numpy(p):
    p=np.asarray(p,dtype=np.float64);order=np.argsort(p,kind='stable');n=len(p)
    ranked=p[order]*n/np.arange(1,n+1,dtype=np.float64)
    adjusted=np.minimum(1,np.minimum.accumulate(ranked[::-1])[::-1]);out=np.empty(n);out[order]=adjusted
    return out

def bh_reference(p):
    n=len(p);order=sorted(range(n),key=lambda i:(float(p[i]),i));out=[0.]*n;current=1.
    for rank in range(n,0,-1):
        i=order[rank-1];current=min(current,float(p[i])*n/rank);out[i]=current
    return out

def run():
    from null171b import plans
    for name,qs,B,seed in plans():
        cp=jr(O/'CHECKPOINTS'/f'{name}.json');assert cp['next_iteration']==B,'Full frozen null not complete: no P/FDR calculation'
    rows,genes,pids,queries,strata=metadata();h=hierarchy(rows,pids)
    obs=np.load(O/'RESULTS/observed_compound_scores.npy');names=jr(O/'RESULTS/query_order.json')
    stats={};audit=[];p_comparisons=[]
    for name,qs,B,seed in plans():
        cp=jr(O/'CHECKPOINTS'/f'{name}.json')
        for q in qs:
            for fold in (['DISCOVERY','VALIDATION'] if q=='G2_R2_CAP8' else ['DISCOVERY']):
                ni=np.load(O/'NULL'/f'{q}_{fold}.npy',mmap_mode='r');assert ni.shape==(B,1792)
                offset=0 if fold=='DISCOVERY' else 1792;o=obs[offset:offset+1792,names.index(q)]
                valid=np.zeros(1792,dtype=np.int64);lower=valid.copy();upper=valid.copy()
                for lo in range(0,B,2048):
                    block=np.asarray(ni[lo:lo+2048]);valid+=np.isfinite(block).sum(axis=0)
                    lower+=(block<=o[None,:]).sum(axis=0);upper+=(block>=o[None,:]).sum(axis=0)
                for i,p in enumerate(pids):
                    good=valid[i]==B and np.isfinite(o[i]);pl=(int(lower[i])+1)/(B+1) if good else None;pu=(int(upper[i])+1)/(B+1) if good else None
                    stats[q,fold,p]=dict(query_id=q,fold=fold,pert_id=p,score=float(o[i]) if np.isfinite(o[i]) else 'NA',B=B,seed=seed,null_valid_iterations=int(valid[i]),lower_extreme_count=int(lower[i]),upper_extreme_count=int(upper[i]),raw_P_lower=pl if good else 'NA',raw_P_upper=pu if good else 'NA',null_evaluable=good)
                    for side,value,count in [('lower',pl,int(lower[i])),('upper',pu,int(upper[i]))]:
                        reference=float(Fraction(1+count,B+1)) if good else 'NA'
                        comparison=compare(value if good else 'NA',reference)
                        assert comparison['numeric_pass'] and comparison['decision_pass']
                        p_comparisons.append(dict(query=q,fold=fold,pert_id=p,tail=side,integer_extreme_count=count,integer_denominator=B+1,**comparison))
                audit.append(dict(query=q,fold=fold,B=B,complete_null_compounds=int(np.sum(valid==B)),incomplete_null_compounds=int(np.sum(valid!=B)),minimum_resolvable_P=1/(B+1),seed=seed))
    registry=csvread(A/'planned_test_family_registry.csv');families=defaultdict(list)
    for r in registry:families[r['family_id']].append(r)
    comparisons=[];familyaudit=[]
    for family,tests in families.items():
        tail=tests[0]['tail'];side={'REVERSAL_LOWER':'lower','MIMETIC_UPPER':'upper'}[tail]
        # Registry tail spelling is asserted, never inferred from observed scores.
        assert tail in ('REVERSAL_LOWER','MIMETIC_UPPER'),tail
        p=[stats[r['query_id'],r['fold'],r['pert_id']]['raw_P_'+side] for r in tests]
        effective=[1. if x=='NA' else float(x) for x in p]
        q=bh_numpy(effective);reference=bh_reference(effective)
        for i,r in enumerate(tests):
            c=compare(float(q[i]),reference[i],'significance');comparisons.append(dict(test_id=r['test_id'],family=family,**c))
            assert c['numeric_pass'] and c['decision_pass'],'BH independent implementation mismatch'
            result=stats[r['query_id'],r['fold'],r['pert_id']];result['BH_FDR_'+side]=float(q[i]);result['family_'+side]=family
        familyaudit.append(dict(family_id=family,planned_tests=len(tests),valid_tests=sum(x!='NA' for x in p),NA_P1_slots=sum(x=='NA' for x in p),BH_method='COMPLETE_FAMILY_BH',alpha=.05,independent_BH_pass=True))
    assert len(registry)==53760 and len(families)==12
    csvwrite(O/'RESULTS/all_query_statistics.csv',list(stats.values()))
    csvwrite(O/'RESULTS/empirical_null_summary.csv',audit)
    csvwrite(O/'RESULTS/multiple_testing_family_audit.csv',familyaudit)
    csvwrite(O/'QA/INDEPENDENT_BH_COMPARISON.csv',comparisons)
    csvwrite(O/'QA/INDEPENDENT_EMPIRICAL_P_COMPARISON.csv',p_comparisons)
    event('P_BH_COMPUTED_FROM_COMPLETE_NULL',families=12,directional_tests=53760,final_classification_released=False)

if __name__=='__main__':run()
