from common171b import *
import numpy as np
import math, ast
from collections import defaultdict
from numeric_rules import compare

def median(values):
    v=sorted(values);n=len(v)
    return v[n//2] if n%2 else (v[n//2-1]+v[n//2])/2 if n else float('nan')

def run():
    assert jr(O/'QA/FULL_INDEPENDENT_OBSERVED_VALIDATION.json')['pass_all']
    meta=[r for r in csvread(A/'discovery_validation_assignment_preview.csv') if r['metadata_eligible']=='True']
    q=jr(O/'RESULTS/query_order.json').index('G2_R2_CAP8')
    score=np.load(O/'QA/FULL_INDEPENDENT_profile_reference.npy',mmap_mode='r')[:,q]
    group=defaultdict(list)
    for i,r in enumerate(meta):group[r['pert_id'],r['fold']].append(i)
    def ev(ids):
        valid=[i for i in ids if math.isfinite(float(score[i]))]
        contexts=defaultdict(list)
        for i in valid:
            r=meta[i];contexts[(r['base_cell_id'],r['cell_id'],float(r['dose_uM']),float(r['time_h']))].append(float(score[i]))
        cells=defaultdict(list)
        for key,values in contexts.items():cells[key[0]].append(median(values))
        bases=[median(v) for v in cells.values()];s=median(bases)
        return dict(valid_profiles=len(valid),contexts=len(contexts),base_cells=len(bases),negative_base_cells=sum(v<0 for v in bases),score=s if math.isfinite(s) else 'NA',support=len(valid)>=3 and len(contexts)>=2,negative=math.isfinite(s) and s<0)
    results=[];audit=[];refs={};by_group=defaultdict(list)
    for old in csvread(O/'RESULTS/robustness_recalculations.csv'):
        p,f,a,lev=(old[k] for k in ['pert_id','fold','analysis','removed_level']);ids=group[p,f]
        if a=='HQ':chosen=[i for i in ids if meta[i]['HQ']=='True']
        elif a=='REMOVE_GE10_UM':chosen=[i for i in ids if float(meta[i]['dose_uM'])<10]
        elif a=='FULL_PROFILE':chosen=ids
        elif a=='TIME_STRATUM':chosen=[i for i in ids if float(meta[i]['time_h'])==float(lev)]
        else:
            kind=a.removeprefix('LEAVE_ONE_').removesuffix('_OUT')
            def key(i):
                r=meta[i]
                if kind=='CONTEXT':return (r['cell_id'],float(r['dose_uM']),float(r['time_h']))
                if kind=='BASE_CELL':return r['base_cell_id']
                return float(r['dose_uM' if kind=='DOSE' else 'time_h'])
            target=ast.literal_eval(lev) if kind=='CONTEXT' else lev if kind=='BASE_CELL' else float(lev)
            chosen=[i for i in ids if key(i)!=target]
        out=dict(pert_id=p,fold=f,analysis=a,removed_level=lev,**ev(chosen));results.append(out);refs[p,f,a,lev]=out;by_group[p,f].append(out)
        for field in ['valid_profiles','contexts','base_cells','negative_base_cells','support','negative']:
            assert str(out[field])==old[field],f'EXACT_ROBUSTNESS_DISCORDANCE {p} {f} {a} {field}'
        x1=float(old['score']) if old['score']!='NA' else None;x2=float(out['score']) if out['score']!='NA' else None
        c=compare(x1,x2,'direction')
        assert c['numeric_pass'] and c['decision_pass'],f'ROBUSTNESS_NUMERICAL_DISCORDANCE {p} {f} {a}'
        audit.append(dict(pert_id=p,fold=f,analysis=a,removed_level=lev,**c))
    refg=[]
    for old in csvread(O/'RESULTS/robustness_gate_inputs.csv'):
        p,f=old['pert_id'],old['fold'];ids=group[p,f]
        rr=by_group[p,f]
        named={r['analysis']:r for r in rr if r['analysis'] in ['FULL_PROFILE','HQ','REMOVE_GE10_UM']}
        full=named['FULL_PROFILE'];hq=named['HQ'];low=named['REMOVE_GE10_UM']
        def loo(kind):return all(r['support'] and r['negative'] for r in rr if r['analysis']=='LEAVE_ONE_'+kind+'_OUT')
        nd=len(set(float(meta[i]['dose_uM']) for i in ids));nt=len(set(float(meta[i]['time_h']) for i in ids))
        out=dict(pert_id=p,fold=f,full_support=full['support'],profiles=full['valid_profiles'],contexts=full['contexts'],base_cells=full['base_cells'],negative_base_cells=full['negative_base_cells'],numerically_invalid_fraction=1-full['valid_profiles']/len(ids),HQ_pass=hq['support'] and hq['negative'],context_LOO_pass=loo('CONTEXT'),base_cell_LOO_pass=loo('BASE_CELL'),dose_LOO_pass=nd>=2 and loo('DOSE'),time_LOO_pass=nt>=2 and loo('TIME') and all(r['negative'] for r in rr if r['analysis']=='TIME_STRATUM'),low_dose_pass=low['support'] and low['negative'],dose_levels=nd,time_levels=nt)
        for field,value in out.items():assert str(value)==old[field],f'GATE_DISCORDANCE {p} {f} {field}'
        refg.append(out)
    csvwrite(O/'QA/INDEPENDENT_ROBUSTNESS_RECALCULATIONS.csv',results)
    csvwrite(O/'QA/INDEPENDENT_ROBUSTNESS_GATE_INPUTS.csv',refg)
    csvwrite(O/'QA/INDEPENDENT_ROBUSTNESS_COMPARISONS.csv',audit)
    js(O/'QA/INDEPENDENT_ROBUSTNESS_VALIDATION.json',dict(utc=now(),pass_all=True,recalculations=len(results),compound_fold_gates=len(refg),reference_source='INDEPENDENT_PROFILE_RECALCULATION_FROM_FROZEN_EXPRESSION',aggregation='SEPARATE_SORTED_MEDIAN_IMPLEMENTATION',counts_directions_and_gate_flags='EXACT_MATCH'))
    event('INDEPENDENT_ROBUSTNESS_VALIDATED',recalculations=len(results),gates=len(refg))

if __name__=='__main__':run()
