from engine171b import *

def run():
    rows,genes,pids,queries,strata=metadata();names=jr(O/'RESULTS/query_order.json')
    scores=np.load(O/'RESULTS/observed_profile_scores.npy');main=scores[:,names.index('G2_R2_CAP8')]
    groups=defaultdict(list)
    for i,r in enumerate(rows):groups[r['pert_id'],r['fold']].append(i)
    def evaluate(ix):
        valid=[i for i in ix if np.isfinite(main[i])]
        cg=defaultdict(list)
        for i in valid:
            r=rows[i];cg[(r['base_cell_id'],r['cell_id'],float(r['dose_uM']),float(r['time_h']))].append(float(main[i]))
        bg=defaultdict(list)
        for k,v in cg.items():bg[k[0]].append(float(np.median(v)))
        base={k:float(np.median(v)) for k,v in bg.items()}
        raw=float(np.median(list(base.values()))) if base else float('nan')
        return dict(valid_profiles=len(valid),contexts=len(cg),base_cells=len(base),negative_base_cells=sum(v<0 for v in base.values()),score=raw if np.isfinite(raw) else 'NA',support=len(valid)>=3 and len(cg)>=2,negative=np.isfinite(raw) and raw<0)
    records=[];gate=[]
    for (p,f),ix in groups.items():
        full=evaluate(ix);hq=evaluate([i for i in ix if rows[i]['HQ']=='True']);low=evaluate([i for i in ix if float(rows[i]['dose_uM'])<10])
        local=[]
        for kind,key in [('CONTEXT',lambda i:(rows[i]['cell_id'],float(rows[i]['dose_uM']),float(rows[i]['time_h']))),('BASE_CELL',lambda i:rows[i]['base_cell_id']),('DOSE',lambda i:float(rows[i]['dose_uM'])),('TIME',lambda i:float(rows[i]['time_h']))]:
            levels=sorted({key(i) for i in ix})
            for level in levels:
                ev=evaluate([i for i in ix if key(i)!=level]);rec=dict(pert_id=p,fold=f,analysis='LEAVE_ONE_'+kind+'_OUT',removed_level=str(level),**ev);records.append(rec);local.append(rec)
                if kind=='TIME':records.append(dict(pert_id=p,fold=f,analysis='TIME_STRATUM',removed_level=str(level),**evaluate([i for i in ix if key(i)==level])))
        for label,ev in [('FULL_PROFILE',full),('HQ',hq),('REMOVE_GE10_UM',low)]:records.append(dict(pert_id=p,fold=f,analysis=label,removed_level='NOT_APPLICABLE',**ev))
        one=lambda kind:all(r['support'] and r['negative'] for r in local if r['analysis']=='LEAVE_ONE_'+kind+'_OUT')
        nt=len({float(rows[i]['time_h']) for i in ix});nd=len({float(rows[i]['dose_uM']) for i in ix})
        time_each=all(evaluate([i for i in ix if float(rows[i]['time_h'])==t])['negative'] for t in {float(rows[i]['time_h']) for i in ix})
        gate.append(dict(pert_id=p,fold=f,full_support=full['support'],profiles=full['valid_profiles'],contexts=full['contexts'],base_cells=full['base_cells'],negative_base_cells=full['negative_base_cells'],numerically_invalid_fraction=1-full['valid_profiles']/len(ix),HQ_pass=hq['support'] and hq['negative'],context_LOO_pass=one('CONTEXT'),base_cell_LOO_pass=one('BASE_CELL'),dose_LOO_pass=nd>=2 and one('DOSE'),time_LOO_pass=nt>=2 and one('TIME') and time_each,low_dose_pass=low['support'] and low['negative'],dose_levels=nd,time_levels=nt))
    csvwrite(O/'RESULTS/robustness_recalculations.csv',records)
    csvwrite(O/'RESULTS/robustness_gate_inputs.csv',gate)
    event('ROBUSTNESS_FROM_FROZEN_RULES_COMPUTED',rows=len(records),compound_fold_gates=len(gate))

if __name__=='__main__':run()
