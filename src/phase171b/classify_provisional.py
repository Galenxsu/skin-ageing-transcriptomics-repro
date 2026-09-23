from engine171b import *
import re

LABELS={k:v for k,v in [('A','ROBUST_EXPLORATORY_REVERSAL_CANDIDATE'),('B','DISCOVERY_ONLY_EXPLORATORY_SIGNAL'),('C','CONTEXT_DEPENDENT_EXPLORATORY_SIGNAL'),('D','RANKED_NONSIGNIFICANT_LEAD'),('E','NO_EVIDENCE_OF_REVERSAL'),('F','MIMETIC_DIRECTION_SIGNAL'),('G','NOT_EVALUABLE')]}

def number(v):return float(v) if v not in ('NA','',None) else float('nan')
def yes(v):return v is True or v=='True'

def run():
    rows,genes,pids,queries,strata=metadata();obs=np.load(O/'RESULTS/observed_compound_scores.npy');names=jr(O/'RESULTS/query_order.json')
    stat={(r['query_id'],r['fold'],r['pert_id']):r for r in csvread(O/'RESULTS/all_query_statistics.csv')}
    gates={(r['pert_id'],r['fold']):r for r in csvread(O/'RESULTS/robustness_gate_inputs.csv')}
    # Only frozen identity validity and historical flags are used, never names or mechanism.
    identity={r['pert_id']:{k:r[k] for k in ['inchi_key','structure_unique','flags']} for r in csvread(A/'eligible_compound_universe.csv')}
    historical={r['pert_id'] for r in rows if r['historically_scored_compound']=='True'}
    ds={p:stat['G2_R2_CAP8','DISCOVERY',p] for p in pids};vs={p:stat['G2_R2_CAP8','VALIDATION',p] for p in pids}
    rank_order=sorted(pids,key=lambda p:(number(ds[p]['raw_P_lower']) if ds[p]['raw_P_lower']!='NA' else 1.,number(ds[p]['score']) if ds[p]['score']!='NA' else math.inf,p))
    ranking={p:i+1 for i,p in enumerate(rank_order)};out=[];critical=set(historical)
    for pi,p in enumerate(pids):
        d,v=ds[p],vs[p];dg,vg=gates.get((p,'DISCOVERY'),{}),gates.get((p,'VALIDATION'),{})
        de=yes(d['null_evaluable']) and yes(dg.get('full_support',False));ve=yes(v['null_evaluable']) and yes(vg.get('full_support',False))
        discovery=de and number(d['score'])<0 and number(d['BH_FDR_lower'])<.05
        validation=ve and number(v['score'])<0 and number(v['BH_FDR_lower'])<.05
        mimic=de and number(d['score'])>0 and number(d['BH_FDR_upper'])<.05
        cell=all(int(g.get('base_cells',0))>=2 and int(g.get('negative_base_cells',0))>=1 and yes(g.get('base_cell_LOO_pass',False)) for g in [dg,vg])
        context=all(yes(g.get('context_LOO_pass',False)) for g in [dg,vg])
        dose=all(yes(g.get('dose_LOO_pass',False)) and yes(g.get('low_dose_pass',False)) for g in [dg,vg])
        timing=all(yes(g.get('time_LOO_pass',False)) for g in [dg,vg])
        hq=all(yes(g.get('HQ_pass',False)) for g in [dg,vg])
        numeric=all(float(g.get('numerically_invalid_fraction',1))<=.05 for g in [dg,vg])
        ident=bool(re.fullmatch(r'[A-Z]{14}-[A-Z]{10}-[A-Z]',identity[p]['inchi_key'])) and yes(identity[p]['structure_unique'])
        sensitivity=all(np.isfinite(obs[ix,names.index(q)]) and obs[ix,names.index(q)]<0 for ix in [pi,1792+pi] for q in ['G2_R1']+[f'G2_E{n}' for n in [25,50,75,100,150,200,300]])
        flags={'MAIN_DISCOVERY_FDR_OR_DIRECTION_NOT_MET':not discovery,'HELDOUT_FDR_SUPPORT_OR_DIRECTION_NOT_MET':not validation,'CELL_LINE_SUPPORT_OR_LEAVE_ONE_CELL_OUT_LIMIT':not cell,'CONTEXT_DEPENDENCE_OR_INSUFFICIENT_LOO_SUPPORT':not context,'DOSE_DEPENDENCE_OR_INSUFFICIENT_LOW_DOSE_SUPPORT':not dose,'TIME_DEPENDENCE_OR_INSUFFICIENT_TIME_SUPPORT':not timing,'HQ_ROBUSTNESS_NOT_MET':not hq,'NUMERICAL_INVALID_PROFILE_FRACTION_LIMIT':not numeric,'IDENTITY_NOT_UNIQUE_OR_NOT_CONFIRMED':not ident,'SIGNATURE_REPRESENTATION_NOT_STABLE':not sensitivity,'HISTORICAL_RESULT_EXPOSURE':p in historical}
        reasons=sorted(k for k,b in flags.items() if b)
        allhigh=not reasons;contextlimit=not(cell and context and dose and timing)
        if not de:code='G'
        elif mimic:code='F'
        elif allhigh:code='A'
        elif discovery and contextlimit:code='C'
        elif discovery:code='B'
        elif number(d['score'])<0 and ranking[p]<=18:code='D'
        else:code='E'
        # Independent rule-table evaluation, not reading the first classification.
        predicate=[('G',not de),('F',mimic),('A',discovery and validation and cell and context and dose and timing and hq and numeric and ident and sensitivity and p not in historical),('C',discovery and (not cell or not context or not dose or not timing)),('B',discovery),('D',number(d['score'])<0 and ranking[p]<=18),('E',True)]
        ref=next(k for k,b in predicate if b);assert ref==code,'CLASSIFICATION_IMPLEMENTATION_DISCORDANCE'
        boundary=any(.04<=number(s['BH_FDR_'+tail])<=.06 for s in [d,v] for tail in ['lower','upper'])
        if code in 'ABCD' or boundary:critical.add(p)
        out.append(dict(pert_id=p,discovery_rank=ranking[p],discovery_score=d['score'],discovery_P=d['raw_P_lower'],discovery_FDR=d['BH_FDR_lower'],heldout_score=v['score'],heldout_P=v['raw_P_lower'],heldout_FDR=v['BH_FDR_lower'],original_statistical_class='MAIN_NOT_EVALUABLE' if not de else 'MAIN_MIMETIC_FDR_PASS' if mimic else 'MAIN_REVERSAL_FDR_PASS' if discovery else 'MAIN_NO_FDR_REVERSAL',discovery_status=discovery,heldout_status=validation,signature_robustness=sensitivity,cell_line_support=cell,context_dependence=not context,dose_dependence=not dose,time_dependence=not timing,HQ_pass=hq,identity_status=ident,toxicity_interpretability='FUNCTIONAL_TOXICITY_NOT_ASSESSABLE',historical_flag='HISTORICALLY_EXAMINED_COMPOUND' if p in historical else 'NOT_HISTORICALLY_EXAMINED',independent_validation_status='INDEPENDENT_VALIDATION_NOT_ELIGIBLE' if p in historical else 'INTERNAL_HELDOUT_ONLY_NOT_INDEPENDENT_EXPERIMENT',provisional_classification=LABELS[code],provisional_code=code,downgrade_reasons='|'.join(reasons) if reasons else 'NONE',independent_decision_code_match=True,publication_status='NOT_RELEASED_PENDING_FULL_INDEPENDENT_NULL_COUNT_AUDIT'))
    csvwrite(O/'RESULTS/provisional_blinded_classification.csv',out)
    csvwrite(O/'QA/CRITICAL_COMPOUND_INDEPENDENT_NULL_AUDIT_SET.csv',[dict(pert_id=p,scope='COMPLETE_MAIN_R2_DISCOVERY_AND_VALIDATION_NULL',reason='HISTORICAL6' if p in historical else 'PROVISIONAL_ABCD_OR_BOUNDARY') for p in sorted(critical)])
    event('PROVISIONAL_CLASSIFICATION_FOR_AUDIT_SELECTION_ONLY',compounds=len(out),critical_audit_compounds=len(critical),name_unblinding=False)

if __name__=='__main__':run()
