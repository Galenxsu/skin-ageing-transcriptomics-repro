from common171b import *
import numpy as np, math, re

def yes(x):return x is True or x=='True'
def num(x):return float(x) if x not in ('NA',None,'') else math.nan

def run():
    assert jr(O/'QA/INDEPENDENT_ROBUSTNESS_VALIDATION.json')['pass_all']
    assert jr(O/'QA/FULL_INDEPENDENT_OBSERVED_VALIDATION.json')['pass_all']
    nullqa=jr(O/'QA/INDEPENDENT_FULL_NULL_VALIDATION.json')
    assert nullqa['status']=='PASS',nullqa
    meta=[r for r in csvread(A/'discovery_validation_assignment_preview.csv') if r['metadata_eligible']=='True']
    pids=sorted({r['pert_id'] for r in meta});assert len(pids)==1792
    historic={r['pert_id'] for r in meta if r['historically_scored_compound']=='True'};assert len(historic)==6
    identity={r['pert_id']:r for r in csvread(A/'eligible_compound_universe.csv')}
    stats={(r['query_id'],r['fold'],r['pert_id']):r for r in csvread(O/'RESULTS/all_query_statistics.csv')}
    gates={(r['pert_id'],r['fold']):r for r in csvread(O/'QA/INDEPENDENT_ROBUSTNESS_GATE_INPUTS.csv')}
    values=np.load(O/'QA/FULL_INDEPENDENT_compound_reference.npy');names=jr(O/'RESULTS/query_order.json')
    prior={r['pert_id']:r for r in csvread(O/'RESULTS/provisional_blinded_classification.csv')}
    def ds(p):return stats['G2_R2_CAP8','DISCOVERY',p]
    ranking=sorted(pids,key=lambda p:(num(ds(p)['raw_P_lower']) if ds(p)['raw_P_lower']!='NA' else 1.,float(values[pids.index(p),0]) if np.isfinite(values[pids.index(p),0]) else math.inf,p))
    ranks={p:i+1 for i,p in enumerate(ranking)}
    output=[];checks=[]
    for pi,p in enumerate(pids):
        d=ds(p);v=stats['G2_R2_CAP8','VALIDATION',p]
        gg=[gates.get((p,f),{}) for f in ['DISCOVERY','VALIDATION']]
        eligible=[yes(s['null_evaluable']) and yes(g.get('full_support',False)) for s,g in zip([d,v],gg)]
        reversal=[e and float(values[pi+1792*j,0])<0 and num(s['BH_FDR_lower'])<.05 for j,(e,s) in enumerate(zip(eligible,[d,v]))]
        mimic=eligible[0] and float(values[pi,0])>0 and num(d['BH_FDR_upper'])<.05
        passes={
          'MAIN_DISCOVERY_FDR_OR_DIRECTION_NOT_MET':reversal[0],
          'HELDOUT_FDR_SUPPORT_OR_DIRECTION_NOT_MET':reversal[1],
          'CELL_LINE_SUPPORT_OR_LEAVE_ONE_CELL_OUT_LIMIT':all(int(g.get('base_cells',0))>=2 and int(g.get('negative_base_cells',0))>=1 and yes(g.get('base_cell_LOO_pass',False)) for g in gg),
          'CONTEXT_DEPENDENCE_OR_INSUFFICIENT_LOO_SUPPORT':all(yes(g.get('context_LOO_pass',False)) for g in gg),
          'DOSE_DEPENDENCE_OR_INSUFFICIENT_LOW_DOSE_SUPPORT':all(yes(g.get('dose_LOO_pass',False)) and yes(g.get('low_dose_pass',False)) for g in gg),
          'TIME_DEPENDENCE_OR_INSUFFICIENT_TIME_SUPPORT':all(yes(g.get('time_LOO_pass',False)) for g in gg),
          'HQ_ROBUSTNESS_NOT_MET':all(yes(g.get('HQ_pass',False)) for g in gg),
          'NUMERICAL_INVALID_PROFILE_FRACTION_LIMIT':all(float(g.get('numerically_invalid_fraction',1))<=.05 for g in gg),
          'IDENTITY_NOT_UNIQUE_OR_NOT_CONFIRMED':bool(re.fullmatch('[A-Z]{14}-[A-Z]{10}-[A-Z]',identity[p]['inchi_key'])) and yes(identity[p]['structure_unique']),
          'SIGNATURE_REPRESENTATION_NOT_STABLE':all(np.isfinite(values[row,names.index(q)]) and values[row,names.index(q)]<0 for row in [pi,pi+1792] for q in ['G2_R1','G2_E25','G2_E50','G2_E75','G2_E100','G2_E150','G2_E200','G2_E300']),
          'HISTORICAL_RESULT_EXPOSURE':p not in historic}
        reason='|'.join(sorted(k for k,value in passes.items() if not value)) or 'NONE'
        restrictions=any(not passes[k] for k in ['CELL_LINE_SUPPORT_OR_LEAVE_ONE_CELL_OUT_LIMIT','CONTEXT_DEPENDENCE_OR_INSUFFICIENT_LOO_SUPPORT','DOSE_DEPENDENCE_OR_INSUFFICIENT_LOW_DOSE_SUPPORT','TIME_DEPENDENCE_OR_INSUFFICIENT_TIME_SUPPORT'])
        code='E'
        if float(values[pi,0])<0 and ranks[p]<=18:code='D'
        if reversal[0]:code='C' if restrictions else 'B'
        if all(passes.values()):code='A'
        if mimic:code='F'
        if not eligible[0]:code='G'
        old=prior[p]
        assert code==old['provisional_code'],'CLASSIFICATION_DISCORDANCE '+p
        assert reason==old['downgrade_reasons'],'DOWNGRADE_REASON_DISCORDANCE '+p
        assert ranks[p]==int(old['discovery_rank']),'RANK_BOUNDARY_DISCORDANCE '+p
        out={k:v for k,v in old.items() if not k.startswith('provisional_') and k!='publication_status'}
        out.update(final_code=code,final_classification=old['provisional_classification'],publication_status='STATISTICALLY_FROZEN_EXPLORATORY_ONLY',independent_classification_and_reasons_exact=True)
        output.append(out);checks.append(dict(pert_id=p,reference_code=code,comparison_code=old['provisional_code'],reference_reasons=reason,comparison_reasons=old['downgrade_reasons'],rank=ranks[p],exact_match=True))
    csvwrite(O/'RESULTS/final_blinded_classification.csv',output)
    csvwrite(O/'QA/INDEPENDENT_FINAL_CLASSIFICATION_COMPARISON.csv',checks)
    files=[O/'RESULTS/all_query_statistics.csv',O/'RESULTS/final_blinded_classification.csv',O/'RESULTS/robustness_recalculations.csv',O/'QA/INDEPENDENT_FULL_NULL_VALIDATION.json',O/'QA/INDEPENDENT_FINAL_CLASSIFICATION_COMPARISON.csv']
    manifest=[dict(file=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in files]
    csvwrite(O/'STATISTICAL_RESULTS_FREEZE_SHA256.csv',manifest)
    for r in manifest:assert sha(r['file'])==r['sha256']
    js(O/'STATISTICAL_RESULTS_FREEZE.json',dict(utc=now(),status='BLINDED_STATISTICAL_RESULTS_FROZEN',compounds=1792,annotation_not_performed=True,manifest_sha256=sha(O/'STATISTICAL_RESULTS_FREEZE_SHA256.csv'),independent_classification_and_reason_matches=len(checks)))
    event('BLINDED_STATISTICAL_RESULTS_FROZEN',compounds=1792,name_and_TCM_annotation=False)

if __name__=='__main__':run()
