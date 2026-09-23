from common171b import *
import numpy as np, math
from collections import defaultdict

def run():
    frozen=jr(O/'STATISTICAL_RESULTS_FREEZE.json');assert frozen['status']=='BLINDED_STATISTICAL_RESULTS_FROZEN'
    for r in csvread(O/'STATISTICAL_RESULTS_FREEZE_SHA256.csv'):assert sha(r['file'])==r['sha256']
    sources=csvread(O/'POST_STATISTICAL_ANNOTATION_SOURCE_REGISTRY.csv')
    for r in sources:assert sha(r['path'])==r['sha256']
    event('POST_STATISTICAL_NAME_AND_STRUCTURE_ANNOTATION_STARTED',statistical_freeze_utc=frozen['utc'])
    universe=csvread(A/'eligible_compound_universe.csv');ids={r['pert_id']:r for r in universe};assert len(ids)==1792
    final=csvread(O/'RESULTS/final_blinded_classification.csv');final.sort(key=lambda r:int(r['discovery_rank']))
    for r in final:r['reported_name']=ids[r['pert_id']]['reported_name']
    stats=csvread(O/'RESULTS/all_query_statistics.csv');stat={(r['query_id'],r['fold'],r['pert_id']):r for r in stats}
    def named(rows):return [dict(r,reported_name=ids[r['pert_id']]['reported_name']) for r in rows]
    def save(name,rows):
        assert rows,'Empty export '+name
        dest=O/name;csvwrite(dest,rows)
        again=csvread(dest);assert len(again)==len(rows)
        for a,b in zip(rows,again):
            for key in a:assert str(a[key])==b[key],(name,key)
        exports[name]=dict(path=str(dest),rows=len(rows),columns=len(rows[0]),sha256=sha(dest),readback='ALL_CELLS_MATCH')
        return rows
    exports={}
    save('final_compound_classification.csv',final);save('complete_primary_ranking.csv',final)
    discovery=named([stat['G2_R2_CAP8','DISCOVERY',r['pert_id']] for r in final]);heldout=named([stat['G2_R2_CAP8','VALIDATION',r['pert_id']] for r in final])
    save('discovery_results.csv',discovery);save('heldout_results.csv',heldout)
    names=jr(O/'RESULTS/query_order.json');scores=np.load(O/'RESULTS/observed_compound_scores.npy');pids=sorted(ids)
    detailed=[]
    for qi,q in enumerate(names):
        for fi,fold in enumerate(['DISCOVERY','VALIDATION']):
            for pi,p in enumerate(pids):
                key=q,fold,p
                if key in stat:r=dict(stat[key],inference_scope='FROZEN_TEST_FAMILY')
                else:
                    s=scores[fi*1792+pi,qi]
                    r=dict(query_id=q,fold=fold,pert_id=p,score=float(s) if np.isfinite(s) else 'NA',B='NOT_APPLICABLE',seed='NOT_APPLICABLE',null_valid_iterations='NOT_RUN',lower_extreme_count='NOT_RUN',upper_extreme_count='NOT_RUN',raw_P_lower='NOT_RUN',raw_P_upper='NOT_RUN',null_evaluable='NOT_APPLICABLE',BH_FDR_lower='NOT_RUN',family_lower='NOT_APPLICABLE',BH_FDR_upper='NOT_RUN',family_upper='NOT_APPLICABLE',inference_scope='DIRECTIONAL_DESCRIPTIVE_ONLY_NO_NEW_TEST')
                r['reported_name']=ids[p]['reported_name'];detailed.append(r)
    save('RESULTS/all_compound_query_fold_results.csv',detailed)
    save('discrete_grid_results.csv',[r for r in detailed if r['query_id'].startswith('G2_E')])
    save('pathway_results.csv',[r for r in detailed if not r['query_id'].startswith('G2_')])
    save('RESULTS/continuous_compound_results.csv',[r for r in detailed if r['query_id'] in ['G2_R2_CAP8','G2_R1']])
    # Full-data scores are recomputed as medians across all base-cell medians,
    # not as the median of the discovery and held-out compound scores.
    bases=csvread(O/'RESULTS/base_cell_order.csv');bs=np.load(O/'RESULTS/observed_base_scores.npy');by=defaultdict(list)
    for i,r in enumerate(bases):by[r['pert_id']].append(i)
    full=[]
    for p in pids:
        v=bs[by[p],0];v=v[np.isfinite(v)]
        full.append(dict(pert_id=p,reported_name=ids[p]['reported_name'],full_data_main_score=float(np.median(v)) if len(v) else 'NA',base_cell_count=len(v),scope='FULL_DATA_DESCRIPTIVE_ONLY',P='NOT_RUN',FDR='NOT_RUN'))
    full.sort(key=lambda r:(r['full_data_main_score'] if r['full_data_main_score']!='NA' else math.inf,r['pert_id']))
    for i,r in enumerate(full):r['descriptive_rank']=i+1
    save('RESULTS/full_data_descriptive_results.csv',full)
    meta=[r for r in csvread(A/'discovery_validation_assignment_preview.csv') if r['metadata_eligible']=='True']
    save('RESULTS/signature_metadata_linkage.csv',meta)
    for level,orderfile,arrayfile in [('profile',None,'observed_profile_scores.npy'),('context','context_order.csv','observed_context_scores.npy')]:
        arr=np.load(O/'RESULTS'/arrayfile,mmap_mode='r');order=meta if level=='profile' else csvread(O/'RESULTS'/orderfile)
        records=[]
        for i,row in enumerate(order):
            r=dict(row)
            for j,q in enumerate(names):r[q]=float(arr[i,j]) if np.isfinite(arr[i,j]) else 'NA'
            records.append(r)
        save('RESULTS/'+level+'_all_query_scores.csv',records)
    # Historical exact identity mapping is never used to upgrade a class.
    hist=csvread(H/'RESULTS/CHEMICAL_IDENTITY_AND_ALIAS_AUDIT.csv');histby={r['pert_id']:r for r in hist};assert len(histby)==6
    historical=[]
    for r in final:
        if r['pert_id'] in histby:
            h=histby[r['pert_id']];historical.append(dict(r,historical_compound=h['compound'],historical_exact_InChIKey=h['target_inchikey'],comparison_scope='HISTORICALLY_EXAMINED_NO_INDEPENDENT_VALIDATION_CLAIM'))
    assert len(historical)==6;save('historical_six_crosswalk.csv',historical)
    chem=csvread(Path(sources[0]['path']));keymap=defaultdict(list)
    for r in chem:
        if r['inchikey'] not in ('','NA','UNKNOWN'):keymap[r['inchikey']].append(r)
    tcm=[]
    for p in pids:
        u=ids[p];matched=keymap.get(u['inchi_key'],[])
        tcm.append(dict(pert_id=p,reported_name=u['reported_name'],full_inchikey=u['inchi_key'],exact_full_key_match=bool(matched),unique_compound_ids='|'.join(sorted({r['unique_compound_id'] for r in matched})) or 'NOT_MATCHED',source_herb_chinese='|'.join(sorted({r['source_herb_chinese'] for r in matched})) or 'UNKNOWN',source_herb_latin='|'.join(sorted({r['source_herb_latin'] for r in matched})) or 'UNKNOWN',source_evidence_level='|'.join(sorted({r['evidence_level'] for r in matched})) or 'UNKNOWN',primary_reference='|'.join(sorted({r['primary_reference'] for r in matched})) or 'UNKNOWN',mapping_rule='FULL_INCHIKEY_EXACT_ONLY_NO_NAME_OR_FIRST_BLOCK_SUBSTITUTION',mapping_time='AFTER_BLINDED_STATISTICAL_FREEZE',source_path=sources[0]['path'],source_sha256=sources[0]['sha256'],classification_modified=False))
    save('tcm_mapping_after_statistical_freeze.csv',tcm)
    identity=[];tox=[]
    for r in final:
        p=r['pert_id'];u=ids[p]
        identity.append(dict(pert_id=p,reported_name=u['reported_name'],reported_aliases=u['reported_aliases'],inchi_key=u['inchi_key'],canonical_smiles=u['canonical_smiles'],structure_unique=u['structure_unique'],shared_full_key_pert_ids=u['shared_full_key_pert_ids'],salt_status=u['salt_status'],stereoisomer_policy=u['stereoisomer_policy'],parent_form_status=u['parent_form_status'],mechanism=u['mechanism'],normal_adult_skin_exact_contexts=0,external_skin_evidence='NOT_ASSESSED_IN_THIS_FROZEN_COMPUTATIONAL_STAGE',source=str(A/'eligible_compound_universe.csv')))
        tox.append(dict(pert_id=p,reported_name=u['reported_name'],dose_robustness_limit=r['dose_dependence'],functional_toxicity_status='FUNCTIONAL_TOXICITY_NOT_ASSESSABLE',toxicity_inference='TOXICITY_INTERPRETATION_NOT_PERMITTED',nonspecific_response='CANNOT_EXCLUDE_NOT_DEMONSTRATED',caution='Dose robustness failure includes insufficient low-dose support; it is not direct evidence of toxicity or demonstrated high-dose-only activity.'))
    save('RESULTS/compound_identity_resolution.csv',identity);save('RESULTS/toxicity_and_nonspecific_audit.csv',tox)
    old=csvread(H/'RESULTS/PRIMARY_36_TEST_RESULTS.csv');comparison=[]
    for r in old:
        matches=[h for h in hist if h['compound_id']==r['compound_id']];assert len(matches)==1
        p=matches[0]['pert_id'];current=next(x for x in final if x['pert_id']==p)
        comparison.append(dict(phase170B_test_id=r['test_id'],pert_id=p,compound=r['compound'],phase170B_query=r['query'],phase170B_layer=r['layer'],phase170B_original_score=r['score'],phase170B_cosine_orientation_score=-float(r['score']),phase170B_reversal_FDR=r['q_upper'],phase171B_main_discovery_score=current['discovery_score'],phase171B_main_discovery_FDR=current['discovery_FDR'],phase171B_heldout_FDR=current['heldout_FDR'],phase171B_class=current['final_classification'],comparability='DIFFERENT_QUERY_UNIVERSE_NULL_AND_FAMILY_DESCRIPTIVE_ONLY',history_preserved=True))
    save('RESULTS/phase170B_vs_phase171B_comparison.csv',comparison)
    robust=csvread(O/'RESULTS/robustness_recalculations.csv')
    for name,analyses in [('leave_one_cell_line_out',['LEAVE_ONE_BASE_CELL_OUT']),('leave_one_context_out',['LEAVE_ONE_CONTEXT_OUT']),('dose_and_time_dependence',['LEAVE_ONE_DOSE_OUT','LEAVE_ONE_TIME_OUT','TIME_STRATUM','REMOVE_GE10_UM']),('HQ_sensitivity',['FULL_PROFILE','HQ'])]:save('RESULTS/'+name+'.csv',[r for r in robust if r['analysis'] in analyses])
    save('RESULTS/signature_method_stability.csv',[{k:r[k] for k in ['pert_id','reported_name','signature_robustness','downgrade_reasons','final_classification']} for r in final])
    logs=[json.loads(line) for line in (O/'NULL/null_execution_log.jsonl').read_text(encoding='utf-8').splitlines() if line.strip()]
    keys=sorted({k for r in logs for k in r});flat=[{k:json.dumps(r[k],sort_keys=True) if isinstance(r.get(k),(dict,list)) else r.get(k,'NA') for k in keys} for r in logs]
    save('null_execution_log.csv',flat)
    boundaries=[dict(topic='Maximum claim',allowed='Exploratory perturbational reversal candidate in the frozen GSE70138 limited universe with internal held-out support only if class A.',forbidden='Validated anti-skin-ageing efficacy; whole LINCS best compound; clinical candidate'),dict(topic='Normal skin',allowed='Normal adult skin exact matched contexts = 0.',forbidden='Effective in normal adult skin'),dict(topic='Historical six',allowed='Historically examined; internal holdout is not independent validation for these compounds.',forbidden='New independent replication'),dict(topic='Toxicity',allowed='Functional toxicity is not assessable from these expression results.',forbidden='Toxic or safe based on reversal or dose alone'),dict(topic='Null results',allowed='No evidence meeting this frozen test and support rule.',forbidden='Compound is ineffective'),dict(topic='Phase170B',allowed='Historical conclusions and numerical outputs retained unchanged.',forbidden='Retroactive replacement of Phase170B conclusions')]
    save('RESULTS/claim_boundary_matrix.csv',boundaries)
    retired=[]
    for file in ['FINAL_REPORT_ZH.md','SUMMARY.json']:
        p=O/file
        if p.exists():retired.append(dict(path=str(p),sha256=sha(p),status='RETAINED_NOT_DELETED_HISTORICAL_PRE_ADDENDUM_NO_GO',current_pointer='FINAL_REPORT_EXECUTED_ZH.md'))
    save('RESULTS/retired_not_deleted_register.csv',retired)
    js(O/'QA/COMPLETE_CSV_EXPORT_INDEX.json',exports)
    event('COMPLETE_MACHINE_READABLE_EXPORTS_WRITTEN',csv_files=len(exports),all_cells_readback=True,compound_count=1792,statistical_classes_unchanged=True)

if __name__=='__main__':run()
