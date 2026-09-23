from common import *
import pandas as pd,numpy as np,collections,math
assert sha(O/'RULES_LOCK.md')==json.loads((O/'RULES_LOCK_RECORD.json').read_text(encoding='utf-8'))['sha256']
branches=['P2_M0','P2_M1','P2_M2','P2_M3','P1_M0','P3_M0']
frames={}
for b in branches+['P0_M0']:
 f=pd.DataFrame(read(phase(163)/f'RESULTS/GENE_{b}.csv')).set_index('id')
 for k in ['logFC','t','P.Value','adj.P.Val']:f[k]=pd.to_numeric(f[k],errors='coerce')
 assert f.index.is_unique
 f['pct_rank']=f.t.abs().rank(method='average',pct=True);frames[b]=f
main=frames['P2_M0'];assert len(main)==11412
tissue=pd.DataFrame(read(phase(135)/'csv/GSE18876_ALL_GENE_RESULTS.csv')).set_index('HGNC_ID')
tissue['rank_pct']=pd.to_numeric(tissue.moderated_t).abs().rank(pct=True,method='average')
bm=read(phase(164)/'RESULTS/226189_STABLE_ID_MAPPING.csv')
bids={r['ensembl']:r['HGNC'] for r in bm if r['status']=='RESOLVED_UNIQUE'}
br=read(R/'22_phase2_5_results/GSE226189_FINAL_SIGNATURE_STABILITY.csv');blookup={bids[r['gene_id']]:r for r in br if r['gene_id'] in bids}
sc=read(phase(137)/'csv/LEADING_EDGE_CELLTYPE_EFFECTS.csv');scdon=read(phase(137)/'csv/LEADING_EDGE_DONOR_EXPRESSION.csv')
donors=collections.defaultdict(lambda:collections.defaultdict(set))
for r in scdon:
 if r['primary_coverage']=='True':donors[(r['HGNC_ID'],r['author_celltype'])][r['age_group']].add(r['donor_id'])
sce=collections.defaultdict(list)
for r in sc:
 if r['lineage']=='Keratinocyte':sce[r['HGNC_ID']].append(r)
out=[];scdetail=[]
for gid,r in main.iterrows():
 sign=np.sign(r.logFC);ranks=[];dirs=[];rec=dict(gene_id=gid,symbol=r.gene_symbol,main_logFC=float(r.logFC),main_t=float(r.t),main_P=float(r['P.Value']),main_FDR=float(r['adj.P.Val']),annotation_status='FROZEN_PHASE131_UNIQUE_REPRESENTATIVE',result_version='TWO_AXIS_RECALCULATED_RESULT')
 for b,f in frames.items():
  present=gid in f.index
  rec[b+'_logFC']=float(f.loc[gid,'logFC']) if present else None
  rec[b+'_FDR']=float(f.loc[gid,'adj.P.Val']) if present else None
  rec[b+'_rank_percentile']=float(f.loc[gid,'pct_rank']) if present else None
  if b in branches:
   dirs.append(present and np.sign(f.loc[gid,'logFC'])==sign)
   if present:ranks.append(float(f.loc[gid,'pct_rank']))
 span=max(ranks)-min(ranks);robust=all(dirs) and len(ranks)==6 and span<=.20
 rec.update(internal_all_six_evaluable=len(ranks)==6,internal_all_six_same_direction=all(dirs),rank_percentile_span=span,internal_robust=bool(robust))
 tr=tissue.loc[gid] if gid in tissue.index else None
 tsame=tr is not None and np.sign(num(tr.age_coefficient))==sign
 rec.update(tissue_same_direction=bool(tsame),tissue_beta=num(tr.age_coefficient) if tr is not None else None,tissue_signed_t=num(tr.moderated_t) if tr is not None else None,tissue_P=num(tr.P_value) if tr is not None else None,tissue_FDR=num(tr.BH_FDR) if tr is not None else None,tissue_rank_percentile=float(tr.rank_pct) if tr is not None else None,tissue_FDR44_member=bool(tr is not None and num(tr.BH_FDR)<.05),tissue_gene_LOO_same_rate='NOT_AVAILABLE_IN_FROZEN_OUTPUT',tissue_strong_support=False,tissue_standardized_effect='NOT_RECOVERED_FROM_FROZEN_TABLE')
 strongsc=[];local=[]
 for sr in sce.get(gid,[]):
  ds=donors[(gid,sr['author_celltype'])];ny=len(ds['YOUNG']);no=len(ds['OLD'])
  eligible=sr['eligible_type']=='True' and num(sr['n_detected_donors'])>=4 and ny>=2 and no>=2 and num(sr['LOO_sign_changes'])==0 and np.sign(num(sr['old_minus_young']))==sign
  local.append(sr['author_celltype']);
  if eligible:strongsc.append(sr['author_celltype'])
  scdetail.append(dict(gene_id=gid,symbol=r.gene_symbol,celltype=sr['author_celltype'],pathway=sr['pathway'],donors_young=ny,donors_old=no,n_detected=sr['n_detected_donors'],old_minus_young=sr['old_minus_young'],LOO_sign_changes=sr['LOO_sign_changes'],support=bool(eligible),inference='DESCRIPTIVE_AGE_BATCH_CONFOUNDED'))
 rec.update(single_cell_keratinocyte_coverage=';'.join(sorted(set(local))) or 'NOT_TESTED_IN_FROZEN_LE_SUBSET',single_cell_direction_support=bool(strongsc),single_cell_support_types=';'.join(sorted(set(strongsc))),single_cell_boundary='5 shared donors; age/batch confounded; only frozen leading-edge subset')
 signal=r['adj.P.Val']<.05
 tier='EPIDERMAL_TIER_A' if signal and robust and strongsc else 'EPIDERMAL_TIER_B' if signal and robust else 'EPIDERMAL_TIER_C' if signal else 'NOT_DISCOVERY_SIGNAL'
 brw=blookup.get(gid);bsame=brw is not None and np.sign(num(brw['beta_per_10_years']))==sign
 rec.update(new_axis_A_tier=tier,new_axis_B_label='NOT_TESTED' if brw is None else 'FIBROBLAST_CONTEXT_SUPPORT' if bsame else 'CONTEXT_DISCORDANT',axis_B_beta=num(brw['beta_per_10_years']) if brw else None,axis_B_FDR=num(brw['FDR']) if brw else None,axis_B_stability=brw['stability_class'] if brw else 'NOT_TESTED',method_sensitive_downgrade=bool(signal and not robust),axis_B_used_in_A_rule=False)
 out.append(rec)
save('AXIS_A_GENE_EVIDENCE_MATRIX',out);save('SC_GENE_SUPPORT_DETAIL',scdetail)
stable=[]
for r in br:
 if r['stability_class'] in ['STABLE_A','STABLE_B']:
  stable.append(dict(gene_id=bids.get(r['gene_id'],'UNRESOLVED'),ensembl=r['gene_id'],symbol=r['gene_symbol'],historical_class=r['stability_class'],axis_B_label='FIBROBLAST_CONTEXT_'+r['stability_class'],beta=r['beta_per_10_years'],FDR=r['FDR'],internal_only=True,query_status='FIBROBLAST_QUERY_NOT_SUPPORTED',reason='No FDR genes; internal stability is not sufficient disease-signature evidence under user gate'))
assert len(stable)==73;save('AXIS_B_FIBROBLAST_CONTEXT_MATRIX',stable)
# Historical signatures kept at record level and unique-symbol level, without reannotation rescue.
by_symbol={r['symbol']:r for r in out};hist=[]
for p in sorted((R/'26_phase3_final_signatures').glob('*FULL_TABLE.csv')):
 for r in read(p):
  sym=r.get('gene_symbol');nr=by_symbol.get(sym);old=r.get('subclass') or r.get('stability_class') or r.get('evidence_tier') or 'UNSPECIFIED'
  hist.append(dict(symbol=sym,historical_file=p.name,historical_tier=old,new_gene_id=nr['gene_id'] if nr else 'UNMAPPED_TO_AXIS_A',new_axis_A_tier=nr['new_axis_A_tier'] if nr else 'NOT_ELIGIBLE',new_axis_B_label=nr['new_axis_B_label'] if nr else 'NOT_ASSESSABLE',historical_status='HISTORICAL_PRE_ROLE_RECLASSIFICATION_RESULT',change_reason='New axis criteria and Phase163 version; no sole causal attribution to removing fibroblast gate'))
save('HISTORICAL_TO_NEW_GENE_TIER_CROSSWALK',hist)
oldD={r['gene_symbol'] for r in read(R/'26_phase3_final_signatures/D_CROSS_COHORT_FULL_TABLE.csv')}
restored=[dict(gene_id=r['gene_id'],symbol=r['symbol'],new_tier=r['new_axis_A_tier'],historical_D_member=r['symbol'] in oldD,rule_counterfactual_discordance=r['new_axis_B_label']=='CONTEXT_DISCORDANT',causal_attribution='MULTIPLE_CHANGES_NOT_SOLELY_GATE_REMOVAL') for r in out if r['new_axis_A_tier'] in ['EPIDERMAL_TIER_A','EPIDERMAL_TIER_B'] and r['new_axis_B_label']=='CONTEXT_DISCORDANT']
save('AXIS_A_RETAINED_DESPITE_FIBROBLAST_DISCORDANCE',restored)
genes=read(R/'90_pre_lincs_compound_expansion/2026-08-15_stage03A_amendment01_full_gene_space/03_full_12328_gene_universe.csv')
symbol_map=collections.defaultdict(list)
for g in genes:symbol_map[g['gene_symbol']].append(g)
hp=track(R/'130_SMALL_DATA_ANNOTATION_PREFLIGHT_2026-09-09_v2/annotations/hgnc_complete_set_2026-09-09.txt')
with hp.open(encoding='utf-8-sig',newline='') as fh:
 hgnc={r['hgnc_id']:r for r in csv.DictReader(fh,delimiter='\t') if r['status']=='Approved'}
entrez_map=collections.defaultdict(list)
for g in genes:entrez_map[g['entrez_id']].append(g)
def stable_maps(r):
 return entrez_map.get(hgnc.get(r['gene_id'],{}).get('entrez_id',''),[])
preds={'Q_A1':lambda r:r['new_axis_A_tier']=='EPIDERMAL_TIER_A','Q_A2':lambda r:r['new_axis_A_tier'] in ['EPIDERMAL_TIER_A','EPIDERMAL_TIER_B'],'Q_A3':lambda r:r['main_FDR']<.05 and r['internal_robust'] and r['tissue_same_direction'],'Q_A4':lambda r:r['new_axis_A_tier'] in ['EPIDERMAL_TIER_A','EPIDERMAL_TIER_B'] and r['single_cell_direction_support']}
queries=[];qgenes=[];full=[]
for q,pred in preds.items():
 candidate=[r for r in out if pred(r)];selected={}
 for sign,label in [(1,'UP'),(-1,'DOWN')]:
  rows=sorted([r for r in candidate if np.sign(r['main_logFC'])==sign],key=lambda r:(-abs(r['main_t']),r['gene_id']));selected[label]=rows[:150]
  for i,r in enumerate(rows):
   maps=stable_maps(r);ok=len(maps)==1
   full.append(dict(query_id=q,direction=label,rank=i+1,gene_id=r['gene_id'],symbol=r['symbol'],selected_before_mapping=i<150,mapping_status='UNIQUE' if ok else 'UNMAPPED_OR_AMBIGUOUS'))
   if i<150 and ok:qgenes.append(dict(query_id=q,direction=label,gene_id=r['gene_id'],symbol=r['symbol'],entrez=maps[0]['entrez_id'],matrix_index=int(maps[0]['matrix_row_index_zero_based']),weight_sign=sign))
 up=[x for x in qgenes if x['query_id']==q and x['direction']=='UP'];dn=[x for x in qgenes if x['query_id']==q and x['direction']=='DOWN']
 nu,nd=len(up),len(dn);ru=nu/len(selected['UP']) if selected['UP'] else 0;rd=nd/len(selected['DOWN']) if selected['DOWN'] else 0
 valid=min(nu,nd)>=10 and max(nu,nd)/max(1,min(nu,nd))<=3 and min(ru,rd)>=.8
 queries.append(dict(query_id=q,uncapped_candidates=len(candidate),selected_up=len(selected['UP']),selected_down=len(selected['DOWN']),mapped_up=nu,mapped_down=nd,mapping_rate_up=ru,mapping_rate_down=rd,status='VALID' if valid else 'QUERY_TOO_SMALL_OR_UNBALANCED' if min(nu,nd)<10 or max(nu,nd)/max(1,min(nu,nd))>3 else 'INSUFFICIENT_MAPPING_COVERAGE',role='PRIMARY' if q=='Q_A1' else 'PRESPECIFIED_SENSITIVITY',no_post_result_rescue=True))
save('NEW_QUERY_SIGNATURES',queries);save('QUERY_GENES_MAPPED',qgenes);save('QUERY_GENES_UNCAPPED',full)
dump('GENE_QUERY_SUMMARY',dict(tiers=dict(collections.Counter(r['new_axis_A_tier'] for r in out)),n_all=len(out),n_signal=sum(r['main_FDR']<.05 for r in out),internal_robust=sum(r['internal_robust'] for r in out),historical_records=len(hist),historical_unique_symbols=len({r['symbol'] for r in hist}),historical_mapped_symbols=len({r['symbol'] for r in hist if r['new_gene_id']!='UNMAPPED_TO_AXIS_A'}),context_discordant_retained=len(restored),queries=queries,LOO_gene_gap='Frozen Phase135 summary and RDS contain pathway LOO only; gene LOO proportions unavailable; no strong tissue support imputed'))
print((D/'GENE_QUERY_SUMMARY.json').read_text(encoding='utf-8'))
