from common171b import *
from numeric_rules import unit_tests
dirs()
assert not (O/'EXPRESSION_READ_STARTED.json').exists()
frozen=csvread(A/'FREEZE_MANIFEST.csv')
assert len(frozen)==55
checks=[dict(path=str(A/r['relative_path']),expected=r['sha256'],actual=sha(A/r['relative_path'])) for r in frozen]
assert all(r['expected']==r['actual'] for r in checks)
csvwrite(O/'QA/55_FROZEN_FILES_RECHECK.csv',checks)
prior=jr(A/'INPUT_BASELINE.json')
protected=[dict(path=r['path'],expected=r['sha256'],actual=sha(Path(r['path']))) for r in prior]
assert all(r['expected']==r['actual'] for r in protected)
for p in [A/'FREEZE_MANIFEST.csv',A/'FREEZE_VERIFICATION.json']+[Path(r['path']) for r in checks]+list(O.glob('*.*')):
    if not any(r['path']==str(p) for r in protected):protected.append(dict(path=str(p),expected=sha(p),actual=sha(p)))
csvwrite(O/'EXECUTION_INPUT_HASH_BEFORE.csv',protected)
assign=csvread(A/'discovery_validation_assignment_preview.csv');elig=[r for r in assign if r['metadata_eligible']=='True']
uu=csvread(A/'eligible_compound_universe.csv')
order=lambda s:hashlib.sha256(('171B_NUMERIC_AUDIT|'+s).encode()).hexdigest()
fixed=sorted(elig,key=lambda r:(order(r['sig_id']),r['sig_id']))[:64]
pids=set(sorted([r['pert_id'] for r in uu],key=lambda x:(order(x),x))[:8])|{r['pert_id'] for r in elig if r['historically_scored_compound']=='True'}
audit=[dict(sig_id=r['sig_id'],pert_id=r['pert_id'],reason='FIRST64_HASHED_SIGNATURE' if r in fixed else 'FULL_FROZEN_AUDIT_COMPOUND') for r in elig if r in fixed or r['pert_id'] in pids]
csvwrite(O/'QA/FROZEN_INDEPENDENT_AUDIT_SAMPLE.csv',audit)
js(O/'QA/SYNTHETIC_COMPARISON_TESTS.json',unit_tests())
def book(name,rows):return dict(name=name,sheets=[dict(name='Rules',columns=['item','value_or_rule'],rows=rows)])
exact=['signature members','gene members and directions','pert_id/sig_id/distil_id','discovery/heldout allocation','profile/context/compound counts','inclusion/exclusion','null seed','valid null iterations','null extreme counts','test-family membership','classification','all downgrade reasons','missing-value states','significance decisions','reversal/mimicry decisions']
rules=[['dtype','float64'],['atol',1e-12],['rtol',1e-10],['comparison','abs(x1-x2) <= 1e-12 + 1e-10*abs(x2); path2 reference'],['record','comparison value; reference value; absolute difference; relative difference'],['continuous scope','profile/context/compound scores; same frozen aggregate; empirical P with identical integer counts; BH floating representation'],['zero reference','relative difference=0 if difference0, otherwise Inf; absolute comparison formula unchanged'],['decision boundary','DECISION_BOUNDARY_DISCORDANCE: opposite sides of raw q<0.05 -> stop final classification'],['direction','DIRECTION_DISCORDANCE: any negative/zero/positive mismatch -> stop final classification'],['missing states','NA,NaN,Inf,-Inf checked separately; exact state match only'],['rounding','No prior rounding/truncation/formatting'],['no scientific change','Only numerical implementation verification; no score/null/FDR/classification rule change']]+[['EXACT '+str(i+1),v] for i,v in enumerate(exact)]
plan=[['path1','NumPy float64 vector scoring; grouped exact medians'],['path2','Independent direct sum-of-products and sorted medians from same frozen inputs; not final scores from path1'],['fixed audit sample',len(audit)],['audit selection','SHA256 171B_NUMERIC_AUDIT|ID: first64 sig_id, first8 pert_id plus all historical6; all members of audit compounds'],['null','Same frozen index/seed; exact effective iteration and extreme-count comparison'],['critical coverage','All final A/B/C/D; main D/V q in[0.04,0.06]; any discordance; historical6; interval only selects QA scope'],['P and BH','Independent implementation for all tests; identical integer counts and denominators required'],['classification','Independent decision implementation for all1792; exact class/reasons/direction/significance'],['boundary tests','10 synthetic unit cases passed before real expression read'],['failure','Stop final classification/publication; preserve evidence; do not adjust seed,B,threshold,signature']]
gate=[['Phase171A frozen hashes','55/55 matched'],['Historical protected hashes',str(len(prior))+' unchanged'],['Matrix expression read','NO at addendum preparation'],['Actual scores/null/ranking','NONE at addendum preparation'],['Numerical comparison synthetic tests','PASS'],['Independent paths','Direct independent implementation plan executable; actual-data agreement remains to be checked'],['Signature/null/family changes','NONE'],['Final startup condition','After these3XLSX+MD readback/hash receipt only: GO_NUMERICAL_TOLERANCE_ADDENDUM_FROZEN'],['Prior stop records','Preserved, not overwritten']]
js(O/'QA/ADDENDUM_WORKBOOK_PAYLOAD.json',[book('NUMERICAL_COMPARISON_RULE_FREEZE',rules),book('INDEPENDENT_IMPLEMENTATION_VALIDATION_PLAN',plan),book('PHASE171B_AMENDED_GO_NO_GO_CHECKLIST',gate)])
event('ADDENDUM_PREPARATION_COMPLETE',expression_read=False,scores_computed=False,null_generated=False,upstream55_match=True)
print('ADDENDUM_PAYLOAD_READY',len(audit),'audit signatures; no expression read')
