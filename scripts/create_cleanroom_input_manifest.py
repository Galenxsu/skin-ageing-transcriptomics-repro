from pathlib import Path
import os,csv,hashlib,json,datetime

required=['CLEANROOM_PHASE171A','CLEANROOM_PHASE171B','CLEANROOM_GCTX']
if any(not os.environ.get(x) for x in required):raise SystemExit('Missing clean-room environment variables')
A=Path(os.environ['CLEANROOM_PHASE171A']);O=Path(os.environ['CLEANROOM_PHASE171B']);M=Path(os.environ['CLEANROOM_GCTX'])
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
paths=[M,A/'FREEZE_MANIFEST.csv',A/'continuous_signature_definition.csv',A/'discovery_validation_assignment_preview.csv',A/'discrete_signature_membership_preview.csv',A/'pathway_signature_membership_preview.csv',A/'eligible_compound_universe.csv',O/'ADDENDUM_SHA256_MANIFEST.csv']
paths+=sorted((O/'SCRIPTS').glob('*.py'))
rows=[{'path':str(p.resolve()),'bytes':p.stat().st_size,'sha256':sha(p),'role':'CLEANROOM_FROZEN_INPUT'} for p in paths]
with (O/'EXECUTION_INPUT_HASH_BEFORE.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
receipt={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'CLEANROOM_INPUT_MANIFEST_FROZEN','files':len(rows),'manifest_sha256':sha(O/'EXECUTION_INPUT_HASH_BEFORE.csv'),'gctx_sha256':sha(M)}
(O/'CLEANROOM_INPUT_MANIFEST_RECEIPT.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps(receipt,indent=2))
