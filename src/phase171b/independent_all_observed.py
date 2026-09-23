from common171b import *
import numpy as np
from independent_full_null import independent_groups,reference_group

rows=[r for r in csvread(A/'discovery_validation_assignment_preview.csv') if r['metadata_eligible']=='True']
genes=csvread(A/'continuous_signature_definition.csv');lookup={r['HGNC_ID']:i for i,r in enumerate(genes)}
queries={'G2_R2_CAP8':np.asarray([float(r['main_R2_L2_weight']) for r in genes]),'G2_R1':np.asarray([float(r['sensitivity_R1_L2_weight']) for r in genes])}
for filename in ['discrete_signature_membership_preview.csv','pathway_signature_membership_preview.csv']:
    for r in csvread(A/filename):
        key=r['query_id']
        if key not in queries:queries[key]=np.zeros(8227,dtype=np.float64)
        queries[key][lookup[r['HGNC_ID']]]=float(r['weight'])
names=jr(O/'RESULTS/query_order.json');assert set(names)==set(queries)
X=np.load(O/'MATRIX/frozen_107201x8227_float64.npy',mmap_mode='r');ref=np.full((len(rows),len(names)),np.nan)
for j,name in enumerate(names):
    w=queries[name];idx=np.flatnonzero(w);v=w[idx];wn=np.sqrt(np.einsum('i,i->',v,v,optimize=False))
    for lo in range(0,len(rows),256):
        block=X[lo:lo+256,idx]
        zn=np.sqrt(np.einsum('ij,ij->i',block,block,optimize=False))
        ref[lo:lo+len(block),j]=np.einsum('ij,j->i',block,v,optimize=False)/(zn*wn)
    print('FULL_INDEPENDENT_OBSERVED',name,flush=True)
X._mmap.close()
ck,pc=independent_groups([(r['pert_id'],r['fold'],r['base_cell_id'],r['cell_id'],float(r['dose_uM']),float(r['time_h'])) for r in rows])
bk,cb=independent_groups([k[:3] for k in ck]);pk,bp=independent_groups([k[:2] for k in bk])
rc=reference_group(ref,*pc);rb=reference_group(rc,*cb);rp=reference_group(rb,*bp)
context_order=[(r['pert_id'],r['fold'],r['base_cell_id'],r['cell_id'],float(r['dose_uM']),float(r['time_h'])) for r in csvread(O/'RESULTS/context_order.csv')]
base_order=[(r['pert_id'],r['fold'],r['base_cell_id']) for r in csvread(O/'RESULTS/base_cell_order.csv')]
compound_order=[(r['pert_id'],r['fold']) for r in csvread(O/'RESULTS/compound_fold_order.csv')]
def align(values,keys,target):
    mapping={k:i for i,k in enumerate(keys)};out=np.full((len(target),len(names)),np.nan)
    for i,k in enumerate(target):
        if k in mapping:out[i]=values[mapping[k]]
    return out
datasets={'profile':ref,'context':align(rc,ck,context_order),'base':align(rb,bk,base_order),'compound':align(rp,pk,compound_order)}
audit=[];files=[]
for level,b in datasets.items():
    primary=O/'RESULTS'/f'observed_{level}_scores.npy';a=np.load(primary)
    assert a.shape==b.shape
    state_ok=np.array_equal(np.isnan(a),np.isnan(b)) and np.array_equal(np.isposinf(a),np.isposinf(b)) and np.array_equal(np.isneginf(a),np.isneginf(b))
    finite=np.isfinite(a)&np.isfinite(b);difference=np.abs(a-b)
    with np.errstate(invalid='ignore',divide='ignore'):relative=np.where(b==0,np.where(difference==0,0,np.inf),difference/np.abs(b))
    sign_ok=np.array_equal(np.sign(a[finite]),np.sign(b[finite]));numeric_ok=bool(np.all(difference[finite]<=1e-12+1e-10*np.abs(b[finite])))
    for label,value in [('reference',b),('absolute_difference',difference),('relative_difference',relative)]:
        path=O/'QA'/f'FULL_INDEPENDENT_{level}_{label}.npy';np.save(path,value);files.append(dict(level=level,role=label,path=str(path),sha256=sha(path)))
    files.append(dict(level=level,role='primary_comparison_value',path=str(primary),sha256=sha(primary)))
    audit.append(dict(level=level,comparisons=int(finite.sum()),missing_states_match=state_ok,directions_match=sign_ok,tolerance_pass=numeric_ok,maximum_absolute_difference=float(np.nanmax(difference))))
csvwrite(O/'QA/FULL_OBSERVED_COMPARISON_ARRAY_INDEX.csv',files)
js(O/'QA/FULL_INDEPENDENT_OBSERVED_VALIDATION.json',dict(utc=now(),pass_all=all(r['missing_states_match'] and r['directions_match'] and r['tolerance_pass'] for r in audit),audit=audit,queries=names))
assert all(r['missing_states_match'] and r['directions_match'] and r['tolerance_pass'] for r in audit),'Independent full observed values disagree; stop final classification'
event('ALL_OBSERVED_PROFILE_CONTEXT_BASE_COMPOUND_INDEPENDENTLY_VERIFIED',comparisons=sum(r['comparisons'] for r in audit))
