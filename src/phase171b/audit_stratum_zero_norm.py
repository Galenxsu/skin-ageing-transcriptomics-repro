from common171b import *
import numpy as np

genes=csvread(A/'continuous_signature_definition.csv');labels=sorted({r['null_stratum'] for r in genes});indices=[np.array([i for i,g in enumerate(genes) if g['null_stratum']==s]) for s in labels]
zs=np.empty((107201,len(labels)),dtype=np.int64)
for lo in range(0,107201,256):
    source=np.load(O/'MATRIX/frozen_107201x8227_float64.npy',mmap_mode='r');zero=source[lo:lo+256]==0
    for j,idx in enumerate(indices):zs[lo:lo+len(zero),j]=np.sum(zero[:,idx],axis=1)
    source._mmap.close()
membership={}
for fn in ['discrete_signature_membership_preview.csv','pathway_signature_membership_preview.csv']:
    for r in csvread(A/fn):
        q=r['query_id']
        if q not in membership:membership[q]=np.zeros(len(labels),dtype=np.int64)
        if float(r['weight'])!=0:membership[q][labels.index(r['null_stratum'])]+=1
rows=[]
for q,n in membership.items():
    at_risk=np.all(zs>=n[None,:],axis=1)
    rows.append(dict(query_id=q,profiles_without_combinatorial_zero_norm_exclusion=int(at_risk.sum()),zero_norm_impossible_for_every_frozen_permutation=not bool(at_risk.any()),proof='For every profile, at least one frozen stratum requires more distinct nonzero-weight query members than that profile has zero expression elements in that stratum.'))
csvwrite(O/'QA/STRATIFIED_NULL_ZERO_NORM_PROOF.csv',rows)
print([(r['query_id'],r['profiles_without_combinatorial_zero_norm_exclusion']) for r in rows],flush=True)
