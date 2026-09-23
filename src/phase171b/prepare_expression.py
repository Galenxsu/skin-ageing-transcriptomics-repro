from common171b import *
import ctypes, shutil, time
from ctypes import wintypes
import numpy as np
sys.path.insert(0,str(H/'runtime_dependencies'))
import h5py

class Memory(ctypes.Structure):
    _fields_=[('length',wintypes.DWORD),('load',wintypes.DWORD)]+[(x,ctypes.c_ulonglong) for x in ['total','avail','tp','ap','tv','av','ae']]

def run():
    dirs()
    assert np.__version__=='2.3.5'
    for r in csvread(O/'ADDENDUM_SHA256_MANIFEST.csv'):
        assert sha(O/r['file'])==r['sha256'],r['file']
    for r in csvread(A/'FREEZE_MANIFEST.csv'):
        assert sha(A/r['relative_path'])==r['sha256'],r['relative_path']
    receipt=jr(O/'ADDENDUM_FREEZE_VERIFICATION.json')
    assert receipt['status']=='GO_NUMERICAL_TOLERANCE_ADDENDUM_FROZEN'
    assert receipt['expression_values_read'] is False
    m=Memory();m.length=ctypes.sizeof(m)
    assert ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
    resource=dict(utc=now(),available_RAM_GiB=m.avail/2**30,free_disk_GiB=shutil.disk_usage(O).free/2**30)
    js(O/'QA/RESOURCE_GATE_BEFORE_EXPRESSION.json',resource)
    assert resource['available_RAM_GiB']>=12 and resource['free_disk_GiB']>=60,'RESOURCE_GATE_FAILED'
    sourcehash=sha(M)
    assert sourcehash=='1f31891cc7138205688f9d1f2282004517c9e61a14c1c3d64a241b7582976a2e'
    rows=[r for r in csvread(A/'discovery_validation_assignment_preview.csv') if r['metadata_eligible']=='True']
    genes=csvread(A/'continuous_signature_definition.csv')
    assert len(rows)==107201 and len(genes)==8227
    assert len({r['sig_id'] for r in rows})==107201
    assert len({r['pert_id'] for r in rows})==1792
    assert len({r['HGNC_ID'] for r in genes})==8227
    assert [r['HGNC_ID'] for r in genes]==sorted(r['HGNC_ID'] for r in genes)
    ds={r['sig_id'] for r in rows if r['fold']=='DISCOVERY'}
    vs={r['sig_id'] for r in rows if r['fold']=='VALIDATION'}
    dd={x for r in rows if r['fold']=='DISCOVERY' for x in r['distil_id'].split('|')}
    vd={x for r in rows if r['fold']=='VALIDATION' for x in r['distil_id'].split('|')}
    assert not(ds&vs) and not(dd&vd)
    with h5py.File(M,'r') as f:
        decode=lambda x:x.decode() if isinstance(x,bytes) else str(x)
        sig=[decode(x) for x in f['0/META/COL/id'][:]]
        gid=[decode(x) for x in f['0/META/ROW/id'][:]]
        assert len(set(sig))==118050 and len(set(gid))==12328
        si={s:i for i,s in enumerate(sig)}
        gi=np.array([int(r['matrix_index']) for r in genes],dtype=np.int64)
        assert all(gid[i]==r['entrez'] for i,r in zip(gi,genes))
        inds=np.array([si[r['sig_id']] for r in rows],dtype=np.int64)
        data=f['0/DATA/0/matrix']
        assert data.shape==(118050,12328)
        start=dict(utc=now(),freeze_utc=receipt['frozen_utc'],freeze_manifest_sha256=sha(O/'ADDENDUM_SHA256_MANIFEST.csv'),source_sha256=sourcehash,phase='STARTING_FIRST_EXPRESSION_READ',dtype='float64',signature_chunk=256)
        assert start['utc']>receipt['frozen_utc']
        dest=O/'MATRIX/frozen_107201x8227_float64.npy'
        assert not dest.exists(),'Existing matrix requires checkpoint recovery'
        if (O/'EXPRESSION_READ_STARTED.json').exists():
            event('PRE_READ_LOGGING_ERROR_RECOVERY',error='TypeError duplicate utc keyword before matrix allocation/read; event writer corrected; original authorization receipt retained',matrix_did_not_exist=True)
        else: js(O/'EXPRESSION_READ_STARTED.json',start)
        event('FIRST_EXPRESSION_READ_AUTHORIZED',**start)
        X=np.lib.format.open_memmap(dest,mode='w+',dtype=np.float64,shape=(len(rows),len(genes)))
        integrity=[]
        for lo in range(0,len(rows),256):
            ix=inds[lo:lo+256]; order=np.argsort(ix)
            block=np.asarray(data[ix[order],:],dtype=np.float64)[:,gi][np.argsort(order)]
            X[lo:lo+len(ix)]=block
            finite=np.isfinite(block)
            for j,z in enumerate(block):
                integrity.append(dict(sig_id=rows[lo+j]['sig_id'],finite_genes=int(finite[j].sum()),NaN=int(np.isnan(z).sum()),positive_inf=int(np.isposinf(z).sum()),negative_inf=int(np.isneginf(z).sum()),constant=bool(np.all(z==z[0])),zero_norm=bool(np.all(z==0))))
            if lo%10240==0: print('READ_SIGNATURES',lo+len(ix),flush=True)
        X.flush()
        fixed={r['sig_id'] for r in csvread(O/'QA/FROZEN_INDEPENDENT_AUDIT_SAMPLE.csv')}
        checks=[]
        for j,r in enumerate(rows):
            if r['sig_id'] in fixed:
                # Independent scalar-row HDF5 indexing, not the chunked saved matrix.
                single=np.array(data[inds[j],:],dtype=np.float64)[gi]
                ok=np.array_equal(single,X[j],equal_nan=True)
                checks.append(dict(sig_id=r['sig_id'],exact_chunk_vs_single_row=ok))
        assert all(c['exact_chunk_vs_single_row'] for c in checks)
        csvwrite(O/'QA/CHUNKED_VS_SINGLE_ROW_INTEGRITY.csv',checks)
    csvwrite(O/'RESULTS/SIGNATURE_MATRIX_INTEGRITY.csv',integrity)
    csvwrite(O/'MATRIX/signature_order.csv',[dict(row=i,sig_id=r['sig_id'],pert_id=r['pert_id'],fold=r['fold']) for i,r in enumerate(rows)])
    csvwrite(O/'MATRIX/gene_order.csv',[dict(column=i,HGNC_ID=r['HGNC_ID'],entrez=r['entrez'],matrix_index=r['matrix_index']) for i,r in enumerate(genes)])
    summary=dict(utc=now(),status='EXPRESSION_INTEGRITY_COMPLETED_NO_SCORES_YET',shape=[107201,8227],dtype='float64',source_shape=[118050,12328],source_dtype='float32',nonfinite_signatures=sum(r['finite_genes']!=8227 for r in integrity),constant_signatures=sum(r['constant'] for r in integrity),zero_norm_signatures=sum(r['zero_norm'] for r in integrity),single_row_crosschecks=len(checks),gene_order_match=True,signature_linkage_one_to_one=True,shared_sig_ids=0,shared_distil_ids=0,matrix_sha256=sha(dest))
    js(O/'RESULTS/EXPRESSION_INTEGRITY_SUMMARY.json',summary);event('EXPRESSION_INTEGRITY_COMPLETED',**summary)
    print(json.dumps(summary),flush=True)

if __name__=='__main__':run()
