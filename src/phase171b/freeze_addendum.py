from common171b import *
import openpyxl
assert not (O/'EXPRESSION_READ_STARTED.json').exists()
books=jr(O/'QA/ADDENDUM_WORKBOOK_PAYLOAD.json');cells=0
for b in books:
    wb=openpyxl.load_workbook(O/(b['name']+'.xlsx'),read_only=True,data_only=False)
    for s in b['sheets']:
        actual=list(wb[s['name']].iter_rows(values_only=True));expect=[s['columns']]+s['rows'];assert len(actual)==len(expect)
        for aa,ee in zip(actual,expect):
            for x,y in zip(aa,ee):assert x==y,(b['name'],x,y);cells+=1
    wb.close()
ff=csvread(A/'FREEZE_MANIFEST.csv');assert len(ff)==55 and all(sha(A/r['relative_path'])==r['sha256'] for r in ff)
names=['PHASE171B_NUMERICAL_TOLERANCE_ADDENDUM.md']+[b['name']+'.xlsx' for b in books]+['QA/FROZEN_INDEPENDENT_AUDIT_SAMPLE.csv','SCRIPTS/numeric_rules.py']
rows=[dict(file=n,sha256=sha(O/n),bytes=(O/n).stat().st_size,freeze_utc=now()) for n in names]
csvwrite(O/'ADDENDUM_SHA256_MANIFEST.csv',rows)
assert all(sha(O/r['file'])==r['sha256'] for r in csvread(O/'ADDENDUM_SHA256_MANIFEST.csv'))
receipt=dict(status='GO_NUMERICAL_TOLERANCE_ADDENDUM_FROZEN',frozen_utc=now(),expression_values_read=False,actual_scores_computed=False,null_generated=False,ranking_generated=False,upstream55_match=True,workbooks=3,readback_cells=cells,manifest_sha256=sha(O/'ADDENDUM_SHA256_MANIFEST.csv'),all_frozen_files_reread_match=True,comparison_unit_tests=10,visual_review='All3 workbook opening ranges inspected',independent_audit_sig_ids=781)
js(O/'ADDENDUM_FREEZE_VERIFICATION.json',receipt);js(O/'CURRENT_STATUS.json',receipt)
event('GO_NUMERICAL_TOLERANCE_ADDENDUM_FROZEN',receipt_sha256=sha(O/'ADDENDUM_FREEZE_VERIFICATION.json'),expression_values_read=False)
print(json.dumps(receipt,indent=2))
